from __future__ import annotations

from collections import deque
import json
import os
import signal
import subprocess
import threading
import time
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path

from . import storage, execution
from .wolfram import parse_result_marker

# Safety bounds. The full step log always lands in runs/<run_id>.log; the
# in-memory event buffer only feeds SSE replay, so capping it bounds server
# memory for chatty binaries without losing information. Finished runs are
# evicted from the engine dict after MAX_FINISHED_RUNS (their runs/*.json
# records stay on disk).
MAX_BUFFERED_EVENTS = 20000
MAX_FINISHED_RUNS = 500

_TERMINAL_STATUSES = ("done", "failed", "cancelled")

# Process-group control is POSIX-only. On Windows the steps still run, but
# cancellation/timeout fall back to killing the direct child only (no group).
POSIX_SESSIONS = hasattr(os, "setsid") and hasattr(os, "killpg")


def kill_process_tree(proc: subprocess.Popen) -> None:
    """Kill a step process and (on POSIX) everything in its group."""
    if POSIX_SESSIONS:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
            return
        except (ProcessLookupError, PermissionError, OSError):
            return
    try:
        proc.kill()
    except OSError:
        pass


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _set_status(run: "Run", status: str) -> bool:
    """Thread-safe status transition.

    A terminal status (done/failed/cancelled) is final: whichever of the
    worker thread or cancel() gets there first wins, the loser skips, and the
    UI sees exactly one terminal status+end pair instead of conflicting ones.
    Returns True when this call performed the transition.
    """
    with run.cond:
        if run.status in _TERMINAL_STATUSES or run.status == status:
            return False
        run.status = status
        run.cond.notify_all()
        return True


class Run:
    def __init__(self, run_id: str, project_id: str, flow_id, label: str, steps: list):
        self.run_id = run_id
        self.project_id = project_id
        self.flow_id = flow_id
        self.label = label
        self.status = "queued"
        self.created_at = _now()
        self.steps = steps
        self.current_step = None
        self.events = deque(maxlen=MAX_BUFFERED_EVENTS)
        self.next_event_id = 0
        self.cancel_requested = False
        self.worker_active = True
        self.cond = threading.Condition()
        self.proc: subprocess.Popen | None = None
        self.on_step_done = None
        self.on_done = None

    def emit(self, event: str, data: dict) -> None:
        with self.cond:
            self.events.append({"id": self.next_event_id, "event": event, "data": data})
            self.next_event_id += 1
            self.cond.notify_all()

    def snapshot(self) -> dict:
        return {
            "run_id": self.run_id,
            "project_id": self.project_id,
            "flow_id": self.flow_id,
            "label": self.label,
            "status": self.status,
            "created_at": self.created_at,
            "current_step": self.current_step,
            "steps": [
                {
                    "id": s["id"],
                    "label": s["label"],
                    "kind": s["kind"],
                    "command": s["command"],
                    "status": s.get("status", "pending"),
                    "returncode": s.get("returncode"),
                }
                for s in self.steps
            ],
        }


class Engine:
    def __init__(self):
        self.runs: dict[str, Run] = {}
        self.lock = threading.Lock()

    def create_run(self, project_id: str, flow_id, label: str, steps: list,
                   on_step_done=None, on_done=None) -> Run:
        # Two concurrent runs writing the same output file would race on the
        # file itself AND poison the .sig fingerprint cache (the second run
        # could bless a half-written file as fresh). Refuse the overlap.
        with self.lock:
            new_outputs = set()
            for step in steps:
                new_outputs.update(step.get("outputs") or [])
            for other in self.runs.values():
                if other.project_id != project_id or not other.worker_active:
                    continue
                clash = new_outputs & {o for s in other.steps for o in (s.get("outputs") or [])}
                if clash or other.worker_active:
                    raise ValueError(
                        "another run is already writing "
                        + (", ".join(sorted(clash)[:3]) or "this project (directory-based steps may share intermediate files)")
                        + f" (run {other.run_id}, status {other.status}) — cancel it or wait for it to finish"
                    )
            run = Run(uuid.uuid4().hex[:12], project_id, flow_id, label, steps)
            run.lease = execution.project_lease(storage.project_dir(project_id))
            run.on_step_done = on_step_done
            run.on_done = on_done
            self.runs[run.run_id] = run
        self._persist(run)
        t = threading.Thread(target=self._execute, args=(run,), daemon=True)
        try:
            t.start()
        except Exception:
            run.worker_active = False
            run.lease.close()
            raise
        return run

    def get(self, run_id: str) -> Run | None:
        with self.lock:
            return self.runs.get(run_id)

    def list_for_project(self, project_id: str) -> list:
        with self.lock:
            out = [r for r in self.runs.values() if r.project_id == project_id]
        out.sort(key=lambda r: r.created_at, reverse=True)
        return out

    def _log_path(self, run: Run) -> Path:
        return storage.project_dir(run.project_id) / "runs" / f"{run.run_id}.log"

    def _meta_path(self, run: Run) -> Path:
        return storage.project_dir(run.project_id) / "runs" / f"{run.run_id}.json"

    def _persist(self, run: Run) -> None:
        try:
            with run.cond:
                path = self._meta_path(run)
                tmp = path.with_suffix(".json.tmp")
                tmp.write_text(json.dumps(run.snapshot(), indent=2))
                os.replace(tmp, path)
        except Exception:
            # The in-memory state and runs/<id>.log still exist, but hiding
            # the failure entirely would let the persisted record silently
            # diverge from reality.
            traceback.print_exc()

    def _append_log(self, run: Run, line: str) -> None:
        try:
            with open(self._log_path(run), "a") as f:
                f.write(line + "\n")
        except Exception:
            traceback.print_exc()

    # ---- smart skipping -------------------------------------------------
    # A step is skipped only when its outputs exist AND a stored fingerprint
    # (command + identity of every existing input file it reads) matches.
    # Any change to the command or to an upstream file invalidates the
    # fingerprint and forces recomputation of the affected chain.
    @staticmethod
    def _step_fingerprint(run: Run, step: dict) -> str | None:
        return execution.fingerprint(step, storage.project_dir(run.project_id))

    def _run_with_outputs(self, run, step):
        return execution.run_isolated(run, step, storage.project_dir(run.project_id), self._run_step)

    @staticmethod
    def _sig_path(run: Run, output: str) -> Path:
        return storage.project_dir(run.project_id) / f"{output}.sig"

    def _step_cached(self, run: Run, step: dict) -> bool:
        return execution.step_cached(step, storage.project_dir(run.project_id))

    def _record_sigs(self, run: Run, step: dict, fp: str | None = None) -> None:
        execution.record_sigs(step, storage.project_dir(run.project_id), fp)

    def _clear_sigs(self, run: Run, step: dict) -> None:
        proj_dir = storage.project_dir(run.project_id)
        for o in step.get("outputs") or []:
            try:
                (proj_dir / f"{o}.sig").unlink(missing_ok=True)
            except OSError:
                pass

    def _record_meaning(self, run: Run, step: dict) -> None:
        # Meaning sidecars: <tensor>.meaning.json records the letters-axis
        # meanings and chain provenance the compiler attached to the output,
        # so Reuse Output nodes can restore them without a graph walk. The
        # recorded size lets a stale sidecar (tensor regenerated outside the
        # engine) be detected and ignored.
        tm = step.get("tensor_meaning") or {}
        if not tm:
            return
        proj_dir = storage.project_dir(run.project_id)
        for out_rel, meaning in tm.items():
            try:
                tensor = proj_dir / out_rel
                if not tensor.exists():
                    continue
                sidecar = proj_dir / (out_rel + ".meaning.json")
                payload = {"size": tensor.stat().st_size, "meaning": meaning}
                tmp = sidecar.with_suffix(sidecar.suffix + ".tmp")
                tmp.write_text(json.dumps(payload))
                os.replace(tmp, sidecar)
            except OSError:
                pass

    def _execute(self, run: Run) -> None:
        try:
            self._execute_steps(run)
        except Exception as exc:
            line = f"Run failed: {type(exc).__name__}: {exc}"
            run.emit("log", {"step_id": run.current_step, "stream": "stderr", "line": line})
            self._append_log(run, line)
            for step in run.steps:
                if step.get("status") == "running":
                    step.update(status="failed", returncode=125)
                    self._clear_sigs(run, step)
        finally:
            with run.cond:
                proc = run.proc
            if proc is not None:
                kill_process_tree(proc)
                proc.wait()
            with run.cond:
                run.proc = None
                run.current_step = None
                if _set_status(run, "cancelled" if run.cancel_requested else "failed"):
                    run.emit("status", {"status": run.status})
                    run.emit("end", {})
            self._persist(run)
            self._evict_finished()
            if run.on_done:
                try:
                    run.on_done(run)
                except Exception:
                    traceback.print_exc()
            with self.lock:
                if getattr(run, 'lease', None):
                    run.lease.close()
                run.worker_active = False

    def _execute_steps(self, run: Run) -> None:
        if run.cancel_requested:
            # The wrapper acknowledges cancellation after cleanup.
            self._persist(run)
            return
        if _set_status(run, "running"):
            run.emit("status", {"status": "running"})
        self._persist(run)
        failed = False
        for step in run.steps:
            if run.cancel_requested:
                break
            step_id = step["id"]
            run.current_step = step_id
            outputs = step.get("outputs") or []
            if step.get("skip_if_exists") and outputs and self._step_cached(run, step):
                step["status"] = "skipped"
                # A cached output is still the tensor this step's meaning
                # metadata describes — write the sidecar so Reuse Output
                # nodes can pick the provenance up even on fully-cached runs.
                self._record_meaning(run, step)
                run.emit("step", {"step_id": step_id, "status": "skipped"})
                self._persist(run)
                continue
            step["status"] = "running"
            run.emit("step", {"step_id": step_id, "status": "running"})
            self._persist(run)
            pre_fp = self._step_fingerprint(run, step)
            rc = self._run_with_outputs(run, step)
            if rc == 0:
                # A zero exit code does not prove the step did its job: every
                # declared output must exist (files nonempty), otherwise the
                # failure only surfaces (maybe) as a confusing file-not-found
                # in a later step. Directory outputs are legal — the compiler
                # declares scratch roots (output/.derived/...) as markers.
                _pd = storage.project_dir(run.project_id)
                missing = [
                    o for o in (step.get("outputs") or [])
                    if not (_pd / o).exists()
                    or ((_pd / o).is_file() and (_pd / o).stat().st_size == 0)
                ]
                if missing:
                    line = (f"Step exited 0 but its declared outputs are missing or empty: "
                            + ", ".join(missing))
                    run.emit("log", {"step_id": step_id, "stream": "stderr", "line": line})
                    self._append_log(run, line)
                    rc = 125
            step["returncode"] = rc
            if rc == 0:
                step["status"] = "done"
                self._record_sigs(run, step, pre_fp)
                self._record_meaning(run, step)
                if run.on_step_done:
                    run.on_step_done(run, step, step.get("result"))
                run.emit("step", {"step_id": step_id, "status": "done"})
            elif run.cancel_requested:
                step["status"] = "failed"
                self._clear_sigs(run, step)
                run.emit("step", {"step_id": step_id, "status": "failed"})
                break
            else:
                step["status"] = "failed"
                self._clear_sigs(run, step)
                run.emit("step", {"step_id": step_id, "status": "failed"})
                failed = True
            self._persist(run)
            if failed:
                break
        run.current_step = None
        # Exactly one terminal transition wins (the other may have been the
        # cancel() thread), so the UI never sees conflicting final statuses.
        with run.cond:
            if _set_status(run, "cancelled" if run.cancel_requested else "failed" if failed else "done"):
                run.emit("status", {"status": run.status})
                run.emit("end", {})
        self._persist(run)
        self._evict_finished()
    def _evict_finished(self) -> None:
        # Bound the engine dict; runs/*.json records on disk are unaffected.
        with self.lock:
            finished = [r for r in self.runs.values() if r.status in _TERMINAL_STATUSES]
            finished.sort(key=lambda r: r.created_at)
            while len(finished) > MAX_FINISHED_RUNS:
                victim = finished.pop(0)
                self.runs.pop(victim.run_id, None)

    def _run_step(self, run: Run, step: dict) -> int:
        argv = step["argv"]
        cwd = step.get("cwd")
        proj_dir = storage.project_dir(run.project_id)
        for o in step.get("outputs") or []:
            try:
                (proj_dir / o).parent.mkdir(parents=True, exist_ok=True)
            except OSError:
                traceback.print_exc()
        env = os.environ.copy()
        try:
            proc = subprocess.Popen(
                argv,
                cwd=cwd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                bufsize=1,
                start_new_session=POSIX_SESSIONS,
            )
        except Exception as exc:
            line = f"Error launching command: {exc}"
            run.emit("log", {"step_id": step["id"], "stream": "stderr", "line": line})
            self._append_log(run, line)
            return 127
        with run.cond:
            run.proc = proc
            if run.cancel_requested:
                kill_process_tree(proc)
        # Optional per-step wall-clock guard: without it a hung binary blocks
        # the run thread forever and nothing short of a manual cancel notices.
        # Disabled by default (exact runs can legitimately take hours/days);
        # set JOBS_STEP_TIMEOUT=<seconds> to arm it.
        try:
            timeout = float(os.environ.get("JOBS_STEP_TIMEOUT", "0") or 0)
        except ValueError:
            timeout = 0.0
        timed_out = threading.Event()

        def _kill_on_timeout():
            timed_out.set()
            kill_process_tree(proc)

        killer = threading.Timer(timeout, _kill_on_timeout) if timeout > 0 else None
        if killer is not None:
            killer.daemon = True
            killer.start()
        result_payload = None
        assert proc.stdout is not None
        try:
            for line in proc.stdout:
                line = line.rstrip("\n")
                marker = parse_result_marker(line)
                if marker is not None:
                    result_payload = marker
                    step["result"] = marker
                run.emit("log", {"step_id": step["id"], "stream": "stdout", "line": line})
                self._append_log(run, line)
            rc = proc.wait()
        finally:
            if killer is not None:
                killer.cancel()
            if proc.poll() is None:
                kill_process_tree(proc)
                proc.wait()
            proc.stdout.close()
            with run.cond:
                run.proc = None
        if timed_out.is_set():
            line = (f"Step exceeded JOBS_STEP_TIMEOUT={timeout:g}s and was killed "
                    "(the run is marked failed; raise or unset JOBS_STEP_TIMEOUT for long steps)")
            run.emit("log", {"step_id": step["id"], "stream": "stderr", "line": line})
            self._append_log(run, line)
            return 124
        return rc

    def cancel(self, run_id: str) -> bool:
        run = self.get(run_id)
        if run is None:
            return False
        with run.cond:
            if run.status in _TERMINAL_STATUSES or run.cancel_requested:
                return False
            run.cancel_requested = True
            proc = run.proc
        if proc is not None:
            kill_process_tree(proc)
        # The worker acknowledges cancellation only after the process exits,
        # output rollback completes, and its reservations can be released.
        return True


engine = Engine()
