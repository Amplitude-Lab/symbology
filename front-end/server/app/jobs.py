from __future__ import annotations

import hashlib
import json
import os
import signal
import subprocess
import threading
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from . import storage
from .wolfram import parse_result_marker


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


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
        self.events = []
        self.cond = threading.Condition()
        self.proc: subprocess.Popen | None = None
        self.on_step_done = None
        self.on_done = None

    def emit(self, event: str, data: dict) -> None:
        with self.cond:
            self.events.append({"event": event, "data": data})
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
        run = Run(uuid.uuid4().hex[:12], project_id, flow_id, label, steps)
        run.on_step_done = on_step_done
        run.on_done = on_done
        with self.lock:
            self.runs[run.run_id] = run
        self._persist(run)
        t = threading.Thread(target=self._execute, args=(run,), daemon=True)
        t.start()
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
            self._meta_path(run).write_text(json.dumps(run.snapshot(), indent=2))
        except Exception:
            pass

    def _append_log(self, run: Run, line: str) -> None:
        try:
            with open(self._log_path(run), "a") as f:
                f.write(line + "\n")
        except Exception:
            pass

    # ---- smart skipping -------------------------------------------------
    # A step is skipped only when its outputs exist AND a stored fingerprint
    # (command + identity of every existing input file it reads) matches.
    # Any change to the command or to an upstream file invalidates the
    # fingerprint and forces recomputation of the affected chain.
    @staticmethod
    def _step_fingerprint(run: Run, step: dict) -> str | None:
        proj_dir = storage.project_dir(run.project_id)
        outputs = set(step.get("outputs") or [])
        cwd = step.get("cwd")
        h = hashlib.sha256()
        h.update(step.get("command", "").encode())
        h.update("\x00".join(str(a) for a in (step.get("argv") or [])).encode())
        for arg in step.get("argv") or []:
            p = Path(arg)
            if not p.is_absolute():
                if cwd is None:
                    continue
                p = Path(cwd) / p
            try:
                rp = p.resolve().relative_to(proj_dir.resolve())
            except ValueError:
                rp = None
            rel = rp.as_posix() if rp is not None else str(p)
            if rel in outputs or not p.exists() or p.is_dir():
                continue
            try:
                st = p.stat()
                h.update(f"{rel}:{st.st_size}:{st.st_mtime_ns}:{st.st_ino};".encode())
            except OSError:
                h.update(f"{rel}:?;".encode())
        return h.hexdigest()

    @staticmethod
    def _sig_path(run: Run, output: str) -> Path:
        return storage.project_dir(run.project_id) / f"{output}.sig"

    def _step_cached(self, run: Run, step: dict) -> bool:
        outputs = step.get("outputs") or []
        if not outputs:
            return False
        proj_dir = storage.project_dir(run.project_id)
        for o in outputs:
            if not (proj_dir / o).exists():
                return False
        fp = self._step_fingerprint(run, step)
        if fp is None:
            return False
        for o in outputs:
            sp = self._sig_path(run, o)
            try:
                if sp.read_text().strip() != fp:
                    return False
            except OSError:
                return False
        return True

    def _record_sigs(self, run: Run, step: dict, fp: str | None = None) -> None:
        if fp is None:
            fp = self._step_fingerprint(run, step)
        if fp is None:
            return
        proj_dir = storage.project_dir(run.project_id)
        for o in step.get("outputs") or []:
            try:
                sp = proj_dir / f"{o}.sig"
                tmp = sp.with_suffix(sp.suffix + ".tmp")
                tmp.write_text(fp)
                os.replace(tmp, sp)
            except OSError:
                pass

    def _clear_sigs(self, run: Run, step: dict) -> None:
        proj_dir = storage.project_dir(run.project_id)
        for o in step.get("outputs") or []:
            try:
                (proj_dir / f"{o}.sig").unlink(missing_ok=True)
            except OSError:
                pass

    def _execute(self, run: Run) -> None:
        if run.status == "cancelled":
            run.emit("end", {})
            self._persist(run)
            return
        run.status = "running"
        run.emit("status", {"status": "running"})
        self._persist(run)
        failed = False
        for step in run.steps:
            if run.status == "cancelled":
                break
            step_id = step["id"]
            run.current_step = step_id
            outputs = step.get("outputs") or []
            if step.get("skip_if_exists") and outputs and self._step_cached(run, step):
                step["status"] = "skipped"
                run.emit("step", {"step_id": step_id, "status": "skipped"})
                self._persist(run)
                continue
            step["status"] = "running"
            run.emit("step", {"step_id": step_id, "status": "running"})
            self._persist(run)
            pre_fp = self._step_fingerprint(run, step)
            rc = self._run_step(run, step)
            step["returncode"] = rc
            if rc == 0:
                step["status"] = "done"
                self._record_sigs(run, step, pre_fp)
                run.emit("step", {"step_id": step_id, "status": "done"})
            elif run.status == "cancelled":
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
        if run.status != "cancelled":
            run.status = "failed" if failed else "done"
        run.emit("status", {"status": run.status})
        run.emit("end", {})
        self._persist(run)
        if run.on_done:
            try:
                run.on_done(run)
            except Exception:
                pass

    def _run_step(self, run: Run, step: dict) -> int:
        argv = step["argv"]
        cwd = step.get("cwd")
        proj_dir = storage.project_dir(run.project_id)
        for o in step.get("outputs") or []:
            try:
                (proj_dir / o).parent.mkdir(parents=True, exist_ok=True)
            except OSError:
                pass
        env = os.environ.copy()
        try:
            proc = subprocess.Popen(
                argv,
                cwd=cwd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                start_new_session=True,
            )
        except Exception as exc:
            line = f"Error launching command: {exc}"
            run.emit("log", {"step_id": step["id"], "stream": "stderr", "line": line})
            self._append_log(run, line)
            return 127
        run.proc = proc
        result_payload = None
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            marker = parse_result_marker(line)
            if marker is not None:
                result_payload = marker
                step["result"] = marker
            run.emit("log", {"step_id": step["id"], "stream": "stdout", "line": line})
            self._append_log(run, line)
        rc = proc.wait()
        run.proc = None
        if rc == 0 and run.on_step_done:
            try:
                run.on_step_done(run, step, result_payload)
            except Exception:
                pass
        return rc

    def cancel(self, run_id: str) -> bool:
        run = self.get(run_id)
        if run is None or run.status not in ("queued", "running"):
            return False
        run.status = "cancelled"
        if run.proc is not None:
            try:
                os.killpg(os.getpgid(run.proc.pid), signal.SIGTERM)
            except Exception:
                pass
            try:
                run.proc.wait(timeout=5)
            except Exception:
                try:
                    os.killpg(os.getpgid(run.proc.pid), signal.SIGKILL)
                except Exception:
                    pass
        run.emit("status", {"status": "cancelled"})
        run.emit("end", {})
        self._persist(run)
        return True


engine = Engine()
