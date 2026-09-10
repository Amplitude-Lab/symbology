"""Standard-library standalone runner, embedded with execution.py on export."""
import json
from contextlib import nullcontext
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid


def run_export(payload):
    root = Path(os.environ.get('PROJ_DIR') or Path(os.environ['SYMBOLOGY_EXPORT_SCRIPT_DIR']).parent).resolve()
    repo = Path(os.environ.get('SYMBOLOGY_ROOT') or payload['repo']).resolve()
    requested_wolfram = os.environ.get('WOLFRAMSCRIPT') or shutil.which('wolframscript') or payload['wolfram']
    wolfram = shutil.which(requested_wolfram) or str((root / requested_wolfram).resolve())
    mode = os.environ.get('WOLFRAM_MODE', 'auto')
    if mode not in ('auto', 'skip', 'fail'):
        raise ValueError('WOLFRAM_MODE must be auto, skip or fail')
    timeout = float(os.environ.get('STEP_TIMEOUT', '0'))
    if not math.isfinite(timeout) or timeout < 0:
        raise ValueError('STEP_TIMEOUT must be a finite nonnegative number of seconds')
    selected = set()
    for part in filter(None, os.environ.get('STEPS', '').split(',')):
        if not re.fullmatch(r'[1-9][0-9]*(?:-[1-9][0-9]*)?', part):
            raise ValueError('STEPS must contain step numbers or ranges, such as 1,3-5')
        ends = list(map(int, part.split('-')))
        lo, hi = ends[0], ends[-1]
        if hi < lo or hi > len(payload['steps']):
            raise ValueError('STEPS contains an invalid range')
        selected.update(range(lo, hi + 1))
    if sys.argv[1:] not in ([], ['--dry-run']):
        raise ValueError('Usage: exported-flow.sh [--dry-run]')
    dry_run = '--dry-run' in sys.argv
    replacements = [(payload['project'], str(root)), (payload['repo'], str(repo))]
    replacements.sort(key=lambda pair: len(pair[0]), reverse=True)

    def relocate(value):
        for old, new in replacements:
            if value == old or value.startswith(old + os.sep):
                return new + value[len(old):]
        return value

    class Attempt:
        run_id = 'export-' + uuid.uuid4().hex[:12]
        cancel_requested = False
        cond = threading.Condition()
        proc = None

    run = Attempt()
    (root / 'runs').mkdir(parents=True, exist_ok=True)
    log_path = root / 'runs' / (run.run_id + '.log')

    def stop(_sig=None, _frame=None):
        run.cancel_requested = True
        with run.cond:
            if run.proc is not None:
                try:
                    if os.name == 'posix':
                        os.killpg(run.proc.pid, signal.SIGKILL)
                    else:
                        run.proc.kill()
                except ProcessLookupError:
                    pass

    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, stop)

    with (nullcontext() if dry_run else execution.project_lease(root)), log_path.open('w') as log:
        def note(line):
            print(line, flush=True)
            log.write(line + '\n')
            log.flush()

        def launch(run, step):
            with run.cond:
                if run.cancel_requested:
                    return 130
                proc = subprocess.Popen(step['argv'], cwd=step['cwd'], stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT, text=True, encoding='utf-8',
                                        errors='replace', start_new_session=os.name == 'posix')
                run.proc = proc
            timer = threading.Timer(timeout, stop) if timeout else None
            if timer:
                timer.start()
            try:
                for line in proc.stdout:
                    note(line.rstrip('\n'))
                return proc.wait()
            finally:
                if timer:
                    timer.cancel()
                    timer.join()
                if proc.poll() is None:
                    stop()
                    proc.wait()
                proc.stdout.close()
                with run.cond:
                    run.proc = None

        note(f'Project: {root}\nRepository: {repo}\nLog: {log_path}')
        for n, original in enumerate(payload['steps'], 1):
            if selected and n not in selected:
                note(f'[{n}] SKIPPED by STEPS filter')
                continue
            step = dict(original)
            step['argv'] = [relocate(a) for a in step['argv']]
            step['cwd'] = str(root)
            for key in ('inputs', 'input_dirs'):
                if key in step:
                    step[key] = [relocate(p) for p in step[key]]
            step['path_replacements'] = replacements
            if step['kind'] == 'wolfram':
                old_exe = step['argv'][0]
                step['argv'][0] = wolfram
                step['inputs'] = [wolfram if p == old_exe else p for p in step.get('inputs', [])]
            label = f'[{n}/{len(payload["steps"])}] {step["label"]}'
            if dry_run:
                note(f'DRY-RUN {label}: {shlex.join(step["argv"])}')
                continue
            if run.cancel_requested:
                raise RuntimeError('Run cancelled')
            available = shutil.which(step['argv'][0])
            if step['kind'] == 'wolfram' and not available:
                trust = original.get('precomputed', {})
                inputs = {str(p): execution.digest(p) if p.is_file() else 'missing'
                          for p in execution.input_files(step, root)
                          if str(p) != str(Path(wolfram).resolve()) and str(p) != wolfram}
                expected = {relocate(p): value for p, value in trust.get('inputs', {}).items()}
                outputs_ok = bool(step.get('outputs')) and all(
                    (root / o).is_file() and (root / o).stat().st_size and
                    execution.digest(root / o) == trust.get('outputs', {}).get(o)
                    for o in step.get('outputs', []))
                if mode != 'fail' and (mode == 'skip' or trust.get('verified')) and outputs_ok and inputs == expected:
                    note(f'{label}: using verified unchanged precomputed files (Wolfram unavailable)')
                    continue
                raise RuntimeError(f'{label}: Wolfram unavailable; unchanged precomputed inputs and outputs are required. Copy the computed project or install Wolfram.')
            if not available:
                raise RuntimeError(f'{label}: executable not found: {step["argv"][0]}')
            if step.get('skip_if_exists') and execution.step_cached(step, root):
                note(f'{label}: CACHED (input and output contents match)')
                continue
            before = execution.fingerprint(step, root)
            note(f'{label}: {shlex.join(step["argv"])}')
            started = time.monotonic()
            rc = execution.run_isolated(run, step, root, launch)
            if rc:
                raise RuntimeError(f'{label}: exit {rc}' + (' (cancelled or timed out)' if run.cancel_requested else ''))
            execution.record_sigs(step, root, before)
            note(f'{label}: done in {time.monotonic() - started:.2f}s')
        note('Selected steps complete.' if selected else 'ALL STEPS COMPLETE')
