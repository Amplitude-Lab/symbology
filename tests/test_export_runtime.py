"""Standalone export must obey the same publication and cache contract as local runs."""
import copy
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from test_robustness import project, step, ROOT
from app import compile as compiler, execution, storage


def export(project, monkeypatch, steps, name="flow"):
    monkeypatch.setattr(compiler, 'compile_flow', lambda *a, **kw: dict(ok=True, _steps_full=copy.deepcopy(steps), flow_outputs=[]))
    result = compiler.export_flow_script(project, {}, name)
    assert result['ok'], result
    return Path(result['path'])


def run(script, **env):
    return subprocess.run(['bash', str(script)], capture_output=True, text=True, timeout=10,
                          env={**os.environ, **env})


def test_export_relocated_spaces_apostrophes_and_changed_seed(project, monkeypatch, tmp_path):
    root=storage.project_dir(project['id']); (root/'data/input with space').write_text('one')
    code="from pathlib import Path; Path('output/result with space').write_text(Path('data/input with space').read_text())"
    s=step(project, code, ['output/result with space']);s['input_dirs']=[str(root/'data')]
    script=export(project,monkeypatch,[s],"flow's name\n$(touch NEVER)")
    dest=tmp_path/"moved project with ' quote";shutil.copytree(root,dest)
    moved=dest/'exported'/script.name
    first=run(moved);assert first.returncode==0,first.stderr
    assert (dest/'output/result with space').read_text()=='one'
    assert 'CACHED' in run(moved).stdout
    inp=dest/'data/input with space';stat=inp.stat();inp.write_text('two');os.utime(inp,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    again=run(moved);assert again.returncode==0,again.stderr
    assert (dest/'output/result with space').read_text()=='two'
    assert not (root/'output/result with space').exists()


@pytest.mark.parametrize('code', ["pass", "from pathlib import Path; Path('output/result').write_text('partial'); raise SystemExit(2)"])
def test_export_failure_preserves_old_output(project, monkeypatch, code):
    root=storage.project_dir(project['id']);(root/'output/result').write_text('old')
    script=export(project,monkeypatch,[step(project,code,['output/result'])])
    result=run(script);assert result.returncode!=0
    assert (root/'output/result').read_text()=='old'
    assert not list((root/'runs').glob('.attempt-*'))


def test_export_timeout_preserves_old_output_and_cleans_stage(project, monkeypatch):
    root=storage.project_dir(project['id']);(root/'output/result').write_text('old')
    s=step(project,"from pathlib import Path; import time; Path('output/result').write_text('partial'); time.sleep(30)",['output/result'])
    result=run(export(project,monkeypatch,[s]),STEP_TIMEOUT='0.15')
    assert result.returncode!=0
    assert 'timed out' in result.stderr
    assert (root/'output/result').read_text()=='old'
    assert not list((root/'runs').glob('.attempt-*'))


def test_export_wolfram_skip_requires_unchanged_complete_files(project, monkeypatch):
    root=storage.project_dir(project['id']);(root/'data/seed').write_text('input')
    (root/'output/result').write_text('ready')
    fake=root/'fake-wolfram';fake.write_text('#!/bin/sh\nexit 0\n');fake.chmod(0o755)
    s=step(project,outputs=['output/result']);s.update(kind='wolfram',argv=[str(fake)],inputs=[str(fake),str(root/'data/seed')])
    execution.record_sigs(s,root)
    script=export(project,monkeypatch,[s])
    env={'WOLFRAMSCRIPT':'/no/such/wolfram'}
    result=run(script,**env);assert result.returncode==0,result.stderr
    assert run(script,**env,WOLFRAM_MODE='fail').returncode!=0
    (root/'data/seed').write_text('changed')
    assert run(script,**env,WOLFRAM_MODE='skip').returncode!=0
    (root/'output/result').unlink()
    assert run(script,**env,WOLFRAM_MODE='skip').returncode!=0


def test_export_wolfram_program_paths_relocate_only_in_stage(project, monkeypatch, tmp_path):
    root=storage.project_dir(project['id'])
    program=root/'data/program.wl';program.write_text(str(root/'data/seed'))
    fake=tmp_path/'fake-wolfram';fake.write_text('#!/usr/bin/env python3\nimport sys\nfrom pathlib import Path\nPath("output/result").write_text(Path(sys.argv[-1]).read_text())\n');fake.chmod(0o755)
    s=step(project,outputs=['output/result']);s.update(kind='wolfram',argv=[str(fake),'-script',str(program)])
    script=export(project,monkeypatch,[s]);dest=tmp_path/'relocated';shutil.copytree(root,dest)
    result=run(dest/'exported'/script.name,WOLFRAMSCRIPT=fake.name,PATH=str(fake.parent)+os.pathsep+os.environ['PATH']);assert result.returncode==0,result.stderr
    assert str(dest/'runs/.attempt-') in (dest/'output/result').read_text()
    assert (dest/'data/program.wl').read_text()==program.read_text()


def test_export_respects_local_project_lease(project, monkeypatch):
    root=storage.project_dir(project['id'])
    s=step(project,"from pathlib import Path; Path('output/result').write_text('new')",['output/result'])
    script=export(project,monkeypatch,[s])
    with execution.project_lease(root):
        result=run(script)
        assert result.returncode!=0 and 'using this project' in result.stderr
        assert not (root/'output/result').exists()
    assert run(script).returncode==0


def test_export_publishes_removal_of_stale_derived_output(project, monkeypatch):
    root=storage.project_dir(project['id']);(root/'output/obsolete').write_text('old')
    code="from pathlib import Path; Path('output/obsolete').unlink(); Path('output/current').write_text('new')"
    s=step(project,code,['output/current']);s['input_dirs']=[str(root/'output')]
    result=run(export(project,monkeypatch,[s]));assert result.returncode==0,result.stderr
    assert not (root/'output/obsolete').exists()
    assert (root/'output/current').read_text()=='new'


def test_removed_outputs_are_restored_after_publication_failure(project, monkeypatch):
    import threading
    from types import SimpleNamespace
    root=storage.project_dir(project['id'])
    for name in ('a-old','b-old'): (root/'output'/name).write_text(name)
    s=step(project,outputs=['output/new']);s['input_dirs']=[str(root/'output')]
    attempt=SimpleNamespace(run_id='removal-test',cancel_requested=False,cond=threading.Condition())
    def launch(run, isolated):
        stage=Path(isolated['cwd'])
        for name in ('a-old','b-old'): (stage/'output'/name).unlink()
        (stage/'output/new').write_text('new')
        return 0
    original=Path.unlink
    def fail_second(path,*args,**kwargs):
        if path == root/'output/b-old': raise OSError('simulated publication failure')
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'unlink',fail_second)
    with pytest.raises(OSError,match='publication failure'):
        execution.run_isolated(attempt,s,root,launch)
    for name in ('a-old','b-old'): assert (root/'output'/name).read_text()==name
    assert not (root/'output/new').exists()
