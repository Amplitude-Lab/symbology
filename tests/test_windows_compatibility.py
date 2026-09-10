"""Portable naming checks plus native-Windows launcher acceptance tests.

The launcher tests require cmd.exe and run only on Windows. They launch a
disposable Python server fixture, never the user's server or project store.
"""
import copy
import os
from pathlib import Path, PureWindowsPath
import shutil
import subprocess
import sys
import venv

import pytest

from test_robustness import ROOT, project, step
from app import compile as compiler, execution, storage


@pytest.mark.parametrize('name', ['NUL', 'con', 'AUX', 'COM1', 'LPT9'])
def test_project_names_can_be_copied_to_windows(project, name):
    created = storage.create_project(name)
    assert not PureWindowsPath(created['id']).is_reserved()
    assert storage.load_project(created['id'])['name'] == name


@pytest.mark.parametrize('name', ['NUL', 'con.extra', 'LPT1'])
def test_export_names_and_line_endings_are_portable(project, monkeypatch, name):
    plan = [step(project)]
    monkeypatch.setattr(compiler, 'compile_flow', lambda *a, **kw:
                        dict(ok=True, _steps_full=copy.deepcopy(plan), flow_outputs=[]))
    result = compiler.export_flow_script(project, {}, name)
    assert result['ok'], result
    script = Path(result['path'])
    assert not PureWindowsPath(script.name).is_reserved()
    content = script.read_bytes()
    assert b'\r\n' not in content  # Bash scripts copied to WSL must keep LF.
    assert content.decode('utf-8').startswith('#!/usr/bin/env bash\n')


def test_project_lease_excludes_another_process_and_releases(project):
    root = storage.project_dir(project['id'])
    code = (
        'import sys; from pathlib import Path; from app.execution import project_lease\n'
        'try:\n'
        '    with project_lease(Path(sys.argv[1])): pass\n'
        'except ValueError:\n'
        '    sys.exit(23)\n'
    )
    def contender():
        return subprocess.run([sys.executable, '-c', code, str(root)],
                              env={**os.environ, 'PYTHONPATH': str(ROOT / 'front-end/server')},
                              capture_output=True, text=True, timeout=10)
    with execution.project_lease(root):
        held = contender()
        assert held.returncode == 23, held.stderr
    released = contender()
    assert released.returncode == 0, released.stderr


@pytest.fixture
def windows_launcher(tmp_path):
    if os.name != 'nt':
        pytest.skip('Requires native Windows cmd.exe and Python')
    front = tmp_path / "Studio (spaces) & user's copy" / 'front-end'
    server = front / 'server'
    server.mkdir(parents=True)
    (front / 'web/dist').mkdir(parents=True)
    (front / 'web/dist/index.html').write_text('<html>fixture</html>', encoding='utf-8')
    shutil.copy2(ROOT / 'front-end/start.bat', front / 'start.bat')
    (server / 'requirements.txt').write_text('', encoding='utf-8')
    (server / 'run.py').write_text(
        'import os, sys\nfrom pathlib import Path\n'
        'Path(os.environ["LAUNCH_MARKER"]).write_text("started")\n'
        'sys.exit(int(os.environ.get("SERVER_EXIT", "0")))\n', encoding='utf-8')
    venv.EnvBuilder(system_site_packages=True, with_pip=False).create(server / '.venv')
    fake_bin = tmp_path / 'fake commands'
    fake_bin.mkdir()
    (fake_bin / 'npm.cmd').write_text(
        '@echo off\n'
        'echo %*>>"%NPM_TRACE%"\n'
        'if "%1"=="ci" exit /b %NPM_CI_EXIT%\n'
        'if "%1"=="install" exit /b %NPM_CI_EXIT%\n'
        'if not "%NPM_BUILD_EXIT%"=="0" exit /b %NPM_BUILD_EXIT%\n'
        'if not exist dist mkdir dist\n'
        'echo fixture>dist\\index.html\n'
        'exit /b 0\n', encoding='ascii', newline='\r\n')
    marker, trace = tmp_path / 'launched', tmp_path / 'npm-trace'
    env = {**os.environ, 'SYMBOLOGY_NO_BROWSER': '1', 'LAUNCH_MARKER': str(marker),
           'NPM_TRACE': str(trace), 'NPM_CI_EXIT': '0', 'NPM_BUILD_EXIT': '0',
           'PIP_NO_INDEX': '1', 'PIP_DISABLE_PIP_VERSION_CHECK': '1',
           'PIP_CONFIG_FILE': os.devnull,
           'PATH': str(fake_bin) + os.pathsep + os.environ['PATH']}
    def launch(**overrides):
        return subprocess.run([os.environ.get('COMSPEC', 'cmd.exe'), '/d', '/c',
                               'call', str(front / 'start.bat')],
                              cwd=tmp_path, env={**env, **overrides},
                              capture_output=True, text=True, timeout=60)
    return front, marker, trace, launch


def test_windows_launcher_handles_spaces_and_server_exit(windows_launcher):
    _, marker, _, launch = windows_launcher
    result = launch(SERVER_EXIT='7')
    assert marker.exists(), result.stdout + result.stderr
    assert result.returncode == 7


def test_windows_launcher_repairs_an_incomplete_venv(windows_launcher):
    front, marker, _, launch = windows_launcher
    shutil.rmtree(front / 'server/.venv')
    (front / 'server/.venv').mkdir()
    result = launch()
    assert result.returncode == 0, result.stdout + result.stderr
    assert marker.exists()


def test_windows_launcher_stops_on_dependency_failure(windows_launcher):
    front, marker, _, launch = windows_launcher
    (front / 'server/requirements.txt').write_text('symbology-nonexistent-test-package==0\n')
    result = launch()
    assert result.returncode != 0, result.stdout + result.stderr
    assert not marker.exists()


@pytest.mark.parametrize('failure', ['NPM_CI_EXIT', 'NPM_BUILD_EXIT'])
def test_windows_launcher_stops_on_web_build_failure(windows_launcher, failure):
    front, marker, trace, launch = windows_launcher
    shutil.rmtree(front / 'web/dist')
    result = launch(**{failure: '9'})
    assert result.returncode != 0, result.stdout + result.stderr
    assert not marker.exists()
    calls = trace.read_text().splitlines()
    assert calls[0] == 'ci --no-audit --no-fund'
    assert len(calls) == (1 if failure == 'NPM_CI_EXIT' else 2)


def test_windows_launcher_builds_missing_index(windows_launcher):
    front, marker, trace, launch = windows_launcher
    (front / 'web/dist/index.html').unlink()
    result = launch()
    assert result.returncode == 0, result.stdout + result.stderr
    assert marker.exists() and (front / 'web/dist/index.html').is_file()
    assert trace.read_text().splitlines() == ['ci --no-audit --no-fund', 'run build']
