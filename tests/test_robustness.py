"""Adversarial audit: assertions describe desired behavior, so defects FAIL.

Run with pytest; requires the server requirements, httpx, pytest.
All project data and executable fixtures are isolated under pytest's tmp_path.
Fake executables test orchestration only, never numerical correctness.
"""
import asyncio
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "front-end/server"))
from app import storage, jobs, compile as compiler, main
from fastapi.testclient import TestClient


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setattr(storage, "PROJECTS_DIR", tmp_path / "projects")
    storage.init_dirs()
    return storage.create_project("Audit")


@pytest.fixture
def client(project):
    with TestClient(main.app, raise_server_exceptions=False) as c:
        yield c


def step(project, code="pass", outputs=None):
    argv = [sys.executable, "-c", code]
    return dict(id="step-1", label="fixture", kind="fixture", argv=argv,
                command=" ".join(argv), cwd=str(storage.project_dir(project["id"])),
                outputs=outputs or [], skip_if_exists=True)


def run_for(project, s):
    return jobs.Run("abcdef123456", project["id"], None, "audit", [s])


def matrix_graph(n=2):
    return {"nodes": [
        {"id": "a", "type": "reuse_output", "data": {"file": "data/M.wxf", "kind": "matrix"}},
        {"id": "b", "type": "matrix_power", "data": {"target": "square", "n": n}},
    ], "edges": [{"id": "e", "source": "a", "sourceHandle": "out", "target": "b", "targetHandle": "matrix"}]}


@pytest.fixture
def matrix_setup(project, tmp_path, monkeypatch):
    (storage.project_dir(project["id"]) / "data/M.wxf").write_bytes(b"input")
    binary = tmp_path / "tensor_ops"
    binary.write_text('#!/bin/sh\nprintf v1 > "$4"\n')
    binary.chmod(0o755)
    monkeypatch.setattr(compiler, "find_tensor_ops", lambda: str(binary))
    return binary


def test_basic_crud(client, project):
    p = client.get(f'/api/projects/{project["id"]}')
    assert p.status_code == 200
    f = client.post(f'/api/projects/{project["id"]}/flows', json={"name": "test"})
    assert f.status_code == 200
    assert client.delete(f'/api/projects/{project["id"]}/flows/{f.json()["id"]}').status_code == 200


@pytest.mark.parametrize("payload", [{"name": 7}, {"name": []}, {"name": {"a": 1}}])
def test_project_invalid_types_have_json_4xx(client, payload):
    r = client.post('/api/projects', json=payload)
    assert 400 <= r.status_code < 500, (r.status_code, r.text)
    assert "detail" in r.json()


def test_alphabet_invalid_letter_type(client, project):
    r = client.post(f'/api/projects/{project["id"]}/alphabets', json={"name": "A", "letters": [1]})
    assert 400 <= r.status_code < 500, (r.status_code, r.text)


def test_bad_graph_rejected_on_save(client, project):
    url = f'/api/projects/{project["id"]}/flows'
    fid = client.post(url, json={"name": "test"}).json()["id"]
    r = client.put(f'{url}/{fid}', json={"revision": 0, "graph": {"nodes": [None], "edges": []}})
    cr = client.post(f'{url}/{fid}/compile', json={})
    assert 400 <= r.status_code < 500, ("save", r.status_code, "compile", cr.status_code, cr.text)


def test_empty_graph_rejected(project):
    assert compiler.compile_flow(project, {})["ok"] is False


def test_graph_cycle_rejected(project, matrix_setup):
    g = matrix_graph()
    g['edges'].append(dict(id="cycle", source="b", sourceHandle="out", target="a", targetHandle="in"))
    assert "cycle" in " ".join(compiler.compile_flow(project, g)["errors"])


def test_duplicate_node_ids_rejected(project, matrix_setup):
    g = matrix_graph()
    g["nodes"].insert(0, {"id": "a", "type": "unknown", "data": {}})
    assert not compiler.compile_flow(project, g)["ok"]


def test_dangling_edge_rejected(project, matrix_setup):
    g = matrix_graph()
    g['edges'].append(dict(id="bad", source="missing", sourceHandle="out", target="b", targetHandle="matrix"))
    assert not compiler.compile_flow(project, g)["ok"]


def test_duplicate_output_paths_rejected(project,matrix_setup):
    g=matrix_graph()
    g['nodes'].append({'id':'c','type':'matrix_power','data':{'target':'square','n':3}})
    g['edges'].append({'id':'e2','source':'a','sourceHandle':'out','target':'c','targetHandle':'matrix'})
    r=compiler.compile_flow(project,g)
    assert not r['ok'],[s['outputs'] for s in r['steps']]


def test_same_input_port_cannot_have_two_sources(project,matrix_setup):
    g=matrix_graph()
    g['nodes'].append({'id':'c','type':'reuse_output','data':{'file':'data/M.wxf','kind':'matrix'}})
    g['edges'].append({'id':'e2','source':'c','sourceHandle':'out','target':'b','targetHandle':'matrix'})
    assert not compiler.compile_flow(project,g)['ok']


@pytest.mark.parametrize("exponent", [0, "0"])
def test_matrix_zero_power(project, matrix_setup, exponent):
    r = compiler.compile_flow(project, matrix_graph(exponent))
    assert r["ok"], r["errors"]


def test_reuse_standard_weight_name(project, monkeypatch):
    p = storage.project_dir(project['id'])
    for f in ['FEC_2.wxf', 'D.wxf']:
        (p/'data'/f).write_bytes(b'fixture')
    monkeypatch.setattr(compiler, 'find_bootstrap', lambda: '/bin/true')
    graph = {'nodes': [
        {'id':'f','type':'reuse_output','data':{'file':'data/FEC_2.wxf','kind':'fec'}},
        {'id':'d','type':'reuse_output','data':{'file':'data/D.wxf','kind':'dlogmat'}},
        {'id':'e','type':'extend','data':{'target_weight':3}},
    ],'edges':[
        {'source':'f','sourceHandle':'out','target':'e','targetHandle':'fec'},
        {'source':'d','sourceHandle':'out','target':'e','targetHandle':'condition'},
    ]}
    r=compiler.compile_flow(project,graph)
    assert r['ok'],r['errors']


def test_missing_output_fails(project):
    s = step(project, outputs=["output/missing.wxf"])
    r = run_for(project, s)
    jobs.Engine()._execute(r)
    assert r.status == "failed" and s['returncode'] == 125


def test_old_output_does_not_mask_failed_recompute(project):
    p = storage.project_dir(project['id'])/'output/old.wxf'
    p.write_bytes(b'old solution')
    s = step(project, outputs=['output/old.wxf'])
    r = run_for(project, s)
    jobs.Engine()._execute(r)
    assert r.status == 'failed', (r.status, s['returncode'], p.read_bytes())


def test_truncated_cached_output_is_not_reused(project):
    p = storage.project_dir(project['id'])/'output/cached.wxf'
    p.write_bytes(b'original')
    s=step(project,outputs=['output/cached.wxf'])
    r=run_for(project,s); e=jobs.Engine()
    e._record_sigs(r,s)
    p.write_bytes(b'')
    assert not e._step_cached(r,s)


def test_changed_direct_input_invalidates_cache(project):
    p=storage.project_dir(project['id'])/'data/input.wxf'
    p.write_bytes(b'a')
    s=step(project); s['argv'].append(str(p))
    r=run_for(project,s); e=jobs.Engine()
    before=e._step_fingerprint(r,s)
    p.write_bytes(b'changed')
    assert before != e._step_fingerprint(r,s)


def test_changed_directory_input_invalidates_cache(project):
    p=storage.project_dir(project['id'])/'data/input.wxf'
    p.write_bytes(b'a')
    s=step(project); s['argv'] += ['--data-dir',str(p.parent)]
    r=run_for(project,s); e=jobs.Engine()
    before=e._step_fingerprint(r,s)
    p.write_bytes(b'changed')
    assert before != e._step_fingerprint(r,s)


def test_invalid_output_encoding_finishes_cleanly(project):
    r=run_for(project,step(project,"import os; os.write(1, bytes([255,10]))"))
    jobs.Engine()._execute(r)
    assert r.status in ('done','failed')


def test_timeout_fails_cleanly(project,monkeypatch):
    monkeypatch.setenv('JOBS_STEP_TIMEOUT','0.1')
    s=step(project,'import time; time.sleep(10)')
    r=run_for(project,s)
    jobs.Engine()._execute(r)
    assert r.status=='failed' and s['returncode']==124


def test_cancellation_during_launch_stops_output(project,monkeypatch):
    spawned=threading.Event();release=threading.Event()
    original=jobs.subprocess.Popen
    marker=storage.project_dir(project['id'])/'output/after-cancel.txt'
    s=step(project,'import time,pathlib; time.sleep(0.25); pathlib.Path("output/after-cancel.txt").write_text("still running")')
    run=run_for(project,s);engine=jobs.Engine();engine.runs[run.run_id]=run
    def delayed_registration(*args,**kwargs):
        proc=original(*args,**kwargs)
        spawned.set();release.wait(timeout=2)
        return proc
    monkeypatch.setattr(jobs.subprocess,'Popen',delayed_registration)
    worker=threading.Thread(target=engine._execute,args=(run,))
    worker.start()
    try:
        assert spawned.wait(timeout=2)
        assert engine.cancel(run.run_id)
    finally:
        release.set();worker.join(timeout=3)
    assert not marker.exists(),{'status':run.status,'output':marker.read_text()}


def test_overlap_check_is_atomic(project,monkeypatch):
    engine=jobs.Engine()
    barrier=threading.Barrier(2)
    original=jobs.Run
    def gated_run(*a,**kw):
        time.sleep(0.05)
        return original(*a,**kw)
    monkeypatch.setattr(jobs,'Run',gated_run)
    monkeypatch.setattr(engine,'_execute',lambda r: None)
    def launch(_):
        barrier.wait(timeout=3)
        try:
            return engine.create_run(project['id'],None,'race',[step(project,outputs=['output/shared.wxf'])]).run_id
        except ValueError:
            return 'rejected'
    with ThreadPoolExecutor(2) as pool:
        outcomes=list(pool.map(launch,range(2)))
    assert outcomes.count('rejected')==1,outcomes


def test_concurrent_project_edits_preserved(project,monkeypatch):
    barrier=threading.Barrier(2)
    original=main._get_project
    def gated_get(pid):
        p=original(pid)
        time.sleep(0.05)
        return p
    monkeypatch.setattr(main,'_get_project',gated_get)
    with ThreadPoolExecutor(2) as pool:
        result=list(pool.map(lambda i: main.api_create_flow(project['id'], {'name':f'flow{i}'}),range(2)))
    saved=storage.load_project(project['id'])['flows']
    assert len(saved)==len(result)==2,([f['name'] for f in saved], [f['name'] for f in result])


def test_sse_delivers_after_full_buffer(project,monkeypatch):
    async def check():
        r=run_for(project,step(project)); r.status='running'
        for i in range(jobs.MAX_BUFFERED_EVENTS): r.emit('log',{'line':str(i)})
        monkeypatch.setattr(main.engine,'get',lambda _:r)
        stream=main.api_run_events(r.run_id).body_iterator
        for i in range(jobs.MAX_BUFFERED_EVENTS): await stream.__anext__()  # anext() needs Python 3.10
        r.emit('log',{'line':'AFTER BUFFER'})
        r.status='done'; r.emit('status',{'status':'done'}); r.emit('end',{})
        tail=[]
        async for item in stream: tail.append(item)
        assert 'AFTER BUFFER' in ''.join(tail) and 'event: end' in ''.join(tail),tail
    asyncio.run(check())


def test_restarted_active_run_reconciled(project,client):
    r=run_for(project,step(project)); r.status='running'
    jobs.Engine()._persist(r)
    res=client.get('/api/runs/'+r.run_id)
    assert res.status_code==200
    assert res.json()['status'] not in ['running','queued'],res.json()


@pytest.mark.skipif(os.name == 'nt', reason='Bash export execution is supported through WSL on Windows')
def test_export_dry_run(project,matrix_setup):
    result=compiler.export_flow_script(project,matrix_graph(),'audit-export')
    assert result['ok'],result
    r=subprocess.run(['bash',result['path'],'--dry-run'],capture_output=True,text=True,timeout=10)
    assert r.returncode==0,r.stdout+r.stderr


@pytest.mark.skipif(os.name == 'nt', reason='Bash export execution is supported through WSL on Windows')
def test_export_wolfram_bash_syntax(project,monkeypatch):
    project['alphabets']=[dict(id='a',name='A',letters=['x','y'],variables=['x','y'],expressions=['x','y'],roots={},
        properties=[dict(id='p',type='integrability',params={},status='pending')])]
    monkeypatch.setattr(compiler,'find_wolframscript',lambda:'/bin/true')
    g=dict(nodes=[dict(id='a',type='alphabet',data=dict(alphabet_id='a',selected_properties=['p']))],edges=[])
    result=compiler.export_flow_script(project,g,'wolfram-audit')
    assert result['ok'],result
    r=subprocess.run(['bash','-n',result['path']],capture_output=True,text=True,timeout=10)
    assert r.returncode==0,r.stderr


@pytest.mark.skipif(os.name == 'nt', reason='Bash export execution is supported through WSL on Windows')
def test_export_binary_rebuild_invalidates_cache(project,matrix_setup):
    result=compiler.export_flow_script(project,matrix_graph(),'audit-export')
    assert result['ok'],result
    first=subprocess.run(['bash',result['path']],capture_output=True,text=True,timeout=10)
    assert first.returncode==0,first.stdout+first.stderr
    matrix_setup.write_text('#!/bin/sh\nprintf v2 > "$4"\n')
    second=subprocess.run(['bash',result['path']],capture_output=True,text=True,timeout=10)
    assert second.returncode==0,second.stdout+second.stderr
    assert (storage.project_dir(project['id'])/'output/square.wxf').read_bytes()==b'v2',second.stdout


def test_changed_alphabet_invalidates_ready_properties(project,client):
    a=client.post(f'/api/projects/{project["id"]}/alphabets',json={'name':'A','letters':['x','y'],'expressions':['x','y']}).json()
    p=storage.load_project(project['id'])
    p['alphabets'][0]['properties']=[dict(id='prop',type='integrability',params={},status='ready',tensor_file='data/D.wxf',summary={'dims':[2,2,1]})]
    (storage.project_dir(project['id'])/'data/D.wxf').write_bytes(b'old')
    storage.save_project(p)
    r=client.put(f'/api/projects/{project["id"]}/alphabets/{a["id"]}',json={'expressions':['x','1-x']})
    assert r.status_code==200
    assert r.json()['properties'][0]['status']!='ready',r.json()


def test_flow_revision_conflict_preserves_newer_edit(project, client):
    base = f'/api/projects/{project["id"]}/flows'
    f = client.post(base, json={'name': 'initial'}).json()
    url = base + '/' + f['id']
    newer = client.put(url, json={'revision': 0, 'name': 'newer'})
    stale = client.put(url, json={'revision': 0, 'name': 'stale'})
    assert newer.status_code == 200 and newer.json()['revision'] == 1
    assert stale.status_code == 409
    assert storage.load_project(project['id'])['flows'][0]['name'] == 'newer'


def test_atomic_output_rollback_on_failure(project):
    p = storage.project_dir(project['id']) / 'output/result.wxf'
    p.write_bytes(b'previous valid result')
    s = step(project, "from pathlib import Path; Path('output/result.wxf').write_bytes(b'partial'); raise SystemExit(1)", ['output/result.wxf'])
    r = run_for(project, s)
    jobs.Engine()._execute(r)
    assert r.status == 'failed'
    assert p.read_bytes() == b'previous valid result'


def test_output_not_visible_until_success(project, monkeypatch):
    p = storage.project_dir(project['id']) / 'output/result.wxf'
    p.write_bytes(b'old')
    s = step(project, "from pathlib import Path; Path('output/result.wxf').write_bytes(b'new')", ['output/result.wxf'])
    r = run_for(project, s); e = jobs.Engine()
    original = e._run_step
    def launch(run, attempt):
        rc = original(run, attempt)
        assert p.read_bytes() == b'old'
        return rc
    monkeypatch.setattr(e, '_run_step', launch)
    e._execute(r)
    assert r.status == 'done' and p.read_bytes() == b'new'


def test_same_size_input_and_output_damage_invalidates_cache(project):
    root = storage.project_dir(project['id'])
    src = root / 'data/input.wxf'; src.write_bytes(b'aaa')
    out = root / 'output/result.wxf'; out.write_bytes(b'bbb')
    s = step(project, outputs=['output/result.wxf']); s['argv'].append(str(src))
    r = run_for(project, s); e = jobs.Engine(); e._record_sigs(r, s)
    assert e._step_cached(r, s)
    stamp = src.stat(); src.write_bytes(b'ccc'); os.utime(src, ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
    assert not e._step_cached(r, s)
    e._record_sigs(r, s); out.write_bytes(b'bad')
    assert not e._step_cached(r, s)


def test_sse_reconnect_cursor(project, monkeypatch):
    from starlette.requests import Request
    async def check():
        r = run_for(project, step(project))
        r.emit('log', {'line':'already seen'}); r.emit('log', {'line':'new'})
        r.status='done'; r.emit('end', {})
        monkeypatch.setattr(main.engine, 'get', lambda _: r)
        req = Request({'type':'http','headers':[(b'last-event-id', b'0')]})
        text = ''.join([item async for item in main.api_run_events(r.run_id, req).body_iterator])
        assert 'already seen' not in text and 'new' in text and 'event: end' in text
    asyncio.run(check())


def test_precomputed_property_cannot_follow_changed_alphabet(project, client):
    p=storage.load_project(project['id'])
    p['alphabets']=[dict(id='a',name='A',letters=['x'],expressions=['x'],properties=[dict(id='p',type='precomputed_tensor',status='ready',tensor_file='data/P.wxf')])]
    (storage.project_dir(project['id'])/'data/P.wxf').write_bytes(b'precomputed')
    storage.save_project(p)
    client.put(f'/api/projects/{project["id"]}/alphabets/a',json={'expressions':['1-x']})
    p=storage.load_project(project['id'])
    result=compiler.compile_flow(p,{'nodes':[dict(id='n',type='alphabet',data={'alphabet_id':'a','selected_properties':['p']})],'edges':[]})
    assert not result['ok']


def test_changed_input_during_run_rejects_publication(project, monkeypatch):
    root=storage.project_dir(project['id'])
    src=root/'data/input.wxf'; src.write_bytes(b'original')
    out=root/'output/result.wxf'; out.write_bytes(b'previous')
    s=step(project,"from pathlib import Path; Path('output/result.wxf').write_bytes(b'new')",['output/result.wxf'])
    s['argv'].append(str(src))
    r=run_for(project,s); e=jobs.Engine(); original=e._run_step
    def launch(run, attempt):
        rc=original(run,attempt); src.write_bytes(b'changed'); return rc
    monkeypatch.setattr(e,'_run_step',launch)
    e._execute(r)
    assert r.status=='failed' and out.read_bytes()==b'previous'


def test_cancelled_job_does_not_release_reservation_before_exit(project, monkeypatch):
    engine=jobs.Engine(); started=threading.Event(); release=threading.Event()
    def worker(run): started.set(); release.wait(timeout=2)
    monkeypatch.setattr(engine,'_execute',worker)
    r=engine.create_run(project['id'],None,'first',[step(project,outputs=['output/shared.wxf'])])
    assert started.wait(timeout=2)
    try:
        assert engine.cancel(r.run_id)
        with pytest.raises(ValueError):
            engine.create_run(project['id'],None,'second',[step(project,outputs=['output/shared.wxf'])])
    finally: release.set()


def test_long_exact_coefficient_is_not_treated_as_a_filename(project):
    s=step(project); s['argv'].append('1/' + '9' * 1000)
    r=run_for(project,s)
    assert jobs.Engine()._step_fingerprint(r,s)
