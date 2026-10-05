"""Exercise process-tree limits and pipeline acceptance with small owned jobs."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
MEASURE = ROOT / 'bench/measure_process.py'
PIPELINE = ROOT / 'bench/run_pipeline.py'


class BudgetTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def measure(self, script, *extra):
        command = [sys.executable, str(MEASURE), '--report', str(self.root/'report.json'),
                   '--log', str(self.root/'child.log'), '--status', str(self.root/'live.json'),
                   '--virtual-gib', '1', '--timeout', '5', *extra,
                   '--', sys.executable, '-c', script]
        result = subprocess.run(command, capture_output=True, text=True, timeout=15)
        data = json.loads((self.root/'report.json').read_text())
        self.assertEqual(json.loads((self.root/'live.json').read_text())['state'], 'finished')
        return result, data

    def test_success_and_status(self):
        result, data = self.measure('print("accepted")')
        self.assertEqual(result.returncode, 0)
        self.assertFalse(data['timeout'])
        self.assertFalse(data['rss_limit_exceeded'])
        self.assertGreaterEqual(data['cpu_user_seconds'], 0)

    def test_deadline(self):
        result, data = self.measure('import time; time.sleep(5)', '--timeout', '.2')
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(data['timeout'])
        self.assertLess(data['wall_seconds'], 3)

    def test_child_rss_limit(self):
        script = ('import subprocess,sys; p=subprocess.Popen([sys.executable,"-c",'
                  '"import time; a=bytearray(80*1024**2); time.sleep(5)"]); '
                  'print(p.pid,flush=True); p.wait()')
        result, data = self.measure(script, '--rss-gib', '.04')
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(data['rss_limit_exceeded'])
        self.assertGreaterEqual(data['max_processes'], 2)
        child = int((self.root/'child.log').read_text().splitlines()[0])
        status = Path(f'/proc/{child}/stat')
        if status.exists():
            self.assertIn(status.read_text().split(') ', 1)[1].split()[0], ('Z', 'X'))

    def test_ignoring_descendant_is_killed_after_leader_exits(self):
        worker = ('import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                  'a=bytearray(80*1024**2); time.sleep(30)')
        script = ('import subprocess,sys; p=subprocess.Popen([sys.executable,"-c",'
                  f'{worker!r}]); print(p.pid,flush=True); p.wait()')
        result, data = self.measure(script, '--rss-gib', '.04')
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(data['rss_limit_exceeded'])
        child = int((self.root/'child.log').read_text().splitlines()[0])
        status = Path(f'/proc/{child}/stat')
        if status.exists():
            self.assertIn(status.read_text().rsplit(') ', 1)[1].split()[0], ('Z', 'X'))

    def configuration(self, stages, seconds=10):
        config = dict(run_directory=str(self.root/'run'), working_directory=str(ROOT),
                      measure_program=str(MEASURE), deadline_unix_seconds=time.time()+seconds,
                      virtual_gib=1, rss_gib=.2, min_available_gib=0, stages=stages)
        path = self.root/'config.json'
        path.write_text(json.dumps(config))
        return path

    def pipeline(self, stages, seconds=10):
        result = subprocess.run([sys.executable, str(PIPELINE), str(self.configuration(stages, seconds))],
                                capture_output=True, text=True, timeout=20)
        return result, json.loads((self.root/'run/status.json').read_text())

    def test_stages_and_output_certificates(self):
        output = self.root/'accepted'
        result, data = self.pipeline([
            dict(name='one', command=[sys.executable, '-c',
                 f'from pathlib import Path; Path({str(output)!r}).write_text("exact"); print("PASS exact")'],
                 required_outputs=[str(output)], success_markers=['PASS exact']),
            dict(name='two', command=[sys.executable, '-c',
                 f'from pathlib import Path; assert Path({str(output)!r}).read_text()=="exact"; print("PASS readback")'],
                 success_markers=['PASS readback'])])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(data['state'], 'complete')
        self.assertEqual(len(data['stages']), 2)
        self.assertEqual(len(data['stages'][0]['outputs'][0]['sha256']), 64)

    def test_missing_certificate_stops_pipeline(self):
        result, data = self.pipeline([
            dict(name='one', command=[sys.executable, '-c', 'print("candidate")'],
                 success_markers=['PASS exact']),
            dict(name='two', command=[sys.executable, '-c', 'raise Exception("must not run")'])])
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(data['state'], 'failed')
        self.assertEqual(len(data['stages']), 1)

    def test_shared_deadline(self):
        result, data = self.pipeline([
            dict(name='one', command=[sys.executable, '-c', 'import time; time.sleep(.3)']),
            dict(name='two', command=[sys.executable, '-c', 'import time; time.sleep(5)'])], seconds=.8)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(data['state'], 'time_limit')
        self.assertLess(data['wall_seconds'], 2)

    def test_cancellation_reaches_child(self):
        config = self.configuration([dict(name='one', command=[sys.executable, '-c',
                                   'import time; time.sleep(10)'])])
        process = subprocess.Popen([sys.executable, str(PIPELINE), str(config)],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        live = self.root/'run/one-live.json'
        try:
            for _ in range(200):
                if live.exists():
                    break
                time.sleep(.01)
            self.assertTrue(live.exists())
            child = json.loads(live.read_text())['pid']
            process.send_signal(signal.SIGTERM)
            process.communicate(timeout=10)
            self.assertEqual(json.loads((self.root/'run/status.json').read_text())['state'], 'cancelled')
            self.assertFalse(Path(f'/proc/{child}').exists())
        finally:
            if process.poll() is None:
                process.terminate()
                process.communicate(timeout=10)


if __name__ == '__main__':
    unittest.main(verbosity=2)
