#!/usr/bin/env python3

import os
import signal
import sys
import tempfile
import textwrap
import time
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster.pipeline import ViOTUclusterPipeline


class TestPipelineModuleTimeoutContract(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = self.tempdir.name
        self.input_dir = os.path.join(self.root, "input")
        self.raw_dir = os.path.join(self.root, "raw")
        self.output_dir = os.path.join(self.root, "output")
        self.db_dir = os.path.join(self.root, "db")
        self.script_dir = os.path.join(self.root, "shell")
        for path in (self.input_dir, self.raw_dir, self.output_dir, self.db_dir, self.script_dir):
            os.makedirs(path, exist_ok=True)
        with open(os.path.join(self.script_dir, "viral_prediction_module.sh"), "w", encoding="utf-8") as handle:
            handle.write("#!/usr/bin/env bash\n")
        self.old_script_dir = os.environ.get("ScriptDir")
        os.environ["ScriptDir"] = self.script_dir

    def tearDown(self):
        if self.old_script_dir is None:
            os.environ.pop("ScriptDir", None)
        else:
            os.environ["ScriptDir"] = self.old_script_dir
        self.tempdir.cleanup()

    def test_run_module_honors_timeout(self):
        script_path = os.path.join(self.script_dir, "sleepy.sh")
        with open(script_path, "w", encoding="utf-8") as handle:
            handle.write(textwrap.dedent("""\
                #!/usr/bin/env bash
                sleep 2
            """))
        os.chmod(script_path, 0o755)

        pipeline = ViOTUclusterPipeline(
            input_dir=self.input_dir,
            raw_seq_dir=self.raw_dir,
            output_dir=self.output_dir,
            database=self.db_dir,
            module_timeout_seconds=1,
        )
        self.assertTrue(pipeline.setup_directories())
        self.assertFalse(pipeline._run_module("Sleepy", script_path, os.environ.copy()))

    def test_run_module_timeout_kills_child_process_group(self):
        child_pid_file = os.path.join(self.root, "child.pid")
        script_path = os.path.join(self.script_dir, "spawner.sh")
        with open(script_path, "w", encoding="utf-8") as handle:
            handle.write(textwrap.dedent(f"""\
                #!/usr/bin/env bash
                python3 - <<'PY'
import os
import signal
import time

pid = os.fork()
if pid == 0:
    os.setsid()
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    time.sleep(30)
else:
    with open("{child_pid_file}", "w", encoding="utf-8") as handle:
        handle.write(str(pid))
    time.sleep(30)
PY
            """))
        os.chmod(script_path, 0o755)

        pipeline = ViOTUclusterPipeline(
            input_dir=self.input_dir,
            raw_seq_dir=self.raw_dir,
            output_dir=self.output_dir,
            database=self.db_dir,
            module_timeout_seconds=1,
        )
        self.assertTrue(pipeline.setup_directories())
        self.assertFalse(pipeline._run_module("Spawner", script_path, os.environ.copy()))

        with open(child_pid_file, "r", encoding="utf-8") as handle:
            child_pid = int(handle.read().strip())

        deadline = time.time() + 3
        while time.time() < deadline:
            try:
                os.kill(child_pid, 0)
            except OSError:
                break
            time.sleep(0.1)
        else:
            self.fail("timed-out module left child process alive")


if __name__ == "__main__":
    unittest.main()
