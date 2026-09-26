#!/usr/bin/env python3

import os
from contextlib import ExitStack
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

from ViOTUcluster import viralprediction


REPO_ROOT = Path(__file__).resolve().parents[1]
DATABASE_INSTALLER = REPO_ROOT / "Modules" / "ViOTUcluster_download-database"


class TestVirSorterDatabaseSetup(unittest.TestCase):
    def test_failed_setup_stops_before_hmm_conversion(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            stub_dir = root / "bin"
            stub_dir.mkdir()
            call_log = root / "calls.log"
            virsorter = stub_dir / "virsorter"
            virsorter.write_text(
                "#!/bin/bash\n"
                'printf "virsorter %s\\n" "$*" >> "$ISSUE32_CALL_LOG"\n'
                "exit 42\n",
                encoding="utf-8",
            )
            virsorter.chmod(0o755)
            hmmconvert = stub_dir / "hmmconvert"
            hmmconvert.write_text(
                "#!/bin/bash\n"
                'printf "hmmconvert %s\\n" "$*" >> "$ISSUE32_CALL_LOG"\n',
                encoding="utf-8",
            )
            hmmconvert.chmod(0o755)

            # A stale prefix with directories but no database files is not complete.
            for part in ("group", "hmm", "rbs"):
                (root / "db" / part).mkdir(parents=True)
            environment = os.environ.copy()
            environment["PATH"] = f"{stub_dir}:{environment['PATH']}"
            environment["ISSUE32_CALL_LOG"] = str(call_log)
            result = subprocess.run(
                ["bash", str(DATABASE_INSTALLER), str(root), "1"],
                env=environment,
                capture_output=True,
                text=True,
            )

            self.assertNotEqual(0, result.returncode)
            self.assertIn("VirSorter2 setup failed", result.stdout + result.stderr)
            self.assertEqual(
                [f"virsorter setup --skip-deps-install -d {root / 'db'} -j 1"],
                call_log.read_text(encoding="utf-8").splitlines(),
            )
            self.assertNotIn("Database setup completed", result.stdout)


class TestVirSorterPrediction(unittest.TestCase):
    def test_prediction_uses_main_environment(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            sample = Path(temp_dir) / "sample.fa"
            sample.write_text(">sample\nACGT\n", encoding="utf-8")
            with ExitStack() as stack:
                stack.enter_context(mock.patch.object(viralprediction, "OUTPUT_DIR", temp_dir))
                stack.enter_context(mock.patch.object(viralprediction, "DATABASE", temp_dir))
                stack.enter_context(mock.patch.object(viralprediction, "Group", "dsDNAphage"))
                stack.enter_context(mock.patch.object(viralprediction, "CONCENTRATION_TYPE", "non-concentration"))
                stack.enter_context(mock.patch.object(viralprediction, "THREADS", 1))
                stack.enter_context(mock.patch.object(viralprediction, "resolve_viralverify_command", return_value=["viralverify"]))
                run_command = stack.enter_context(mock.patch.object(viralprediction, "run_command", return_value=""))
                viralprediction.process_file(str(sample))

            virsorter_calls = [
                call.args[0]
                for call in run_command.call_args_list
                if call.args[0][:2] == ["virsorter", "run"]
            ]
            self.assertEqual(1, len(virsorter_calls))
            self.assertIn("--use-conda-off", virsorter_calls[0])


if __name__ == "__main__":
    unittest.main()
