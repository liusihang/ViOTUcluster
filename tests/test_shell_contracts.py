#!/usr/bin/env python3

import os
import unittest


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_text(*parts):
    with open(os.path.join(REPO_ROOT, *parts), "r", encoding="utf-8") as handle:
        return handle.read()


class TestShellContracts(unittest.TestCase):
    def test_drep_module_uses_min_length_for_bins(self):
        content = read_text("Modules", "drep_module.sh")
        self.assertIn('-l "${MIN_LENGTH}"', content)

    def test_dependency_check_covers_runtime_tools(self):
        content = read_text("Modules", "ViOTUcluster_Check")
        for tool_name in ("conda", "viralverify", "sambamba", "coverm", "parallel", "makeblastdb", "blastn"):
            self.assertIn(f'"{tool_name}"', content)


if __name__ == "__main__":
    unittest.main()
