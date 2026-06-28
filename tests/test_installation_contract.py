#!/usr/bin/env python3

import os
import unittest


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def read_text(*parts):
    with open(os.path.join(REPO_ROOT, *parts), "r", encoding="utf-8") as handle:
        return handle.read()


class TestInstallationContract(unittest.TestCase):
    def test_main_environment_manifest_exists_and_carries_default_pipeline(self):
        manifest_path = os.path.join(REPO_ROOT, "environment.yml")
        self.assertTrue(os.path.exists(manifest_path), "environment.yml should exist at repo root")

        content = read_text("environment.yml")
        self.assertIn("name: ViOTUcluster", content)
        self.assertIn("- conda-forge", content)
        self.assertIn("- bioconda", content)

        for package_name in (
            "python=3.10",
            "fastp",
            "megahit",
            "spades",
            "viralverify",
            "virsorter",
            "genomad",
            "checkv",
            "drep",
            "vrhyme",
            "bwa",
            "sambamba",
            "coverm",
            "blast",
            "parallel",
        ):
            self.assertIn(package_name, content)

        self.assertNotIn("iphop", content, "iPhop should stay out of the main environment")

    def test_optional_satellite_manifests_exist(self):
        dram_manifest = read_text("environments", "dram.yml")
        iphop_manifest = read_text("environments", "iphop.yml")

        self.assertIn("name: DRAM", dram_manifest)
        self.assertIn("dram", dram_manifest)

        self.assertIn("name: iPhop", iphop_manifest)
        self.assertIn("python=3.8", iphop_manifest)
        self.assertIn("iphop", iphop_manifest)

    def test_setup_script_uses_conda_manifests_instead_of_prepacked_tarballs(self):
        content = read_text("setup_ViOTUcluster.sh")
        self.assertIn("environment.yml", content)
        self.assertNotIn("ViOTUcluster.tar.gz", content)
        self.assertNotIn("vRhyme.tar.gz", content)
        self.assertNotIn("zenodo.org", content)

    def test_dependency_check_matches_revised_runtime_contract(self):
        content = read_text("Modules", "ViOTUcluster_Check")
        self.assertIn('"metaspades.py"', content)
        self.assertIn('"vRhyme"', content)
        self.assertNotIn('"checkm"', content)


if __name__ == "__main__":
    unittest.main()
