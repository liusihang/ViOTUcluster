#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster import viotucluster_allinone_cli


class TestAllInOneContract(unittest.TestCase):
    def test_preprocessing_outputs_ready_requires_fasta_payload(self):
        with tempfile.TemporaryDirectory() as output_dir:
            contigs_dir = os.path.join(output_dir, "Contigs")
            os.makedirs(contigs_dir, exist_ok=True)
            self.assertFalse(
                viotucluster_allinone_cli.preprocessing_outputs_ready(output_dir)
            )

    def test_preprocessing_outputs_ready_accepts_fasta_payload(self):
        with tempfile.TemporaryDirectory() as output_dir:
            contigs_dir = os.path.join(output_dir, "Contigs")
            os.makedirs(contigs_dir, exist_ok=True)
            with open(os.path.join(contigs_dir, "sample1.fasta"), "w", encoding="utf-8") as handle:
                handle.write(">s1\nACGT\n")
            self.assertTrue(
                viotucluster_allinone_cli.preprocessing_outputs_ready(output_dir)
            )


if __name__ == "__main__":
    unittest.main()
