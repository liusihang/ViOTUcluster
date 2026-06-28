#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster import ContigsPreprocess


class TestContigsPreprocessContract(unittest.TestCase):
    def test_default_assembly_threads_respect_concurrency(self):
        self.assertEqual(
            8,
            ContigsPreprocess.resolve_assembly_threads(32, 4, None),
        )

    def test_explicit_assembly_threads_override_default(self):
        self.assertEqual(
            5,
            ContigsPreprocess.resolve_assembly_threads(32, 4, 5),
        )

    def test_main_returns_failure_when_no_pairs_found(self):
        with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as output_dir:
            exit_code = ContigsPreprocess.main(
                ["-i", input_dir, "-o", output_dir, "-a", "megahit", "-c", "4"]
            )
        self.assertEqual(1, exit_code)

    def test_main_returns_failure_when_orphan_reads_exist(self):
        with tempfile.TemporaryDirectory() as input_dir, tempfile.TemporaryDirectory() as output_dir:
            with open(os.path.join(input_dir, "sample1_R1.fastq"), "w", encoding="utf-8") as handle:
                handle.write("@r1\nACGT\n+\n####\n")
            with open(os.path.join(input_dir, "sample1_R2.fastq"), "w", encoding="utf-8") as handle:
                handle.write("@r2\nTGCA\n+\n####\n")
            with open(os.path.join(input_dir, "sample2_R1.fastq"), "w", encoding="utf-8") as handle:
                handle.write("@r3\nNNNN\n+\n####\n")

            exit_code = ContigsPreprocess.main(
                ["-i", input_dir, "-o", output_dir, "-a", "megahit", "-c", "4"]
            )
        self.assertEqual(1, exit_code)

    def test_scan_fastq_pairs_reports_orphans(self):
        with tempfile.TemporaryDirectory() as input_dir:
            with open(os.path.join(input_dir, "sample1_R1.fastq"), "w", encoding="utf-8") as handle:
                handle.write("@r1\nACGT\n+\n####\n")
            with open(os.path.join(input_dir, "sample1_R2.fastq"), "w", encoding="utf-8") as handle:
                handle.write("@r2\nTGCA\n+\n####\n")
            with open(os.path.join(input_dir, "sample2_R1.fastq"), "w", encoding="utf-8") as handle:
                handle.write("@r3\nNNNN\n+\n####\n")

            pairs, orphans = ContigsPreprocess.scan_fastq_pairs(input_dir)

        self.assertEqual(1, len(pairs))
        self.assertEqual(1, len(orphans))

    def test_determine_exit_code_fails_when_assembly_failed(self):
        with tempfile.TemporaryDirectory() as output_dir:
            self.assertEqual(
                1,
                ContigsPreprocess.determine_exit_code([], ["sample1"], output_dir),
            )

    def test_determine_exit_code_succeeds_with_generated_contigs(self):
        with tempfile.TemporaryDirectory() as output_dir:
            contigs_dir = os.path.join(output_dir, "Contigs")
            os.makedirs(contigs_dir, exist_ok=True)
            with open(os.path.join(contigs_dir, "sample1.fasta"), "w", encoding="utf-8") as handle:
                handle.write(">s1\nACGT\n")
            self.assertEqual(
                0,
                ContigsPreprocess.determine_exit_code([], [], output_dir),
            )


if __name__ == "__main__":
    unittest.main()
