#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster.compat import fastaparser as fastaparser_compat


class TestFastaparserCompat(unittest.TestCase):
    def test_read_fasta_yields_header_sequence_tuples(self):
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            handle.write(">seq1 desc\nACGT\nTGCA\n>seq2\nNNNN\n")
            fasta_path = handle.name

        try:
            records = list(fastaparser_compat.read_fasta(fasta_path))
        finally:
            os.unlink(fasta_path)

        self.assertEqual(
            [
                (">seq1 desc", "ACGTTGCA"),
                (">seq2", "NNNN"),
            ],
            records,
        )
