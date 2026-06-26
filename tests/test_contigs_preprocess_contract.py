#!/usr/bin/env python3

import os
import sys
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


if __name__ == "__main__":
    unittest.main()
