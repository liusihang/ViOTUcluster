#!/usr/bin/env python3

import os
import sys
import unittest
from concurrent.futures import Future

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ViOTUcluster import task_utils


def done_future(result=None, exc=None):
    future = Future()
    if exc is not None:
        future.set_exception(exc)
    else:
        future.set_result(result)
    return future


class TestTaskUtilsContract(unittest.TestCase):
    def test_ensure_futures_succeeded_returns_when_all_ok(self):
        task_utils.ensure_futures_succeeded(
            [done_future(result=1), done_future(result=2)],
            context="unit-test",
        )

    def test_ensure_futures_succeeded_raises_aggregated_error(self):
        with self.assertRaises(RuntimeError) as ctx:
            task_utils.ensure_futures_succeeded(
                [done_future(result=1), done_future(exc=ValueError("boom"))],
                context="unit-test",
            )
        self.assertIn("unit-test", str(ctx.exception))
        self.assertIn("boom", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
