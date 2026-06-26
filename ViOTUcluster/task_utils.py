#!/usr/bin/env python3
"""
Helpers for orchestration tasks that need consistent failure propagation.
"""

from concurrent.futures import as_completed


def ensure_futures_succeeded(futures, context: str) -> None:
    """Wait for futures and raise a single error if any of them failed."""
    errors = []
    for future in as_completed(list(futures)):
        try:
            future.result()
        except Exception as exc:
            errors.append(str(exc))

    if errors:
        raise RuntimeError(f"{context} failed: {'; '.join(errors)}")
