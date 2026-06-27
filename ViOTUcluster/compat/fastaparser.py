#!/usr/bin/env python3
"""Minimal compat shim for viralverify/check_circular expectations."""


def read_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield (header, "".join(seq_parts))
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield (header, "".join(seq_parts))
