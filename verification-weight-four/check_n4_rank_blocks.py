#!/usr/bin/env python3
"""Check the identity minors of the exported reduced N4 matrices.

This validates independence over the stated prime field. Membership in N4
is supplied by the reproducible generation script, not by this rank check.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def prime(n):
    if n < 2:
        return False
    d = 2
    while d*d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def check(directory):
    summary = json.loads((directory / "N4_rank_summary.json").read_text())
    matrices = np.load(directory / "N4_rank_blocks.npz", allow_pickle=False)
    names = json.loads((directory / "N4_rank_block_names.json").read_text())
    p = summary["p"]
    if not prime(p):
        raise ValueError("Modulus is not prime")
    if set(matrices.files) != set(names):
        raise ValueError("Incomplete block inventory")
    checked = {}
    total = 0
    for name, weight in names.items():
        matrix = matrices[name]
        if matrix.ndim != 2 or not np.issubdtype(matrix.dtype, np.integer):
            raise ValueError(f"Invalid matrix type: {name}")
        if np.any(matrix < 0) or np.any(matrix >= p):
            raise ValueError(f"Noncanonical field entries: {name}")
        pivots = []
        for row in matrix:
            indices = np.flatnonzero(row)
            if len(indices) == 0:
                raise ValueError(f"Zero certificate row: {name}")
            pivot = int(indices[0])
            if int(row[pivot]) != 1:
                raise ValueError(f"Pivot not normalized: {name}")
            pivots.append(pivot)
        if pivots != sorted(set(pivots)):
            raise ValueError(f"Invalid pivot order: {name}")
        # The pivot-column minor must be exactly the identity modulo p.
        minor = matrix[:, pivots]
        if not np.array_equal(minor, np.eye(len(pivots), dtype=matrix.dtype)):
            raise ValueError(f"Identity minor failed: {name}")
        rank = len(pivots)
        if rank != summary["block_ranks"][weight]:
            raise ValueError(f"Block rank differs from manifest: {name}")
        checked[weight] = rank
        total += rank
    if set(checked) != set(summary["block_ranks"]):
        raise ValueError("Manifest has unrepresented blocks")
    if total != summary["N4_rank"] or total != 6010:
        raise ValueError(f"Wrong total rank: {total}")
    return {"status": "PASS", "prime": p, "blocks": len(checked),
            "rank": total, "method": "identity minors in reduced coordinate matrices",
            "scope": "independence; generation in N4 is checked by the producer"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, nargs="?", default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    print(json.dumps(check(args.directory.resolve()), indent=2))
