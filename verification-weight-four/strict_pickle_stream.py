"""Fail-closed helpers for a future W4 certificate regeneration.

These helpers do not certify the mathematical contents of a record. In
particular, zero_labels must be produced by exact product recomputation, not
inferred from absent files or absent records. A manifest must distinguish a
certified empty block from a missing block before calling read_stream.
"""
from pathlib import Path
import pickle


def read_stream(path):
    """Read a required pickle stream, rejecting absence and truncated records."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Required certificate stream absent: {path}")
    records = []
    with path.open("rb") as stream:
        while True:
            offset = stream.tell()
            if stream.read(1) == b"":
                return records
            stream.seek(offset)
            try:
                records.append(pickle.load(stream))
            except (EOFError, pickle.UnpicklingError) as exc:
                raise ValueError(
                    f"Corrupt or truncated certificate at {path}:{offset}"
                ) from exc


def check_coverage(records, expected_labels, zero_labels):
    """Check exact label coverage, without silently discarding conflicting copies."""
    expected = set(expected_labels)
    zeros = set(zero_labels)
    nonzero = {}
    duplicates = 0
    for label, vector in records:
        if not vector:
            raise ValueError(f"Unexpected empty vector in nonzero stream: {label}")
        if label in nonzero:
            if nonzero[label] != vector:
                raise ValueError(f"Conflicting duplicate certificate: {label}")
            duplicates += 1
        else:
            nonzero[label] = vector
    if set(nonzero) & zeros:
        raise ValueError("A product is simultaneously marked zero and nonzero")
    observed = set(nonzero) | zeros
    if observed != expected:
        missing = sorted(expected - observed, key=repr)
        extra = sorted(observed - expected, key=repr)
        raise ValueError(f"Incomplete coverage: missing={missing}; extra={extra}")
    return {
        "expected": len(expected),
        "nonzero": len(nonzero),
        "zero": len(zeros),
        "identical_duplicate_records": duplicates,
    }
