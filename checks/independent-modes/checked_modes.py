"""Opt-in guard for the historical norm<=8 engine; original files unchanged."""
from independent_modes import weight

def checked_mode(engine, u, n, v):
    if weight(u) > 4 or weight(v) > 4:
        raise ValueError("Historical engine input outside the audited weight window")
    target = weight(u) + weight(v) - n - 1
    if target > 4:
        raise ValueError("Historical engine output outside the audited weight window")
    return engine.mode(u, n, v)
