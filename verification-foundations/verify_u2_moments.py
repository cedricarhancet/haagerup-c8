#!/usr/bin/env python3
"""Recompute the nine U_2 joint moments in each case and compare the saved data.
Requires sympy and mpmath, besides the supplied classification sources.
This does not rewrite any original output file.
"""
import argparse
import json
import os
import sys
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--sources', type=Path, default=Path(__file__).resolve().parents[1]/'calculs-sources'/'classification')
args = parser.parse_args()
sources = args.sources.resolve()
sys.path.insert(0, str(sources))
os.chdir(sources)
import u2_moments as UM
saved = json.loads(Path('u2_moments.json').read_text())
for name in UM.WM.CANDS:
    values, k1, k2 = UM.solve_U2(name)
    assert {str(k): str(v) for k, v in values.items()} == saved[name]
    print('MATCHES STORED EXACT MOMENTS:', name, flush=True)
print('All nine joint moments agree in both cases.', flush=True)
