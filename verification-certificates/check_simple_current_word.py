"""Verify an explicit word for the case-II simple current, without piv_w4 files."""
import argparse
import json
import os
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--source-dir', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
source = args.source_dir.resolve()
outpath = args.output.resolve()
sys.path.insert(0, str(source))
os.chdir(source)
import exact_engine as E
import Kfield as K
import lattice_voa as V
import cases as C
from gmpy2 import mpq

setup = E.setup(0)
u = setup['hw']['(1, 7)']
v = E.kv_from_terms(setup['hw']['(2, 6)'])
product = E.kprod(u, -1, v)
weight = (2, 2, 4, 5, 4, 2, 1, 1)
key = (weight, ())
coefficient = (mpq(0), mpq(507, 5), mpq(3549, 5), mpq(0))
assert product == {key: coefficient}
assert V.ip(weight, weight) == 8
assert C.weight_of_key(setup['T'], key) == (3, 13)
checks = []
for name in ('e1', 'e2'):
    value = E.kprod([(K.ONE, setup['T'][name])], 0, {key: K.ONE})
    checks.append({'operator': name + '(0)', 'zero': not value})
for i, current in enumerate(setup['U1']):
    for mode in range(1, 5):
        value = E.kprod([(K.ONE, current)], mode, {key: K.ONE})
        checks.append({'operator': f'y{i}({mode})', 'zero': not value})
assert len(checks) == 26 and all(check['zero'] for check in checks)
result = {'status': 'PASS', 'choice': setup['choice'],
          'word': 'u_(1,7)_{(-1)} u_(2,6)',
          'normalization': 'exact_engine.setup(0) highest-weight vectors',
          'lattice_vector': list(weight), 'norm_squared': 8,
          'Cartan_weight': [3, 13],
          'coefficient_basis': ['1', 'sqrt(13)', 'sqrt(-3)', 'sqrt(-39)'],
          'coefficient': [str(c) for c in coefficient],
          'singularity_checks': checks, 'uses_weight_four_pivot_files': False}
outpath.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': 'PASS', 'coefficient': result['coefficient'], 'singularity_checks': len(checks), 'output': str(outpath)}))
