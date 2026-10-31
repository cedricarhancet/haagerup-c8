#!/usr/bin/env python3
"""Recheck the 39 stored rational polynomial systems without the VOA engine."""
import json
from pathlib import Path
import sympy as sp

path = Path(__file__).with_name('primary33-exclusions.json')
data = json.loads(path.read_text())
assert data['coefficient_field'] == 'Q'
assert data['primary_type'] == [3, 3]
assert sorted(map(int, data['candidates'])) == [0,1,2,3,4,5,6,7,8,11,12,13,14]
count = 0
for candidate in data['candidates'].values():
    assert [chart['first_nonzero_coordinate'] for chart in candidate['charts']] == [0,1,2]
    for chart in candidate['charts']:
        variables = [sp.Symbol(name) for name in chart['variables']]
        polys = [sp.sympify(text) for text in chart['polynomials']]
        assert list(sp.groebner(polys, *variables, order='grevlex', domain=sp.QQ)) == [1]
        count += 1
assert count == data['verified_projective_charts'] == 39
print(json.dumps({'unit_ideals_recomputed': count, 'VOA_engine_used': False,
                  'coefficient_field': 'Q', 'result': 'PASS'}, indent=2))
