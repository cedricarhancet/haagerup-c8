#!/usr/bin/env python3
"""Check all 256 line choices against the explicitly certified forbidden pairs."""
import ast
import itertools
import json
from pathlib import Path

here = Path(__file__).resolve().parent
cert = json.loads((here / 'closure-certify-exact-intervals.json').read_text())
assert cert['arithmetic'] == 'recursive exact rational and radical Arb evaluation'
forbidden = set()
for c in cert['certificates']:
    assert c['annihilator'] and c['mode'] in (1, 2)
    forbidden.add((ast.literal_eval(c['tau']), c['line'], ast.literal_eval(c['sigma']), c['other_line']))
assert len(forbidden) == cert['certified'] == 62
assert len(cert['not_certified']) == 10
assert cert['rejected'] == 72
kinds = sorted({t for t, _, _, _ in forbidden} | {s for _, _, s, _ in forbidden})
assert kinds == [(0,4), (1,1), (1,3), (1,5), (1,7), (2,2), (2,6), (3,3)]
surviving = []
exclusions = []
for bits in itertools.product(range(2), repeat=len(kinds)):
    choice = dict(zip(kinds, bits))
    witness = next((pair for pair in sorted(forbidden) if choice[pair[0]] == pair[1] and choice[pair[2]] == pair[3]), None)
    if witness is None:
        surviving.append(list(bits))
    else:
        exclusions.append({'choice': list(bits), 'forbidden_pair': witness})
assert surviving == [[0,1,1,0,0,1,0,1], [1,0,0,1,1,0,1,0]]
assert len(exclusions) == 254
print(json.dumps({'types': kinds, 'number_choices': 256, 'excluded': len(exclusions), 'surviving': surviving, 'exclusion_witnesses': exclusions, 'result':'PASS'}, indent=2))
