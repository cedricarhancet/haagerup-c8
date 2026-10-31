#!/usr/bin/env python3
"""Check the rank of the stored 120 exact vectors, using only Python stdlib.

The generating words and their exact vectors are produced by i6_generation_certificate.py.
This independent checker certifies the modular rank, not the VOA product engine itself.
"""
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path

path = Path(__file__).with_name('i6_rank_certificate.json')
raw = path.read_bytes()
data = json.loads(raw)
p, r, q = data['prime'], data['sqrt13'], data['sqrt_minus3']
assert p == 8388571 and r*r % p == 13 and q*q % p == p-3
assert len(data['basis_vectors']) == len(data['words']) == data['rank'] == 120
for i, word in enumerate(data['words']):
    assert word[0] in ('seed', 'mode1')
    if word[0] == 'seed': assert 0 <= word[1] < len(data['seed_vectors'])
    else: assert len(word) == 3 and all(0 <= j < i for j in word[1:])
def residue(text):
    x = Fraction(text)
    assert x.denominator % p, 'Nonintegral coordinate'
    return x.numerator * pow(x.denominator, -1, p) % p
basis = {}
for vector in data['basis_vectors']:
    row = {}
    for key, coordinates in vector:
        state = ast.literal_eval(key)
        a,b,c,d = map(residue, coordinates)
        z = (a+b*r+c*q+d*r*q) % p
        if z: row[state] = z
    while row:
        pivot = min(row)
        if pivot not in basis:
            factor = pow(row[pivot], -1, p)
            basis[pivot] = {key:value*factor % p for key,value in row.items()}
            break
        factor = row[pivot]
        for key,value in basis[pivot].items():
            new = (row.get(key,0)-factor*value) % p
            if new: row[key] = new
            else: row.pop(key,None)
assert len(basis) == 120
print(json.dumps({'certificate_sha256':hashlib.sha256(raw).hexdigest(), 'prime':p,
                  'rank':len(basis), 'all_denominators_invertible':True,
                  'word_DAG_well_founded':True, 'result':'PASS'},indent=2))
