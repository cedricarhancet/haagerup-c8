"""Axis D: dimension of the space of s-submodules U_2 of V'_2 with the forced decomposition and containing the
fixed part (L_s)_2 + C omega (omega the conformal vector of V'), counted as a product of Grassmannians."""
import json
from collections import Counter
d = json.load(open('u2_embedding_test.json'))
def parse(x):
    return {tuple(int(t) for t in k.strip('()').split(',')): v for k, v in x.items()}
for name, rec in d.items():
    V2 = parse(rec['V2']); U2 = parse(rec['solutions'][0])
    k1 = 1 if name.startswith('A11') else 3
    fixed = Counter({(2, 0): 1, (0, 2): 1, (0, 4): 1, (2, 2): 1, (0, 0): 3})
    if k1 >= 2: fixed[(4, 0)] += 1
    tot = 0; rows = []
    for t in sorted(U2):
        u, v, l = U2[t], V2[t], fixed.get(t, 0)
        g = (u - l)*(v - u)
        tot += g; rows.append((t, u, v, l, g))
    print(name, ': total Grassmannian dimension', tot)
    for r in rows: print('    type %s: mult in U_2 %d, in V\'_2 %d, fixed %d -> %d parameters' % r)
