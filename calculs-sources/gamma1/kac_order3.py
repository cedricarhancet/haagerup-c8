"""Order-3 elements of E8 (Kac coordinates) and the dimensions of their fixed subalgebras.
Kac: s_0 + sum a_i s_i = 3, marks a = (2,3,4,6,5,4,3,2) (Bourbaki labelling), s_i >= 0, gcd = 1.
Fixed subalgebra dimension = 8 + #{roots alpha = sum c_i alpha_i : sum c_i s_i = 0 mod 3}.
Usage: python3 kac_order3.py
"""
from itertools import product
import numpy as np
C = np.array([[2,0,-1,0,0,0,0,0],[0,2,0,-1,0,0,0,0],[-1,0,2,-1,0,0,0,0],[0,-1,-1,2,-1,0,0,0],
              [0,0,0,-1,2,-1,0,0],[0,0,0,0,-1,2,-1,0],[0,0,0,0,0,-1,2,-1],[0,0,0,0,0,0,-1,2]])
# positive roots in simple-root coordinates by the standard string algorithm
simple = [tuple(int(i == j) for j in range(8)) for i in range(8)]
pos = set(simple); frontier = list(simple)
while frontier:
    new = []
    for r in frontier:
        rv = np.array(r)
        for i in range(8):
            # alpha + alpha_i is a root iff <alpha, alpha_i^vee> < 0 or the string continues: use p - q rule
            s = tuple(rv + np.eye(8, dtype=int)[i])
            if s in pos: continue
            # p = how far we can go down from r along alpha_i
            p = 0; t = rv.copy()
            while True:
                t = t - np.eye(8, dtype=int)[i]
                if tuple(t) in pos: p += 1
                else: break
            q = p - int(rv @ C[:, i])
            if q > 0:
                pos.add(s); new.append(s)
    frontier = new
assert len(pos) == 120
marks = [2,3,4,6,5,4,3,2]
res = []
for s in product(range(4), repeat=8):
    tot = sum(a*x for a, x in zip(marks, s))
    if tot > 3: continue
    s0 = 3 - tot
    if np.gcd.reduce([s0] + list(s)) != 1: continue
    fixed = sum(1 for r in pos if sum(c*x for c, x in zip(r, s)) % 3 == 0)
    res.append(((s0,) + s, 8 + 2*fixed))
for k, dim in sorted(res, key=lambda t: t[1]):
    print('Kac coordinates (s0; s1..s8) =', k, ' dim of fixed subalgebra =', dim)
dims = [dim for _, dim in res]
assert sorted(dims) == [80, 86, 92, 134], dims  # s0 = 3 is the identity (gcd 3), excluded
print('unique class with fixed dimension 86:', [k for k, dim in res if dim == 86])
