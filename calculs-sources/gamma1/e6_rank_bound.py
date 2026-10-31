"""Exact check: no subspace of dimension >= 5 of a Cartan subalgebra of e6 makes the multiset of
weights of the 27-dimensional representation symmetric (W|h' = -W|h').  For a 5-dimensional h' with
annihilator C v in h*, symmetry forces w + w' in C v for some weights w, w'; we test all such lines
exactly (rational arithmetic).  Hence a reductive subalgebra t of e6 with 27 = 27* as t-modules has rank <= 4.
Usage: python3 e6_rank_bound.py
"""
from fractions import Fraction as Fr
import itertools
C = [[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]
def orbit(w0):
    seen = {tuple(w0)}; st = [tuple(w0)]
    while st:
        w = st.pop()
        for i in range(6):
            if w[i] != 0:
                v = tuple(w[k] - w[i]*C[i][k] for k in range(6))
                if v not in seen: seen.add(v); st.append(v)
    return sorted(seen)
W = orbit((1,0,0,0,0,0)); assert len(W) == 27
assert set(W) != {tuple(-x for x in w) for w in W}
# work in fundamental-weight coordinates; a class modulo C v is detected by reducing along v exactly
def key(w, v):
    # choose pivot index p with v[p] != 0 and subtract (w[p]/v[p]) v
    p = next(i for i in range(6) if v[i] != 0)
    t = Fr(w[p], v[p])
    return tuple(Fr(w[i]) - t*v[i] for i in range(6))
lines = set()
for a, b in itertools.product(W, W):
    s = tuple(x + y for x, y in zip(a, b))
    if any(s):
        from math import gcd
        g = 0
        for x in s: g = gcd(g, abs(x))
        s = tuple(x // g for x in s)
        if next(x for x in s if x != 0) < 0: s = tuple(-x for x in s)
        lines.add(s)
good = []
for v in lines:
    A = sorted(key(w, v) for w in W); B = sorted(key(tuple(-x for x in w), v) for w in W)
    if A == B: good.append(v)
print('candidate lines tested:', len(lines), ' lines giving a symmetric multiset:', good)
assert not good
print('=> rank of a reductive subalgebra of e6 with self-dual 27 is at most 4')
