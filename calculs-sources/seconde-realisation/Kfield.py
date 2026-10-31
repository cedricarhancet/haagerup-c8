"""Exact arithmetic in K = Q(r, q), r^2 = 13, q^2 = -3 (r = sqrt13 > 0, q = sqrt(-3) = i sqrt3), basis 1, r, q, rq.
Complex conjugation: q -> -q (r fixed).  Also: sqrt(6) sqrt(-7 + sqrt13) = q (r - 1) (principal branches)."""
from gmpy2 import mpq
import sympy as sp
Z = (mpq(0),)*4
ONE = (mpq(1), mpq(0), mpq(0), mpq(0))
R = (mpq(0), mpq(1), mpq(0), mpq(0))
Q = (mpq(0), mpq(0), mpq(1), mpq(0))
RQ = (mpq(0), mpq(0), mpq(0), mpq(1))
# products of basis elements: index pairs -> (coefficient, index)
_T = {(0, 0): (1, 0), (0, 1): (1, 1), (0, 2): (1, 2), (0, 3): (1, 3),
      (1, 1): (13, 0), (1, 2): (1, 3), (1, 3): (13, 2),
      (2, 2): (-3, 0), (2, 3): (-3, 1), (3, 3): (-39, 0)}
def mul(a, b):
    out = [mpq(0)]*4
    for i in range(4):
        if not a[i]: continue
        for j in range(4):
            if not b[j]: continue
            c, k = _T[(i, j) if i <= j else (j, i)]
            out[k] += c*a[i]*b[j]
    return tuple(out)
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def scal(c, a): return tuple(c*x for x in a)
def conj(a): return (a[0], a[1], -a[2], -a[3])
def iszero(a): return not any(a)
def from_sympy(expr):
    r, q = sp.symbols('rr qq')
    e = sp.sympify(expr)
    e = e.subs(sp.sqrt(6)*sp.sqrt(-7 + sp.sqrt(13)), q*(r - 1))
    e = e.subs(sp.sqrt(39), -sp.I*r*q).subs(sp.sqrt(3), -sp.I*q).subs(sp.sqrt(13), r)
    e = sp.expand(e)
    assert not e.has(sp.I), e
    P = sp.Poly(e, r, q)
    out = [mpq(0)]*4
    for (a, b), c in P.terms():
        c = sp.Rational(c)
        # r^a q^b reduce
        coef = mpq(int(c.p), int(c.q)) * (13**(a//2)) * ((-3)**(b//2))
        idx = (a % 2) + 2*(b % 2)       # 0:1, 1:r, 2:q, 3:rq
        out[idx] += coef
    return tuple(out)
def to_complex(a):
    import math
    r = math.sqrt(13); q = 1j*math.sqrt(3)
    return float(a[0]) + float(a[1])*r + float(a[2])*q + float(a[3])*r*q
