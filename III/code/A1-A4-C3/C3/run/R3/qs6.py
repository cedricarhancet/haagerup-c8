"""Exact arithmetic in Q(s), s = sqrt(-6), and W-vectors with coefficients in Q(s) (pairs of rational dicts).
An element is a pair (a, b) of mpq meaning a + b s ; complex conjugation: (a, -b)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gmpy2 import mpq
import e8e8_voa as V
D6 = -6
Z = (mpq(0), mpq(0)); ONE = (mpq(1), mpq(0)); S = (mpq(0), mpq(1))
def add(x, y): return (x[0] + y[0], x[1] + y[1])
def sub(x, y): return (x[0] - y[0], x[1] - y[1])
def mul(x, y): return (x[0]*y[0] + D6*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def conj(x): return (x[0], -x[1])
def norm(x): return x[0]*x[0] - D6*x[1]*x[1]
def inv(x):
    n = norm(x); return (x[0]/n, -x[1]/n)
def div(x, y): return mul(x, inv(y))
def scal(q, x): return (q*x[0], q*x[1])
def isreal(x): return x[1] == 0
def iszero(x): return x[0] == 0 and x[1] == 0
def tofloat(x): return complex(float(x[0]), float(x[1])*6**0.5)
def s(x): return '%s + %s*sqrt(-6)' % (x[0], x[1])
# vectors: (re, im) dicts
def vzero(): return ({}, {})
def vreal(v): return (dict(v), {})
def vadd(acc, v, c=ONE):
    a, b = c
    if a: V.vadd(acc[0], v[0], a); V.vadd(acc[1], v[1], a)
    if b: V.vadd(acc[0], v[1], D6*b); V.vadd(acc[1], v[0], b)
    return acc
def vclean(v): return (V.clean(v[0]), V.clean(v[1]))
def vconj(v): return (dict(v[0]), {k: -c for k, c in v[1].items()})
def vscale(v, c): return vclean(vadd(vzero(), v, c))
def mode(u, n, v):
    """u_(n) v for Q(s)-vectors"""
    rr = V.mode_vec(u[0], n, v[0]); ii = V.mode_vec(u[1], n, v[1])
    ri = V.mode_vec(u[0], n, v[1]); ir = V.mode_vec(u[1], n, v[0])
    re = V.vadd(rr, ii, D6); im = V.vadd(ri, ir)
    return (V.clean(re), V.clean(im))
def mode_rq(x, n, v):
    """x rational dict, v Q(s)-vector"""
    return (V.clean(V.mode_vec(x, n, v[0])), V.clean(V.mode_vec(x, n, v[1])))
def vac(v): return (v[0].get(V.VAC, mpq(0)), v[1].get(V.VAC, mpq(0)))
def iszero_v(v): return not V.clean(v[0]) and not V.clean(v[1])
def eq(u, v): return iszero_v(vadd(vadd(vzero(), u), v, (mpq(-1), mpq(0))))
def coeff(v, key): return (v[0].get(key, mpq(0)), v[1].get(key, mpq(0)))
