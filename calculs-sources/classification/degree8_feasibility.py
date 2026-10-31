"""Degree-eight test: existence of sl2-module structures of R_{d2}, R_{d3} with the power sums of degrees
2, 4, 6 (weights 2, 4, 6) and 8 (weight-eight relations, weightn_single_check.py), highest weights <= level.
"""
import json
from functools import lru_cache
import sympy as sp
import weightn_single_check as WN
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
@lru_cache(maxsize=None)
def feas(dim, t2, t4, t6, t8, nmax):
    if t2 == 0:
        return t4 == 0 and t6 == 0 and t8 == 0
    if dim <= 0 or nmax == 0 or not (t2 <= t4 <= t6 <= t8):
        return False
    for nn in range(min(nmax, dim - 1), 0, -1):
        if M(nn, 2) <= t2 and M(nn, 4) <= t4 and M(nn, 6) <= t6 and M(nn, 8) <= t8:
            if feas(dim - nn - 1, t2 - M(nn, 2), t4 - M(nn, 4), t6 - M(nn, 6), t8 - M(nn, 8), nn):
                return True
    return False
def wit(dim, t2, t4, t6, t8, nmax):
    if t2 == 0:
        return []
    for nn in range(min(nmax, dim - 1), 0, -1):
        if M(nn, 2) <= t2 and M(nn, 4) <= t4 and M(nn, 6) <= t6 and M(nn, 8) <= t8 and feas(dim - nn - 1, t2 - M(nn, 2), t4 - M(nn, 4), t6 - M(nn, 6), t8 - M(nn, 8), nn):
            return [nn] + wit(dim - nn - 1, t2 - M(nn, 2), t4 - M(nn, 4), t6 - M(nn, 6), t8 - M(nn, 8), nn)
out = {}
for key in [(1, 0), (3, 0), (9, 4), (35, 4)]:
    D, k, st, t2, t4, t6 = WN.ideal_structures(key)
    ok, eqs = WN.check(key, 8, verbose=False)
    sol = sp.solve(eqs, sp.symbols('U2 U3'), dict=True)[0]
    t8 = [sp.Rational(sol[sp.Symbol('U2')])*(2*k)**4, sp.Rational(sol[sp.Symbol('U3')])*(2*k)**4]
    res = []
    for jj, j in enumerate((1, 2)):
        if t8[jj].q != 1:
            res.append('d%d: Tr H^8 = %s not an integer' % (j + 1, t8[jj])); continue
        f = feas(WN.DIMS[j], t2[j], t4[j], t6[j], int(t8[jj]), min(k, WN.DIMS[j] - 1))
        res.append('d%d: (Tr H^2,4,6,8) = (%d,%d,%d,%d) %s %s' % (j + 1, t2[j], t4[j], t6[j], t8[jj], 'FEASIBLE' if f else 'INFEASIBLE',
                   wit(WN.DIMS[j], t2[j], t4[j], t6[j], int(t8[jj]), min(k, WN.DIMS[j] - 1)) if f else ''))
    print('ideal', key, 'level', k)
    for r in res: print('   ', r)
    out[str(key)] = res
json.dump(out, open('degree8_feasibility.json', 'w'), indent=1)
