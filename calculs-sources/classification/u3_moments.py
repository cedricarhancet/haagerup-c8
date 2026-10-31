"""Moments of U_3 (dimension 956) as an s-module from the exact relations of weights 2, 4, 6 with the extra
positions (0,1), (0,2) of the vacuum component (same method as u2_moments.py)."""
import json
from math import comb
import sympy as sp
import obstruction as O, vvmf_exact as V, u2_moments as UM, weight6_mixed as WM
def dfact(m):
    r = 1
    while m > 1:
        r *= m; m -= 2
    return r
names = V.names; d = V.d
A = [(0, -1), (0, 0)] + [(z, 0) for z in range(1, d)] + [(0, 1), (0, 2)]
PQ = UM.PQ
E2 = [1, -24, -72, -96]
def E2pow(k):
    s = [1, 0, 0, 0]
    for _ in range(k):
        s = [sum(s[i]*E2[j - i] for i in range(j + 1)) for j in range(4)]
    return s
def solve_U3(name):
    tab, k1, k2 = UM.candidate_tables(name)
    u2, _, _ = UM.solve_U2(name)
    nrm = lambda i, j: sp.Integer(2*k1)**sp.Rational(i, 2)*sp.Integer(2*k2)**sp.Rational(j, 2)
    def tr(z, i, j):
        if (i, j) == (0, 0): return UM.DIMS[z]
        if i % 2 or j % 2: return 0
        return sp.Rational(tab[z][(i, j)])/nrm(i, j)
    def trU(g, i, j, U3):
        if (i, j) == (0, 0): return [1, 6, 120, 956][g]
        if i % 2 or j % 2: return 0
        if g == 0: return 0
        if g == 1: return sp.Rational(tab['s'][(i, j)])/nrm(i, j)
        if g == 2: return sp.Rational(u2[(i, j)])/nrm(i, j)
        return U3[(i, j)]
    U3 = {pq: sp.Symbol('V_%d_%d' % pq) for pq in PQ}
    known = {}; res = {}
    for (p, q) in PQ:
        n = p + q
        R = O.relations(n, [0], A, 2)
        def coeff(trf):
            tot = 0
            for kk in range(n//2 + 1):
                ck = comb(n, 2*kk)*dfact(2*kk - 1)*sp.Rational(-1, 12)**kk
                m = n - 2*kk
                for a in range(kk + 1):
                    b = kk - a; i, j = p - 2*a, q - 2*b
                    if i < 0 or j < 0 or i + j != m: continue
                    tot += ck*comb(kk, a)*comb(m, i)*trf(i, j, kk)
            return tot
        w = {}
        for (z, mm) in A:
            if z == 0:
                g = mm + 1
                w[(z, mm)] = coeff(lambda i, j, kk, g=g: sum(E2pow(kk)[t]*trU(g - t, i, j, U3) for t in range(g + 1)))
            else:
                w[(z, mm)] = coeff(lambda i, j, kk, zn=names[z]: tr(zn, i, j))
        eqs = [sp.expand(sum(c*w[pos] for c, pos in zip(row, A)).subs(known)) for row in R]
        unknown = [u for u in U3.values() if any(e.has(u) for e in eqs)]
        sol = sp.solve(eqs, unknown, dict=True)
        if not sol:
            print(name, (p, q), 'INCONSISTENT'); res[(p, q)] = None; continue
        val = sol[0].get(U3[(p, q)])
        if val is not None: known[U3[(p, q)]] = val
        res[(p, q)] = sp.nsimplify(val*nrm(p, q)) if val is not None else None
        print('%s: Tr_{U_3} H\'^%d H\'\'^%d = %s' % (name, p, q, res[(p, q)]))
    return res, u2, k1, k2
if __name__ == '__main__':
    out = {}
    for name in WM.CANDS:
        r, u2, k1, k2 = solve_U3(name)
        out[name] = {'U3': {str(k): str(v) for k, v in r.items()}, 'U2': {str(k): str(v) for k, v in u2.items()}}
    json.dump(out, open('u3_moments.json', 'w'), indent=1)
