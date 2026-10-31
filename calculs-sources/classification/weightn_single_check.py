"""Higher-weight checks (n = 8, 10, 12, 14) for the four surviving simple ideals (Axis C, final step).

W_n(h) = sum_k C(n,2k) (2k-1)!! (-E2/12)^k <h,h>^k Tr h(0)^{n-2k} q^{L0-c/24}  (Heisenberg recursion) is a
holomorphic weight-n form for rho.  Leading coefficients: sum_k c_k Tr_{R_z} h^{n-2k}, c_k = C(n,2k)(2k-1)!!(-1/12)^k;
vacuum: w0 = c_{n/2}, w0p = sum_{k<n/2} c_k T_{n-2k} + c_{n/2}(6 - 12 n), T_m = Tr_s (ad h)^m.
For each ideal, R_{d1}, R_{d4}, R_{d5}, R_{d6} have uniquely determined module structures (weights 2-6), so their
traces of all degrees are known; R_{d2}, R_{d3} have known traces up to degree 6 and unknown degree-n traces.
Each exact relation of weight n then gives either a value for the unknowns or a CHECK.
"""
import json, sys
from math import comb
from fractions import Fraction as Fr
import sympy as sp
import obstruction as O, vvmf_exact as V, weight6_single as W6, sextic_joint as SJ

def dfact(m):
    r = 1
    while m > 1:
        r *= m; m -= 2
    return r
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
names = V.names; d = V.d
labels = ['w0', 'w0p'] + ['w_' + names[z] for z in range(1, d)]
A = [(0, -1), (0, 0)] + [(z, 0) for z in range(1, d)]
DIMS = [7, 42, 119, 5, 13, 14]

def structures_full(dim, t2, t4, t6, nmax):
    out = []
    def rec(rem, r2, r4, r6, mx, cur):
        if r2 == 0 and r4 == 0 and r6 == 0:
            out.append(tuple(cur + [0]*rem)); return
        for nn in range(min(mx, rem - 1), 0, -1):
            if M(nn, 2) <= r2 and M(nn, 4) <= r4 and M(nn, 6) <= r6:
                rec(rem - nn - 1, r2 - M(nn, 2), r4 - M(nn, 4), r6 - M(nn, 6), nn, cur + [nn])
    rec(dim, t2, t4, t6, nmax, [])
    return out

def ideal_structures(key):
    """unique structures of R_d1, R_d4, R_d5, R_d6 and degree <= 6 traces of R_d2, R_d3 (integer normalization)"""
    D, k, sol, Ss, t2, t4 = SJ.single_sextic(key)
    st = {}
    for j in (3, 4):
        s = SJ.structures(DIMS[j], (t2[j], t4[j]), min(k, DIMS[j] - 1)); assert len(s) == 1; st[j] = s[0]
    st[5] = tuple(D['R_d6'])
    sub = {Ss[j]: sp.Rational(sum(M(nn, 6) for nn in st[j]), (2*k)**3) for j in (3, 4, 5)}
    t6 = [int(sp.Rational(sol[Ss[j]].subs(sub))*(2*k)**3) if j < 3 else int(sub[Ss[j]]*(2*k)**3) for j in range(6)]
    s1 = structures_full(7, t2[0], t4[0], t6[0], min(k, 6)); assert len(s1) == 1; st[0] = s1[0]
    return D, k, st, t2, t4, t6

def check(key, n, verbose=True):
    D, k, st, t2, t4, t6 = ideal_structures(key)
    n2 = 2*k
    R = O.relations(n, [0], A, 1)
    c = [sp.Rational(comb(n, 2*kk)*dfact(2*kk - 1)) * sp.Rational(-1, 12)**kk for kk in range(n//2 + 1)]
    ps = {m: W6.sector_power_sums(D['lab'], key[1], m) for m in range(2, n + 1, 2)}
    T = {m: sp.Rational(2*2**m, n2**(m//2)) for m in range(2, n + 1, 2)}
    w = {'w0': c[n//2], 'w0p': sum(c[kk]*T[n - 2*kk] for kk in range(n//2)) + c[n//2]*(6 - 12*n)}
    dims = {'b': 80, 'a': 81, 'c0': 81, 'c1': 3, 'c2': 27}
    for z, dz in dims.items():
        w['w_' + z] = sum(c[kk]*sp.Rational(ps[n - 2*kk][z], n2**((n - 2*kk)//2)) for kk in range(n//2)) + c[n//2]*dz
    U = sp.symbols('U2 U3')
    for j in range(6):
        if j in (1, 2):
            low = {2: t2[j], 4: t4[j], 6: t6[j]}
            def trm(m, j=j, low=low):
                if m in low: return sp.Rational(low[m], n2**(m//2))
                if m == n: return U[j - 1]
                return None
        else:
            def trm(m, j=j):
                return sp.Rational(sum(M(nn, m) for nn in st[j]), n2**(m//2))
        terms = []
        ok = True
        for kk in range(n//2):
            v = trm(n - 2*kk)
            if v is None:
                ok = False; break
            terms.append(c[kk]*v)
        if not ok:
            # degree strictly between 6 and n unknown for d2, d3: introduce symbols
            terms = []
            for kk in range(n//2):
                m = n - 2*kk
                v = trm(m)
                if v is None:
                    v = sp.Symbol('U%d_%d' % (j + 1, m))
                terms.append(c[kk]*v)
        w['w_d%d' % (j + 1)] = sum(terms) + c[n//2]*DIMS[j]
    eqs = [sp.expand(sum(cc*w[l] for cc, l in zip(row, labels))) for row in R]
    unknowns = sorted(set().union(*[e.free_symbols for e in eqs]), key=str)
    sol = sp.solve(eqs, unknowns, dict=True)
    if verbose:
        print('ideal %s, weight %d: %d relations, unknowns %s ->' % (key, n, len(R), [str(u) for u in unknowns]),
              'CONSISTENT' if sol else 'INCONSISTENT', '' if not sol else {str(a): str(b) for a, b in sol[0].items()})
    return bool(sol), eqs

if __name__ == '__main__':
    keys = [(1, 0), (3, 0), (9, 4), (35, 4)]
    res = {}
    for n in [int(x) for x in sys.argv[1:]] or [8, 10]:
        for key in keys:
            ok, _ = check(key, n)
            res['%s n=%d' % (key, n)] = ok
    json.dump({k: v for k, v in res.items()}, open('weightn_single_check.json', 'w'), indent=1)
