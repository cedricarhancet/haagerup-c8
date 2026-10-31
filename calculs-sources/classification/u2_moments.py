"""Moments of U_2 as an s-module from the exact relations of weights 2, 4, 6 (Axis C4/D).

For h = x h' + y h'' (h', h'' orthonormal Cartan elements of the two ideals), the coefficient of x^p y^q in
   W_n(h) = sum_k C(n,2k)(2k-1)!! (-E2/12)^k <h,h>^k Tr h(0)^{n-2k} q^{L0-1/3}
is a holomorphic weight-n form (n = p + q).  Its leading coefficients involve the joint power sums
Tr_{R_z} h'^i h''^j (known for all z by the previous steps); its coefficient at q^{5/3} in the vacuum component
involves Tr_{U_2} h'^p h''^q.  The relations with the extra position (0, 1) determine the latter.
Test: for the announced algebra the result must agree with the decomposition of U_2 announced by Riesen.
"""
import json, itertools
from math import comb
import sympy as sp
import obstruction as O, vvmf_exact as V, e8_joint as J, sl2_pair_feasibility as S, sextic_joint as SJ, weight6_mixed as WM

def dfact(m):
    r = 1
    while m > 1:
        r *= m; m -= 2
    return r
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
names = V.names; d = V.d
A0 = [(0, -1), (0, 0)] + [(z, 0) for z in range(1, d)]
A = A0 + [(0, 1)]
DIMS = {'b': 80, 'a': 81, 'c0': 81, 'c1': 3, 'c2': 27, 'd1': 7, 'd2': 42, 'd3': 119, 'd4': 5, 'd5': 13, 'd6': 14}
PQ = [(2, 0), (0, 2), (4, 0), (0, 4), (2, 2), (6, 0), (0, 6), (4, 2), (2, 4)]

def candidate_tables(name):
    p_key, q_key, lab = WM.CANDS[name]
    js = J.joint_sums(lab, p_key[0], degs=PQ)
    d1, d2 = S.single[p_key], S.single[q_key]
    k1, k2 = d1['k'], d2['k']
    tab = {}
    for z in ('b', 'a', 'c0', 'c1', 'c2'):
        tab[z] = {}
    for pq, v in js.items():
        tab['b'][pq] = v['e6'] + v['sl3'] - v['s']; tab['a'][pq] = v['27x3']; tab['c0'][pq] = v['27x3']
        tab['c1'][pq] = v['1x3']; tab['c2'][pq] = v['27']
    # s itself
    tab['s'] = {pq: (2*2**pq[0] if pq[1] == 0 else 0) + (2*2**pq[1] if pq[0] == 0 else 0) for pq in PQ}
    # d-sectors: rerun the sextic joint step and collect the targets
    found = SJ.run(name, p_key, q_key, lab)
    ok = [f for f in found if f[0]]
    assert len(ok) == 1
    _, st, wits = ok[0]
    for j in (3, 4, 5):
        tab['d%d' % (j + 1)] = {pq: sum(M(a, pq[0])*M(b, pq[1]) for a, b in st[j]) for pq in PQ}
    for j in (0, 1, 2):
        w = wits['d%d' % (j + 1)]
        # reconstruct the targets from the witness (n_a, T2_a, T4T6_a): they realize the targets exactly
        A_ = range(len(w['n_a']))
        t = {}
        t[(2, 0)] = sum(M(a, 2)*w['n_a'][a] for a in A_); t[(4, 0)] = sum(M(a, 4)*w['n_a'][a] for a in A_)
        t[(6, 0)] = sum(M(a, 6)*w['n_a'][a] for a in A_)
        t[(0, 2)] = sum((a + 1)*w['T2_a'][a] for a in A_); t[(2, 2)] = sum(M(a, 2)*w['T2_a'][a] for a in A_)
        t[(4, 2)] = sum(M(a, 4)*w['T2_a'][a] for a in A_)
        t[(0, 4)] = sum((a + 1)*w['T4T6_a'][a][0] for a in A_); t[(2, 4)] = sum(M(a, 2)*w['T4T6_a'][a][0] for a in A_)
        t[(0, 6)] = sum((a + 1)*w['T4T6_a'][a][1] for a in A_)
        tab['d%d' % (j + 1)] = t
    return tab, k1, k2

def E2pow(k, N=3):
    e2 = [1, -24, -72]
    s = [1] + [0]*(N - 1)
    for _ in range(k):
        s = [sum(s[i]*e2[j - i] for i in range(j + 1)) for j in range(N)]
    return s

def solve_U2(name):
    tab, k1, k2 = candidate_tables(name)
    nrm = lambda i, j: sp.Integer(2*k1)**sp.Rational(i, 2)*sp.Integer(2*k2)**sp.Rational(j, 2)
    def tr(z, i, j):
        if (i, j) == (0, 0):
            return DIMS[z]
        if i % 2 or j % 2:
            return 0
        return sp.Rational(tab[z][(i, j)])/nrm(i, j)
    U2 = {pq: sp.Symbol('U_%d_%d' % pq) for pq in PQ}
    results = {}
    known = {}
    for (p, q) in PQ:
        n = p + q
        R = O.relations(n, [0], A, 2)
        # coefficient of x^p y^q in sum_k c_k rho^k Tr h^{n-2k}
        def coeff(trf, Eser_index):
            tot = 0
            for kk in range(n//2 + 1):
                ck = comb(n, 2*kk)*dfact(2*kk - 1)*sp.Rational(-1, 12)**kk
                m = n - 2*kk
                for a in range(kk + 1):
                    b = kk - a
                    i, j = p - 2*a, q - 2*b
                    if i < 0 or j < 0 or i + j != m:
                        continue
                    tot += ck*comb(kk, a)*comb(m, i)*trf(i, j, kk)
            return tot
        w = {}
        for (z, mm) in A:
            if z == 0 and mm == -1:
                w[(z, mm)] = coeff(lambda i, j, kk: 1 if (i, j) == (0, 0) else 0, 0)
            elif z == 0 and mm == 0:
                def f0(i, j, kk):
                    e = E2pow(kk)
                    if (i, j) == (0, 0):
                        return 6 + e[1]
                    return tr('s', i, j) if 's' in tab else 0
                tab['s'] = tab['s']
                w[(z, mm)] = coeff(lambda i, j, kk: (6 + E2pow(kk)[1]) if (i, j) == (0, 0) else (0 if (i % 2 or j % 2) else sp.Rational(tab['s'][(i, j)])/nrm(i, j)), 0)
            elif z == 0 and mm == 1:
                def f1(i, j, kk):
                    e = E2pow(kk)
                    if (i, j) == (0, 0):
                        return 120 + 6*e[1] + e[2]
                    if i % 2 or j % 2:
                        return 0
                    sv = sp.Rational(tab['s'][(i, j)])/nrm(i, j)
                    return U2[(i, j)] + e[1]*sv
                w[(z, mm)] = coeff(f1, 0)
            else:
                zn = names[z]
                w[(z, mm)] = coeff(lambda i, j, kk, zn=zn: tr(zn, i, j), 0)
        eqs = [sp.expand(sum(c*w[pos] for c, pos in zip(row, A)).subs(known)) for row in R]
        unknown = [u for u in U2.values() if any(e.has(u) for e in eqs)]
        sol = sp.solve(eqs, unknown, dict=True)
        if not sol:
            print(name, (p, q), 'INCONSISTENT'); results[(p, q)] = None; continue
        val = sol[0].get(U2[(p, q)])
        if val is not None:
            known[U2[(p, q)]] = val
        results[(p, q)] = sp.nsimplify(val*nrm(p, q)) if val is not None else None
        print('%s: Tr_{U_2} H\'^%d H\'\'^%d = %s' % (name, p, q, results[(p, q)]))
    return results, k1, k2

if __name__ == '__main__':
    riesenU2 = [(2, 0)] + [(0, 0)]*4 + [(0, 2), (2, 2), (0, 4), (0, 4), (1, 13), (1, 7), (1, 5), (1, 1), (0, 12), (0, 8), (0, 8)]
    out = {}
    for name in WM.CANDS:
        res, k1, k2 = solve_U2(name)
        out[name] = {str(k): str(v) for k, v in res.items()}
        if name.startswith('A11'):
            obs = {pq: sum(M(a, pq[0])*M(b, pq[1]) for a, b in riesenU2) for pq in PQ}
            print('   test against the announced decomposition of U_2:', 'AGREE' if all(res[pq] == obs[pq] for pq in PQ) else 'DIFFER %s' % obs)
    json.dump(out, open('u2_moments.json', 'w'), indent=1)
