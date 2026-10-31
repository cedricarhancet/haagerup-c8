"""Joint feasibility of the two candidate weight-one algebras with all relations of weights 2, 4, 6 (Axis C2).

For s = g' + g'' (g' of level k1, g'' of level k2) and each R_{dj} = sum_a V(a) x N_a (a <= k1, N_a a g''-module
with highest weights <= k2), the following power sums are determined by the exact relations of weights 2, 4, 6
(single and mixed), once the module structures of the small sectors R_{d4}, R_{d5}, R_{d6} are fixed:
   Tr H'^p (p = 2,4,6), Tr H''^q (q = 2,4,6), Tr H'^2 H''^2, Tr H'^4 H''^2, Tr H'^2 H''^4.
We enumerate the joint structures of the small sectors, propagate, and decide the large sectors exactly.
As a test, the announced decompositions (Riesen) of the candidate A11+A139 must satisfy all of them.
"""
import json, itertools, sys
from functools import lru_cache
import sympy as sp
import sl2_pair_feasibility as S, weight6_single as W6, weight6_mixed as WM

DIMS = [7, 42, 119, 5, 13, 14]
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))

def single_sextic(key):
    """single-ideal data with sextic traces (unique compatible configuration, weight6_single pass 2)"""
    D = W6.single_data()[key]
    k = D['k']; n2 = 2*k
    R6 = WM.R6; labels = WM.labels
    ps = {p: W6.sector_power_sums(D['lab'], key[1], p) for p in (2, 4, 6)}
    T, T4, T6 = sp.Rational(4, k), sp.Rational(8, k*k), sp.Rational(128, n2**3)
    Ssym = sp.symbols('S1:7')
    w = {'w0': sp.Rational(-5, 576), 'w0p': T6 - 5*T4/4 + 5*T/16 + sp.Rational(5*66, 576)}
    for z, dz in (('b', 80), ('a', 81), ('c0', 81), ('c1', 3), ('c2', 27)):
        w['w_' + z] = sp.Rational(ps[6][z], n2**3) - sp.Rational(5, 4)*sp.Rational(ps[4][z], n2**2) + sp.Rational(5, 16)*sp.Rational(ps[2][z], n2) - sp.Rational(5*dz, 576)
    for j in range(6):
        w['w_d%d' % (j + 1)] = Ssym[j] - sp.Rational(5, 4)*D['Q'][j] + sp.Rational(5, 16)*D['tau'][j] - sp.Rational(5*DIMS[j], 576)
    eqs = [sp.expand(sum(c*w[l] for c, l in zip(row, labels))) for row in R6]
    sol = sp.solve(eqs, list(Ssym[:3]), dict=True)[0]
    t2 = [int(D['tau'][j]*n2) for j in range(6)]; t4 = [int(D['Q'][j]*n2**2) for j in range(6)]
    return D, k, sol, Ssym, t2, t4

def structures(dim, mom, nmax):
    """sl2-modules of dimension dim with prescribed (Tr H^2, Tr H^4), highest weights <= nmax"""
    out = []
    def rec(rem, r2, r4, mx, cur):
        if r2 == 0 and r4 == 0:
            out.append(tuple(cur + [0]*rem)); return
        for nn in range(min(mx, rem - 1), 0, -1):
            if M(nn, 2) <= r2 and M(nn, 4) <= r4:
                rec(rem - nn - 1, r2 - M(nn, 2), r4 - M(nn, 4), nn, cur + [nn])
    rec(dim, mom[0], mom[1], nmax, [])
    return out

def joint_small(dim, k1, k2, tg1, tg2, X22):
    """joint structures sum V(a) x V(b) with marginal (Tr H'^2, H'^4), (Tr H''^2, H''^4), mixed X22"""
    return [st for st, mx in S.joint_structures(dim, tg1, tg2, ) if mx == X22 and max(a for a, b in st) <= k1 and max(b for a, b in st) <= k2]

def mom_joint(st, p, q):
    return sum(M(a, p)*M(b, q) for a, b in st)

@lru_cache(maxsize=None)
def t46set(dim, t2, bmax):
    """set of (Tr H^4, Tr H^6) over sl2-modules with dimension dim, Tr H^2 = t2, highest weights <= bmax"""
    if t2 == 0:
        return frozenset([(0, 0)])
    if dim <= 0 or bmax == 0:
        return frozenset()
    out = set()
    for b in range(min(bmax, dim - 1), 0, -1):
        if M(b, 2) <= t2:
            for (x, y) in t46set(dim - b - 1, t2 - M(b, 2), b):
                out.add((x + M(b, 4), y + M(b, 6)))
    return frozenset(out)

def decide(dim, k1, k2, tgt):
    """tgt: dict with keys (p,q) for the joint power sums Tr H'^p H''^q; exact decision with N_a method"""
    A = range(0, min(k1, dim - 1) + 1)
    # n_a from (dim, p = 2, 4, 6)
    sols = []
    rngs = [range(0, tgt[(2, 0)]//M(a, 2) + 1) for a in A if a > 0]
    for nr in itertools.product(*rngs):
        na = [None] + list(nr)
        if any(sum(M(a, p)*na[a] for a in A if a > 0) != tgt[(p, 0)] for p in (2, 4, 6)):
            continue
        n0 = dim - sum((a + 1)*na[a] for a in A if a > 0)
        if n0 < 0:
            continue
        na[0] = n0
        # T2_a: sum (a+1) T2_a = tgt(0,2), sum M2(a) T2_a = tgt(2,2), sum M4(a) T2_a = tgt(4,2)
        r2 = [range(0, (tgt[(2, 2)]//M(a, 2)) + 1) for a in A if a > 0]
        for t2r in itertools.product(*r2):
            T2 = [None] + list(t2r)
            if sum(M(a, 2)*T2[a] for a in A if a > 0) != tgt[(2, 2)] or sum(M(a, 4)*T2[a] for a in A if a > 0) != tgt[(4, 2)]:
                continue
            T2[0] = tgt[(0, 2)] - sum((a + 1)*T2[a] for a in A if a > 0)
            if T2[0] < 0:
                continue
            sets = [t46set(na[a], T2[a], k2) if na[a] > 0 else (frozenset([(0, 0)]) if T2[a] == 0 else frozenset()) for a in A]
            if any(not s for s in sets):
                continue
            # combine: sum (a+1) T4_a = tgt(0,4), sum M2(a) T4_a = tgt(2,4), sum (a+1) T6_a = tgt(0,6)
            reach = {(0, 0, 0): ()}
            for a in A:
                new = {}
                for (s4, x24, s6), path in reach.items():
                    for (t4, t6) in sets[a]:
                        key = (s4 + (a + 1)*t4, x24 + M(a, 2)*t4, s6 + (a + 1)*t6)
                        if key[0] <= tgt[(0, 4)] and key[1] <= tgt[(2, 4)] and key[2] <= tgt[(0, 6)] and key not in new:
                            new[key] = path + ((t4, t6),)
                reach = new
            goal = (tgt[(0, 4)], tgt[(2, 4)], tgt[(0, 6)])
            if goal in reach:
                return {'n_a': na, 'T2_a': T2, 'T4T6_a': reach[goal]}
    return None

def run(name, p_key, q_key, lab, riesen=None):
    print('==', name)
    D1, k1, sol1, Ss1, t2a, t4a = single_sextic(p_key)
    D2, k2, sol2, Ss2, t2b, t4b = single_sextic(q_key)
    res, d1, d2, Xn = WM.analyse(name, p_key, q_key, lab)
    X22 = [int(Xn[j]*4*k1*k2) for j in range(6)]
    small = {}
    for j in (3, 4, 5):
        small[j] = joint_small(DIMS[j], k1, k2, (t2a[j], t4a[j]), (t2b[j], t4b[j]), X22[j])
        print('   R_d%d: %d joint structures' % (j + 1, len(small[j])))
    P = sp.symbols('P1:7')
    found = []
    for combo in itertools.product(*[small[j] for j in (3, 4, 5)]):
        st = dict(zip((3, 4, 5), combo))
        sub1 = {Ss1[j]: sp.Rational(sum(M(a, 6)*(b + 1) for a, b in st[j]), (2*k1)**3) for j in (3, 4, 5)}
        sub2 = {Ss2[j]: sp.Rational(sum(M(b, 6)*(a + 1) for a, b in st[j]), (2*k2)**3) for j in (3, 4, 5)}
        t6a = [int(sp.Rational(sol1[Ss1[j]].subs(sub1))*(2*k1)**3) if j < 3 else int(sub1[Ss1[j]]*(2*k1)**3) for j in range(6)]
        t6b = [int(sp.Rational(sol2[Ss2[j]].subs(sub2))*(2*k2)**3) if j < 3 else int(sub2[Ss2[j]]*(2*k2)**3) for j in range(6)]
        subP42 = {P[j]: sp.Rational(mom_joint(st[j], 4, 2), (2*k1)**2*(2*k2)) for j in (3, 4, 5)}
        subP24 = {P[j]: sp.Rational(mom_joint(st[j], 2, 4), (2*k1)*(2*k2)**2) for j in (3, 4, 5)}
        s42 = res['x4y2']['solution']; s24 = res['x2y4']['solution']
        x42 = [sp.Rational(s42[P[j]].subs(subP42))*(2*k1)**2*(2*k2) if j < 3 else subP42[P[j]]*(2*k1)**2*(2*k2) for j in range(6)]
        x24 = [sp.Rational(s24[P[j]].subs(subP24))*(2*k1)*(2*k2)**2 if j < 3 else subP24[P[j]]*(2*k1)*(2*k2)**2 for j in range(6)]
        ok = True; wits = {}
        for j in (0, 1, 2):
            vals = [t6a[j], t6b[j], x42[j], x24[j]]
            if any(sp.Rational(v).q != 1 or v < 0 for v in vals):
                ok = False; wits['d%d' % (j + 1)] = 'non-integral or negative: %s' % [str(v) for v in vals]; break
            tgt = {(2, 0): t2a[j], (4, 0): t4a[j], (6, 0): t6a[j], (0, 2): t2b[j], (0, 4): t4b[j], (0, 6): t6b[j],
                   (2, 2): X22[j], (4, 2): int(x42[j]), (2, 4): int(x24[j])}
            if riesen is not None:
                rb = riesen['d%d' % (j + 1)]
                obs = {pq: sum(M(a, pq[0])*M(b, pq[1]) for a, b in rb) for pq in tgt}
                print('   test against the announced R_d%d:' % (j + 1), 'AGREE' if obs == tgt else 'DIFFER %s vs %s' % (obs, tgt))
            w = decide(DIMS[j], k1, k2, tgt)
            if w is None:
                ok = False; wits['d%d' % (j + 1)] = 'no module structure with targets %s' % tgt; break
            wits['d%d' % (j + 1)] = w
        found.append((ok, st, wits))
        print('   small sectors %s ->' % {('d%d' % (j + 1)): v for j, v in st.items()}, 'COMPATIBLE' if ok else 'EXCLUDED', '' if ok else wits)
    return found

if __name__ == '__main__':
    riesen = {'d1': [(0, 6)], 'd2': [(0, 0), (0, 4), (0, 6), (0, 8), (1, 7), (1, 1)],
              'd3': [(0, 0), (0, 2), (0, 12), (1, 9), (1, 7), (1, 5), (1, 3), (1, 1), (0, 4), (0, 4), (0, 6), (0, 6), (0, 8), (0, 8)],
              'd4': [(0, 4)], 'd5': [(0, 4), (1, 3)], 'd6': [(0, 2), (0, 6), (1, 1)]}
    out = {}
    for name, (p, q, lab) in WM.CANDS.items():
        f = run(name, p, q, lab, riesen if p == (1, 0) else None)
        out[name] = [{'ok': ok, 'small': {str(k): str(v) for k, v in st.items()}, 'detail': {k: str(v) for k, v in w.items()}} for ok, st, w in f]
    json.dump(out, open('sextic_joint.json', 'w'), indent=1)
