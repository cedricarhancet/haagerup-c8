"""Exact joint analysis of U_1 = g' + g'' (gamma = 0), with the level bound (Axis C2, final form).

Input: the surviving simple ideals with the level bound (sl2_ideal_feasibility_level.json); for each ideal the
power sums Tr H^2, Tr H^4 on every R_{dj} are determined (unique R_{d6} structure).  Mixed weight-four
relations determine X_j = Tr_{R_dj} H'^2 H''^2 (integer normalization) from X_6.
Each R_{dj} is written R = sum_{a <= k'} V(a) x N_a with N_a a g''-module with highest weights <= k''.
Unknowns: n_a = dim N_a, T2_a = Tr_{N_a} H''^2, T4_a = Tr_{N_a} H''^4.  Linear constraints:
   sum (a+1) n_a = dim R,  sum M2(a) n_a = S2',  sum M4(a) n_a = S4',
   sum (a+1) T2_a = S2'',  sum M2(a) T2_a = X,  sum (a+1) T4_a = S4'',
and each (n_a, T2_a, T4_a) must be realized by a g''-module.  Exhaustive exact search.
"""
import json, itertools, sys
from functools import lru_cache
import sympy as sp
import sl2_pair_feasibility as S

def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))

@lru_cache(maxsize=None)
def t4set(dim, t2, bmax):
    """set of Tr H^4 over sl2-modules of dimension dim, Tr H^2 = t2, highest weights <= bmax"""
    if t2 == 0:
        return frozenset([0])
    if dim <= 0 or bmax == 0:
        return frozenset()
    out = set()
    for b in range(min(bmax, dim - 1), 0, -1):
        if M(b, 2) <= t2:
            for x in t4set(dim - b - 1, t2 - M(b, 2), b):
                out.add(x + M(b, 4))
    return frozenset(out)

def sector(n, S2a, S4a, S2b, S4b, X, ka, kb):
    A = list(range(0, min(ka, n - 1) + 1))
    sols = []
    # n_a for a >= 1 bounded by S2a; n_0 determined
    ranges = [range(0, S2a//M(a, 2) + 1) for a in A[1:]]
    for na_rest in itertools.product(*ranges):
        na = (None,) + na_rest
        s2 = sum(M(a, 2)*na[a] for a in A[1:]); s4 = sum(M(a, 4)*na[a] for a in A[1:])
        if s2 != S2a or s4 != S4a:
            continue
        n0 = n - sum((a + 1)*na[a] for a in A[1:])
        if n0 < 0:
            continue
        na = (n0,) + na_rest
        # T2_a for a >= 1: sum M2(a) T2_a = X ; T2_0 = S2b - sum (a+1) T2_a
        r2 = [range(0, X//M(a, 2) + 1) for a in A[1:]]
        for t2rest in itertools.product(*r2):
            if sum(M(a, 2)*t2rest[a - 1] for a in A[1:]) != X:
                continue
            t20 = S2b - sum((a + 1)*t2rest[a - 1] for a in A[1:])
            if t20 < 0:
                continue
            T2 = (t20,) + t2rest
            sets = [t4set(na[a], T2[a], kb) if na[a] > 0 else (frozenset([0]) if T2[a] == 0 else frozenset()) for a in A]
            if any(len(s) == 0 for s in sets):
                continue
            # subset-sum: sum (a+1) T4_a = S4b
            reach = {0: ()}
            for a in A:
                new = {}
                for v, path in reach.items():
                    for t in sets[a]:
                        w = v + (a + 1)*t
                        if w <= S4b and w not in new:
                            new[w] = path + (t,)
                reach = new
            if S4b in reach:
                sols.append({'n_a': na, 'T2_a': T2, 'T4_a': reach[S4b]})
                return sols          # one witness suffices
    return sols

if __name__ == '__main__':
    surv = json.load(open('sl2_ideal_feasibility_level.json'))
    keys = sorted(set((s['i1'], s['i2']) for s in surv))
    pairs = [(p, q) for p, q in itertools.product(keys, keys) if p <= q and not (p[1] != 0 and q[1] != 0)]
    # orient so that the first factor g' has the smaller level (enumeration over a <= k')
    out = {}
    for p, q in pairs:
        d1, d2 = S.single[p], S.single[q]
        if d1['k'] > d2['k']:
            d1, d2 = d2, d1
        sol = S.mixed_solution(d1, d2)
        kk = 4*d1['k']*d2['k']
        res6 = []
        for struct6, mix6 in S.joint_structures(14, (d1['t2'][5], d1['t4'][5]), (d2['t2'][5], d2['t4'][5])):
            if max(a for a, b in struct6) > d1['k'] or max(b for a, b in struct6) > d2['k']:
                continue
            X6 = sp.Rational(mix6, kk)
            ok, detail = True, {'d6': struct6}
            for j in (3, 0, 4, 1, 2):
                mj = sp.Rational(sol[S.X[j]].subs(S.X[5], X6))*kk
                if mj.q != 1 or mj < 0:
                    ok = False; detail['fail'] = 'd%d: mixed trace %s' % (j + 1, mj); break
                sc = sector(S.DIMS[j], d1['t2'][j], d1['t4'][j], d2['t2'][j], d2['t4'][j], int(mj), d1['k'], d2['k'])
                if not sc:
                    ok = False; detail['fail'] = 'd%d: no (sl2+sl2)-module structure' % (j + 1); break
                detail['d%d' % (j + 1)] = sc[0]
            res6.append((ok, detail))
        good = [dd for ok, dd in res6 if ok]
        print('pair %s + %s (levels %d, %d): R_d6 structures %d, compatible %d' % (p, q, d1['k'], d2['k'], len(res6), len(good)))
        for ok, dd in res6:
            print('    ', 'OK ' if ok else 'NO ', dd if not ok else {kx: v for kx, v in dd.items()})
        out[str((p, q))] = {'compatible': len(good), 'details': [str(dd) for ok, dd in res6]}
    json.dump(out, open('sl2_pair_exact.json', 'w'), indent=1)
