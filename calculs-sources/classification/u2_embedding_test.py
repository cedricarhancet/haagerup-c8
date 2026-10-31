"""Is there an s-module U_2 with the computed dimension and moments, (L_s)_2 + C omega_C <= U_2 <= V'_2 ?
V' = V_{E6+A2}; V'_2 = (27 x 27bar, 1) + (78, 8) + (1, 1 + 8 + 8) as e6+sl3-modules (L_1(e6)_2 = 1+78+650 = 27x27bar,
L_1(sl3)_2 = 1 + 8 + 8).  Joint weights of 27 and 78 under (H', H''_e6) are computed exactly in e8_joint.py.
(L_s)_2 = (L_{k1}(sl2))_2 x 1 + 1 x (L_{k2}(sl2))_2 + V(2) x V(2), with L_k(sl2)_2 = V(2)+V(4)+V(0) (k >= 2),
V(2)+V(0) (k = 1); the coset conformal vector gives one more trivial summand (c(coset) != 0).
"""
import json, itertools
from collections import Counter
import e8_joint as J, weight6_mixed as WM
import u2_moments as UM
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
def decompose2(wts):
    cnt = Counter(wts); dec = Counter()
    while cnt:
        a = max(k[0] for k in cnt if cnt[k] > 0)
        b = max(k[1] for k in cnt if cnt[k] > 0 and k[0] == a)
        m = cnt[(a, b)]
        assert m > 0
        dec[(a, b)] += m
        for i in range(a + 1):
            for j in range(b + 1):
                key = (a - 2*i, b - 2*j)
                cnt[key] -= m
                assert cnt[key] >= 0
                if cnt[key] == 0:
                    del cnt[key]
        cnt = Counter({k: v for k, v in cnt.items() if v > 0})
    return dec
sl3ad = [2, 2, -2, -2, 4, -4, 0, 0]
def V2_decomposition(lab, i1p):
    ad6, w27 = J.joint_weights(lab, i1p)
    wts = [(x[0] - y[0], x[1] - y[1]) for x in w27 for y in w27]
    wts += [(x[0], x[1] + s) for x in ad6 for s in sl3ad]
    wts += [(0, 0)] + [(0, s) for s in sl3ad]*2
    assert len(wts) == 1370
    return decompose2(wts)
def Ls2(k1, k2):
    L = Counter()
    for (k, side) in ((k1, 0), (k2, 1)):
        parts = [2, 0] + ([4] if k >= 2 else [])
        for p in parts:
            L[(p, 0) if side == 0 else (0, p)] += 1
    L[(2, 2)] += 1
    L[(0, 0)] += 1           # coset conformal vector
    return L
def search(upper, lower, dim, mom):
    items = sorted(upper, key=lambda t: -((t[0] + 1)*(t[1] + 1)))
    PQ = list(mom)
    sols = []
    base = {pq: sum(lower[it]*M(it[0], pq[0])*M(it[1], pq[1]) for it in lower) for pq in PQ}
    rem_dim = dim - sum(lower[it]*(it[0] + 1)*(it[1] + 1) for it in lower)
    target = {pq: mom[pq] - base[pq] for pq in PQ}
    if rem_dim < 0 or any(v < 0 for v in target.values()):
        return []
    def rec(i, rd, tg, cur):
        if rd == 0 and all(v == 0 for v in tg.values()):
            sols.append(dict(cur)); return
        if i == len(items) or rd <= 0:
            return
        it = items[i]
        cap = upper[it] - lower.get(it, 0)
        dd = (it[0] + 1)*(it[1] + 1)
        cs = {pq: M(it[0], pq[0])*M(it[1], pq[1]) for pq in PQ}
        for m in range(min(cap, rd//dd), -1, -1):
            ntg = {pq: tg[pq] - m*cs[pq] for pq in PQ}
            if any(v < 0 for v in ntg.values()):
                continue
            if m: cur[it] = cur.get(it, 0) + m
            rec(i + 1, rd - m*dd, ntg, cur)
            if m:
                cur[it] -= m
                if cur[it] == 0: del cur[it]
            if len(sols) > 50: return
    rec(0, rem_dim, target, {})
    return sols
if __name__ == '__main__':
    out = {}
    for name, (pk, qk, lab) in WM.CANDS.items():
        k1 = pk[0] + pk[1]; k2 = qk[0] + qk[1]
        V2 = V2_decomposition(lab, pk[0])
        res, _, _ = UM.solve_U2(name)
        mom = {pq: int(v) for pq, v in res.items()}
        L = Ls2(k1, k2)
        ok_lower = all(V2.get(k, 0) >= v for k, v in L.items())
        sols = search(V2, L, 120, mom) if ok_lower else []
        print('%s: V\'_2 = %s' % (name, dict(sorted(V2.items()))))
        print('   lower bound (L_s)_2 + omega_C fits in V\'_2:', ok_lower, '; U_2 decompositions found:', len(sols))
        for s in sols[:5]:
            full = Counter(L); full.update(s)
            print('     U_2 =', dict(sorted(full.items())))
        out[name] = {'V2': {str(k): v for k, v in V2.items()}, 'n_solutions': len(sols),
                     'solutions': [{str(k): v for k, v in (Counter(L) + Counter(s)).items()} for s in sols[:20]]}
    json.dump(out, open('u2_embedding_test.json', 'w'), indent=1)
