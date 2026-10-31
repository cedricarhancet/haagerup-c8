"""Uniqueness of the (sl2+sl2)-module structures of R_{d1}, R_{d2}, R_{d3} given all determined joint power sums
(degrees <= 6, single and mixed) -- exhaustive enumeration with the multiplicity-module method; for R_{d3} the
enumeration is capped (the cap is reported)."""
import itertools, json, sys
from functools import lru_cache
import sextic_joint as SJ, weight6_mixed as WM, u2_moments as UM
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
@lru_cache(maxsize=None)
def mods(dim, t2, t4, t6, bmax):
    """all sl2-modules (tuples of highest weights >= 1) of dimension dim with the given power sums"""
    if t2 == 0:
        return ((),) if t4 == 0 and t6 == 0 and dim >= 0 else ()
    if dim <= 0 or bmax == 0:
        return ()
    out = []
    for b in range(min(bmax, dim - 1), 0, -1):
        if M(b, 2) <= t2 and M(b, 4) <= t4 and M(b, 6) <= t6:
            for rest in mods(dim - b - 1, t2 - M(b, 2), t4 - M(b, 4), t6 - M(b, 6), b):
                out.append((b,) + rest)
    return tuple(out)
def enumerate_all(dim, k1, k2, tgt, cap=2000):
    A = range(0, min(k1, dim - 1) + 1)
    sols = []
    rngs = [range(0, tgt[(2, 0)]//M(a, 2) + 1) for a in A if a > 0]
    for nr in itertools.product(*rngs):
        na = [None] + list(nr)
        if any(sum(M(a, p)*na[a] for a in A if a > 0) != tgt[(p, 0)] for p in (2, 4, 6)):
            continue
        na[0] = dim - sum((a + 1)*na[a] for a in A if a > 0)
        if na[0] < 0: continue
        r2 = [range(0, tgt[(2, 2)]//M(a, 2) + 1) for a in A if a > 0]
        for t2r in itertools.product(*r2):
            T2 = [None] + list(t2r)
            if sum(M(a, 2)*T2[a] for a in A if a > 0) != tgt[(2, 2)] or sum(M(a, 4)*T2[a] for a in A if a > 0) != tgt[(4, 2)]:
                continue
            T2[0] = tgt[(0, 2)] - sum((a + 1)*T2[a] for a in A if a > 0)
            if T2[0] < 0: continue
            # enumerate (T4_a, T6_a) choices via module lists
            lists = []
            for a in A:
                L = {}
                for t4 in range(T2[a], (tgt[(0, 4)]//(a + 1)) + 1) if na[a] else [0]:
                    pass
                lists.append(None)
            # direct: enumerate modules N_a with (n_a, T2_a) and any T4, T6, then filter sums
            cand = []
            for a in A:
                opts = []
                if na[a] == 0:
                    opts = [((), 0, 0)] if T2[a] == 0 else []
                else:
                    # enumerate modules with given dim and T2 (T4, T6 free) -- bounded search
                    def gen(dim_, t2_, bmax):
                        if t2_ == 0:
                            yield (); return
                        if dim_ <= 0 or bmax == 0: return
                        for b in range(min(bmax, dim_ - 1), 0, -1):
                            if M(b, 2) <= t2_:
                                for rest in gen(dim_ - b - 1, t2_ - M(b, 2), b):
                                    yield (b,) + rest
                    for m in gen(na[a], T2[a], k2):
                        opts.append((m, sum(M(b, 4) for b in m), sum(M(b, 6) for b in m)))
                cand.append(opts)
            for choice in itertools.product(*cand):
                if sum((a + 1)*choice[a][1] for a in A) != tgt[(0, 4)] or sum(M(a, 2)*choice[a][1] for a in A) != tgt[(2, 4)] \
                   or sum((a + 1)*choice[a][2] for a in A) != tgt[(0, 6)]:
                    continue
                st = []
                for a in A:
                    m = choice[a][0]
                    st += [(a, b) for b in m] + [(a, 0)]*(na[a] - sum(b + 1 for b in m))
                sols.append(tuple(sorted(st, reverse=True)))
                if len(sols) >= cap:
                    return sols, True
    return sols, False
if __name__ == '__main__':
    out = {}
    for name, (pk, qk, lab) in WM.CANDS.items():
        tab, k1, k2 = UM.candidate_tables(name)
        for j in (1, 2, 3):
            tg = {pq: int(tab['d%d' % j][pq]) for pq in tab['d%d' % j]}
            dim = [7, 42, 119][j - 1]
            if j == 3 and '--d3' not in sys.argv:
                continue
            sols, capped = enumerate_all(dim, k1, k2, tg, cap=50)
            print('%s R_d%d: %d structure(s)%s' % (name, j, len(sols), ' (cap reached)' if capped else ''))
            for s in sols[:3]:
                from collections import Counter
                print('     ', dict(Counter(s)))
            out['%s d%d' % (name, j)] = {'n': len(sols), 'capped': capped, 'examples': [str(s) for s in sols[:5]]}
    json.dump(out, open('joint_uniqueness.json', 'w'), indent=1)
