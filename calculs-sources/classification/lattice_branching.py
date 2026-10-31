"""Joint (H', H'')-character of the weight-n spaces of V' = V_{E6+A2} (n <= 3) from the lattice, and its
decomposition into sl2+sl2-modules V(a) x V(b).  H' is represented in the standard Cartan (cartan_conjugate.py)."""
import itertools, json
from collections import Counter
import numpy as np
import e6_tools as E
def p_rank(r, N):
    a = [1] + [0]*N
    for n in range(1, N + 1):
        for _ in range(r):
            for m in range(n, N + 1):
                a[m] += a[m - n]
    return a
P8 = p_rank(8, 4)
def e6_vectors(maxnorm):
    C = np.array(E.C)
    rng = range(-6, 7)
    out = []
    pts = np.array(list(itertools.product(rng, repeat=6)), dtype=np.int64)
    nr = np.einsum('ij,jk,ik->i', pts, C, pts)
    sel = pts[nr <= maxnorm]
    return sel, nr[nr <= maxnorm]
def a2_vectors(maxnorm):
    C = np.array([[2, -1], [-1, 2]])
    pts = np.array(list(itertools.product(range(-6, 7), repeat=2)), dtype=np.int64)
    nr = np.einsum('ij,jk,ik->i', pts, C, pts)
    return pts[nr <= maxnorm], nr[nr <= maxnorm]
def decompose2(wts):
    cnt = Counter(wts); dec = Counter()
    while cnt:
        a = max(k[0] for k in cnt)
        b = max(k[1] for k in cnt if k[0] == a)
        m = cnt[(a, b)]
        dec[(a, b)] += m
        for i in range(a + 1):
            for j in range(b + 1):
                key = (a - 2*i, b - 2*j)
                cnt[key] -= m
                assert cnt[key] >= 0, (key, a, b)
                if cnt[key] == 0: del cnt[key]
    return dec
def V_char(lab, lt, nmax=3):
    e6, n6 = e6_vectors(2*nmax)
    a2, n2 = a2_vectors(2*nmax)
    h1 = e6 @ np.array(lt); h2 = e6 @ np.array(lab)
    s2 = 2*a2[:, 0] + 2*a2[:, 1]
    res = {}
    for n in range(0, nmax + 1):
        cnt = Counter()
        for m6 in range(0, 2*n + 1, 2):
            sel6 = n6 == m6
            c6 = Counter(zip(h1[sel6].tolist(), h2[sel6].tolist()))
            for m2 in range(0, 2*n + 1 - m6, 2):
                k = n - m6//2 - m2//2
                if k < 0: continue
                sel2 = n2 == m2
                c2 = Counter(s2[sel2].tolist())
                for (x, y), u in c6.items():
                    for z, v in c2.items():
                        cnt[(x, y + z)] += u*v*P8[k]
        res[n] = cnt
    return res
if __name__ == '__main__':
    cc = json.load(open('cartan_conjugate.json'))
    out = {}
    for name, lab in (('A11+A139', (2, 1, 1, 0, 1, 2)), ('A13+A113', (1, 0, 0, 1, 0, 1))):
        lt = cc[name][0]
        ch = V_char(lab, lt, 3)
        dims = [sum(ch[n].values()) for n in range(4)]
        print(name, 'dims V\'_n:', dims)
        out[name] = {}
        for n in (1, 2, 3):
            dec = decompose2(ch[n])
            out[name][n] = {str(k): v for k, v in sorted(dec.items())}
            print('  n=%d: %d types' % (n, len(dec)), dict(sorted(dec.items())) if n <= 2 else '')
    json.dump(out, open('lattice_branching.json', 'w'), indent=1)
