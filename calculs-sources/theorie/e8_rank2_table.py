"""For each reduced positive-definite even binary lattice K = [[2a,b],[b,2c]] (0 <= b <= a <= c)
with det(K) <= DMAX, compute the maximal number of E8 roots orthogonal to a primitive embedding
K -> E8. Gluing K + K' = E8 shows that this is the maximal number of roots of an even
positive-definite rank-6 lattice K' with discriminant form -q_K.
The gluing obstruction applies as soon as this number is >= 14 (then 6 + #roots >= 20 > 19).
Usage: python3 e8_rank2_table.py DMAX
"""
import math, sys, time
import numpy as np
from fractions import Fraction

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
t0 = time.time()
# Bourbaki simple roots of E8 in (doubled) orthonormal coordinates
simple = np.array([
 [1,-1,-1,-1,-1,-1,-1,1],
 [2,2,0,0,0,0,0,0],
 [-2,2,0,0,0,0,0,0],
 [0,-2,2,0,0,0,0,0],
 [0,0,-2,2,0,0,0,0],
 [0,0,0,-2,2,0,0,0],
 [0,0,0,0,-2,2,0,0],
 [0,0,0,0,0,-2,2,0]], dtype=np.int64)
G = (simple @ simple.T)//4          # Cartan matrix (Gram of simple roots)
assert all(G[i, i] == 2 for i in range(8))
Ginv = np.rint(np.linalg.inv(G)).astype(np.int64)   # E8 unimodular: integral
assert (G @ Ginv == np.eye(8, dtype=np.int64)).all()
fund = Ginv @ simple                  # fundamental weights (rows), doubled coordinates

maxnorm = 2*((2*DMAX)//3 // 2 + 1)
# all vectors up to maxnorm (doubled coords: x even-type or odd-type, sum = 0 mod 4)
def gen(maxnorm):
    lim = 4*maxnorm
    out = []
    for parity in (0, 1):
        stack = [((), 0)]
        while stack:
            cur, s2 = stack.pop()
            i = len(cur)
            if i == 8:
                if s2 > 0 and sum(cur) % 4 == 0:
                    out.append(cur)
                continue
            r = math.isqrt(lim - s2)
            for x in range(-r, r+1):
                if (x - parity) % 2 == 0:
                    stack.append((cur + (x,), s2 + x*x))
    return np.array(out, dtype=np.int64)
V = gen(maxnorm)
norms = (V*V).sum(1)//4
roots = V[norms == 2]
assert len(roots) == 240
print('vectors up to norm', maxnorm, ':', len(V), ' (%.1fs)' % (time.time()-t0))
# dominant representatives of W-orbits: m in Z>=0^8, norm = m^T Ginv m
reps = {}
def rec(i, m):
    if i == 8:
        if any(m):
            mm = np.array(m)
            n = int(mm @ Ginv @ mm)
            if n <= maxnorm:
                reps.setdefault(n, []).append(mm @ fund)
        return
    for k in range(0, 4):
        mm = np.array(m + [k] + [0]*(7-i))
        if int(mm @ Ginv @ mm) > maxnorm: break
        rec(i+1, m + [k])
rec(0, [])
print('number of W-orbits by norm:', {n: len(r) for n, r in sorted(reps.items())})
coords_inv = np.linalg.inv(simple.T.astype(float))
def coords(v):
    c = coords_inv @ v.astype(float); ci = np.rint(c).astype(np.int64)
    assert np.allclose(c, ci); return ci
def primitive(u, v):
    A = np.stack([coords(u), coords(v)], 1); g = 0
    for i in range(8):
        for j in range(i+1, 8):
            g = math.gcd(g, int(A[i,0]*A[j,1] - A[j,0]*A[i,1]))
    return g == 1
by_norm = {n: V[norms == n] for n in range(2, maxnorm+1, 2)}
results = []; certs = []
for a in range(1, DMAX):
    for c in range(a, DMAX):
        if 2*c > maxnorm: continue
        for b in range(0, a+1):
            det = 4*a*c - b*b
            if det > DMAX or det <= 0: continue
            best = -1; pair = None
            Vc = by_norm[2*c]
            for u in reps.get(2*a, []):
                ipu = (Vc @ u)//4
                cand = Vc[ipu == b]
                if len(cand) == 0: continue
                ru = roots[(roots @ u)//4 == 0]
                counts = (((ru @ cand.T)//4) == 0).sum(0) if len(ru) else np.zeros(len(cand), int)
                for k in np.argsort(-counts):
                    if counts[k] <= best: break
                    if primitive(u, cand[k]):
                        best = int(counts[k]); pair = (u.copy(), cand[k].copy()); break
            results.append((det, 2*a, b, 2*c, best))
            certs.append((det, 2*a, b, 2*c, best, pair))
results.sort()
lines = ['det  K (Gram)             max #roots(K^perp)  gluing obstruction applies']
for det, A, b, C, best in results:
    lines.append(f'{det:4d}  [[{A},{b}],[{b},{C}]]'.ljust(30) + str(best).rjust(6) + '      ' + ('yes' if best >= 14 else 'no'))
print('\n'.join(lines))
print('time %.1fs' % (time.time()-t0))
open('e8_rank2_table_DMAX%d.txt' % DMAX, 'w').write('\n'.join(lines) + '\n')
# certificates: an explicit primitive pair (u, v) in doubled coordinates realizing the maximum,
# checked independently by check_e8_certificates.py
certs.sort(key=lambda t: t[:4])
with open('e8_rank2_certificates_DMAX%d.txt' % DMAX, 'w') as fh:
    fh.write('# det Gram(2a,b,2c) roots | u (doubled coords) | v (doubled coords)\n')
    for det, A, b, C, best, pair in certs:
        u, v = pair
        fh.write(f'{det} {A} {b} {C} {best} | ' + ' '.join(map(str, u)) + ' | ' + ' '.join(map(str, v)) + '\n')
