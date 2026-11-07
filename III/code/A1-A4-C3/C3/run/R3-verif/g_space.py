# Step 2: G = Com_W(A x A) cap W_2 ; Sugawara vectors ; omega_CC
import sys, pickle, time
sys.path.insert(0, '.')
from voa import *
from lat import *
import flint
t0 = time.time()
D = pickle.load(open('sl2_4.pkl', 'rb'))
E, F, H = D['E'], D['F'], D['H']
Hs = [Hvec(X) for X in range(4)]
N0 = [a for a in NORM4 if all(ip(a, h) == 0 for h in Hs)]
R0 = [a for a in ROOTS if all(ip(a, h) == 0 for h in Hs)]
assert not R0
basis = []
for i in range(16):
    for j in range(i, 16):
        basis.append((((1, i), (1, j)), ZERO))
for i in range(16):
    basis.append((((2, i),), ZERO))
for b in N0:
    basis.append(((), b))
print('dim of H-weight-0 subspace of W_2:', len(basis))
gens = E + F + H
cols = []
for bk in basis:
    st = {bk: Fr(1)}
    col = {}
    for gi, x in enumerate(gens):
        for n in (0, 1, 2):
            r = mode_vec(x, n, st)
            for k, v in r.items():
                col[(gi, n, k)] = v
    cols.append(col)
rowkeys = sorted(set(k for c in cols for k in c), key=str)
ridx = {k: i for i, k in enumerate(rowkeys)}
print('equations (rows):', len(rowkeys))
import math
den = 1
for c in cols:
    for v in c.values():
        den = den*v.denominator//math.gcd(den, v.denominator)
M = flint.fmpz_mat(len(rowkeys), len(basis))
for j, c in enumerate(cols):
    for k, v in c.items():
        M[ridx[k], j] = int(v*den)
NS, nul = M.nullspace()
print('dim G =', nul)
# nicer basis: rref of the kernel (rows)
K = flint.fmpq_mat(nul, len(basis))
for a in range(nul):
    for j in range(len(basis)):
        K[a, j] = NS[j, a]
R, rank = K.rref()
pivots = []
Gb = []
for a in range(rank):
    row = {basis[j]: Fr(int(R[a, j].p), int(R[a, j].q)) for j in range(len(basis)) if R[a, j] != 0}
    piv = next(j for j in range(len(basis)) if R[a, j] != 0)
    pivots.append(basis[piv])
    Gb.append(row)
print('pivot coordinates:', pivots)
print('support sizes of G basis vectors:', [len(g) for g in Gb])
# Sugawara vectors
def sug(X):
    k = LEVELS[X]
    s = {}
    sadd(s, mode_vec(E[X], -1, F[X]))
    sadd(s, mode_vec(F[X], -1, E[X]))
    sadd(s, mode_vec(H[X], -1, H[X]), Fr(1, 2))
    return {kk: v/(2*(k+2)) for kk, v in s.items()}
oms = [sug(X) for X in range(4)]
for X in range(4):
    cc = mode_vec(oms[X], 3, oms[X])
    print('Sugawara #%d: omega_(3)omega = %s (c/2 expected %s);  L0 omega = 2 omega: %s' % (
        X, cc, Fr(3*LEVELS[X], 2*(LEVELS[X]+2)), mode_vec(oms[X], 1, oms[X]) == {k: 2*v for k, v in oms[X].items()}))
omAA = {}
for o in oms: sadd(omAA, o)
omW = {(((1, i), (1, i)), ZERO): Fr(1, 2) for i in range(16)}
omCC = sadd(dict(omW), omAA, -1)
print('omega_W_(3)omega_W =', mode_vec(omW, 3, omW))
cAA = 2*mode_vec(omAA, 3, omAA)[((), ZERO)]
cCC = 2*mode_vec(omCC, 3, omCC)[((), ZERO)]
print('c(AxA) =', cAA, '  c(CC) =', cCC)
# express omega_CC in G basis (pivot coords) and check
def coords_in_G(v):
    co = [v.get(p, Fr(0)) for p in pivots]
    rec = {}
    for c, g in zip(co, Gb):
        if c: sadd(rec, g, c)
    return co, clean(rec) == clean(v)
co, ok = coords_in_G(omCC)
print('omega_CC in G:', ok, ' coords:', [str(x) for x in co])
print('omega_CC is Virasoro with L0 = 2: ', mode_vec(omCC, 1, omCC) == {k: 2*v for k, v in omCC.items()})
# direct annihilation check of omega_CC by x_(n)
okk = all(not mode_vec(x, n, omCC) for x in gens for n in (0, 1, 2))
print('omega_CC annihilated by x_(n), n=0,1,2, x in sl2^4:', okk)
pickle.dump({'basis': basis, 'Gb': Gb, 'pivots': pivots, 'omAA': omAA, 'omW': omW, 'omCC': omCC, 'oms': oms}, open('G.pkl', 'wb'))
print('time %.1fs' % (time.time()-t0))
