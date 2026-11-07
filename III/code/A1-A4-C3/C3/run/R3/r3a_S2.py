"""R3a step 3: S_2 x 1 inside W = V_{E8+E8}, exact over Q(sqrt(-6)), for each candidate omega_{C x 1}.
For each sl2^2-type tau = (a, b): P_tau = A x A-primary vectors of W_2 of weight (a, b, 0, 0), hw for (g'_1, g''_1),
invariant for (g'_2, g''_2) (rational basis).  o(G) preserves P_tau; the S_2 x 1 primaries of type tau are the
eigenvectors of o(omega_{C x 1}) on P_tau with eigenvalue 2 - h_A(tau) (= kernel of o(omega_{1 x C})).
Then S_2 x 1 = primaries + lowering descendants + affine descendants: dimension and sl2^2 decomposition vs "U2",
and the commutant conditions (omega_{1 x C})_(n) v = 0 (n = 1, 2, 3) on the primaries.
Output: S2.pkl"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time
from fractions import Fraction as Fr
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import e8e8_voa as V
import common as CM
import qs6 as Q
t0 = time.time()
CI = int(sys.argv[1]) if len(sys.argv) > 1 else 0
DG = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
OC = pickle.load(open(os.path.join(CM.HERE, 'omegaC1.pkl'), 'rb'))
u, wC1 = OC['cands'][CI]
G = DG['G']
om1C = Q.vclean(Q.vadd(Q.vreal(DG['omCC']), wC1, (mpq(-1), mpq(0))))
def qkernel(M):
    """kernel of a k x k matrix over Q(s) (list of rows of Q(s) elements) as list of Q(s) vectors (Weil restriction)"""
    r, c = len(M), len(M[0])
    A = fmpq_mat(2*r, 2*c)
    for i in range(r):
        for j in range(c):
            a, b = M[i][j]
            # (a + b s)(x + y s) = (a x - 6 b y) + (a y + b x) s
            A[2*i, 2*j] = CM.to_fq(a); A[2*i, 2*j + 1] = CM.to_fq(-6*b)
            A[2*i + 1, 2*j] = CM.to_fq(b); A[2*i + 1, 2*j + 1] = CM.to_fq(a)
    Rr, rk = A.rref()
    n = 2*c
    piv = []
    for i in range(rk):
        piv.append(next(j for j in range(n) if Rr[i, j] != 0))
    free = [j for j in range(n) if j not in piv]
    vecs = []
    for fj in free:
        x = [fmpq(0)]*n; x[fj] = fmpq(1)
        for i, pj in enumerate(piv): x[pj] = -Rr[i, fj]
        vecs.append([(CM.to_mq(x[2*j]), CM.to_mq(x[2*j + 1])) for j in range(c)])
    # the Q-kernel is a Q(s)-space: keep a Q(s)-basis (dimension over Q(s) = len/2)
    basis = []
    for v in vecs:
        cand = basis + [v]
        if qrank(cand) == len(cand): basis = cand
    return basis
def qrank(vecs):
    if not vecs: return 0
    c = len(vecs[0])
    A = fmpq_mat(2*len(vecs), 2*c)
    for i, v in enumerate(vecs):
        for j, (a, b) in enumerate(v):
            # rows: v and s v
            A[2*i, 2*j] = CM.to_fq(a); A[2*i, 2*j + 1] = CM.to_fq(b)
            A[2*i + 1, 2*j] = CM.to_fq(-6*b); A[2*i + 1, 2*j + 1] = CM.to_fq(a)
    return A.rref()[1]//2
U2 = {(0, 0): 4, (0, 2): 1, (0, 4): 2, (1, 1): 1, (1, 3): 1, (1, 5): 1, (1, 7): 1, (2, 0): 1, (2, 2): 2, (2, 6): 1, (3, 3): 1, (4, 0): 1}
AFF = {(2, 0): 1, (4, 0): 1, (0, 0): 2, (0, 2): 1, (0, 4): 1, (2, 2): 1}       # affine descendants in weight 2
PRIM_EXPECTED = {t: U2.get(t, 0) - AFF.get(t, 0) for t in set(U2) | set(AFF)}
types = sorted(set(PRIM_EXPECTED) | {(0, 6), (2, 4), (3, 1), (0, 8), (1, 9), (3, 5), (4, 2)})
prims = {}
summary = {}
for tau in types:
    a, b = tau
    P = CM.cond_kernel(2, (a, b, 0, 0), CM.primary_ops()) if tau != (0, 0) else G
    if not P:
        prims[tau] = []; summary[str(tau)] = (0, 0, PRIM_EXPECTED.get(tau, 0)); continue
    # matrix of o(omega_{C x 1}) on P
    cols = []
    for p in P:
        img = Q.mode(wC1, 1, Q.vreal(p))
        cr = CM.span_coords(P, img[0]); ci = CM.span_coords(P, img[1])
        assert cr is not None and ci is not None, 'P_tau not o(G)-stable'
        cols.append([(CM.to_mq(x), CM.to_mq(y)) for x, y in zip(cr, ci)])
    k = len(P)
    lam = 2 - CM.hA(a, b)
    M = [[Q.sub(cols[j][i], (mpq(lam.numerator, lam.denominator), mpq(0)) if i == j else Q.Z) for j in range(k)] for i in range(k)]
    ker = qkernel(M)
    vecs = []
    for c in ker:
        v = Q.vzero()
        for p, x in zip(P, c):
            if not Q.iszero(x): Q.vadd(v, Q.vreal(p), x)
        vecs.append(Q.vclean(v))
    prims[tau] = vecs
    summary[str(tau)] = (len(P), len(vecs), PRIM_EXPECTED.get(tau, 0))
    print('type', tau, ': dim P =', len(P), '; S_2 x 1 primaries:', len(vecs), '(expected', PRIM_EXPECTED.get(tau, 0), ')',
          '(%.0fs)' % (time.time() - t0), flush=True)
# commutant checks: omega_{1 x C}_(n) p = 0 for n = 1, 2, 3 ; (x) A x 1-primary already imposed
bad = 0
for tau, vecs in prims.items():
    for p in vecs:
        for nn in (1, 2, 3):
            if not Q.iszero_v(Q.mode(om1C, nn, p)): bad += 1
print('omega_{1 x C}_(n) p = 0 (n = 1,2,3) for all primaries: failures', bad, flush=True)
pickle.dump({'CI': CI, 'prims': prims, 'summary': summary}, open(os.path.join(CM.HERE, 'S2_%d.pkl' % CI), 'wb'))
print('time %.0fs' % (time.time() - t0))
