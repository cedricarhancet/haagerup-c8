"""Exact linear algebra over K = Q(sqrt13, sqrt-3) (basis 1, r, q, rq; r^2 = 13, q^2 = -3) by restriction of scalars.
Elements: 4-tuples of gmpy2.mpq (Kfield convention).  Reduction modulo the degree-one prime p3 of the
coherent embedding of the exact engine (r -> s13, q -> s3).  Principle: equalities are checked exactly; ranks and
pivots are found modulo p (a nonzero minor modulo p is nonzero over K)."""
import sys
sys.path.insert(0, '/home/claude/hj/work/Calculs/generation-forte')
from gmpy2 import mpq
from flint import fmpq_mat, fmpz_mat, nmod_mat, fmpq
import Kfield as K
import modp
F = modp.Fp(2**23); P = F.p
_T = K._T
def mulmat(a):
    """4x4 rational matrix of multiplication by a (columns = images of 1, r, q, rq)"""
    cols = []
    for j, bj in enumerate((K.ONE, K.R, K.Q, K.RQ)):
        cols.append(K.mul(a, bj))
    return [[cols[j][i] for j in range(4)] for i in range(4)]
def kred(a):
    s13, s3 = F.s13, F.s3
    return (F.q(a[0]) + F.q(a[1])*s13 + F.q(a[2])*s3 + F.q(a[3])*s13 % P*s3) % P
def den_ok(a):
    """True if the denominators of a are prime to P (a is then integral at the prime of the embedding)"""
    return all(int(x.denominator) % P for x in a if x)
def kinv(a):
    M = fmpq_mat(4, 4, [fmpq(int(x.numerator), int(x.denominator)) for row in mulmat(a) for x in row])
    b = M.solve(fmpq_mat(4, 1, [1, 0, 0, 0]))
    return tuple(mpq(int(b[i, 0].p), int(b[i, 0].q)) for i in range(4))
def to_fmpq(x): return fmpq(int(x.numerator), int(x.denominator))
def weil(Mrows):
    """K-matrix (list of rows of K elements, n x m) -> rational 4n x 4m matrix W with W . w(v) = w(M v)"""
    n = len(Mrows); m = len(Mrows[0]) if n else 0
    W = [[fmpq(0)] * (4 * m) for _ in range(4 * n)]
    for i in range(n):
        for j in range(m):
            a = Mrows[i][j]
            if K.iszero(a): continue
            mm = mulmat(a)
            for s in range(4):
                for t in range(4):
                    if mm[s][t]: W[4 * i + s][4 * j + t] = to_fmpq(mm[s][t])
    return fmpq_mat(4 * n, 4 * m, [x for row in W for x in row])
def solve_K(A, B):
    """A: K-matrix n x n (rows), B: K-matrix n x k ; returns X (n x k, K) with A X = B, exact."""
    n = len(A); k = len(B[0]) if B else 0
    WA = weil(A)
    cols = []
    for j in range(k):
        rhs = fmpq_mat(4 * n, 1, [to_fmpq(B[i][j][s]) for i in range(n) for s in range(4)])
        cols.append(rhs)
    if not cols: return [[] for _ in range(n)]
    RHS = fmpq_mat(4 * n, k, [to_fmpq(B[i][j][s]) for i in range(n) for s in range(4) for j in range(k)]) if False else None
    big = fmpq_mat(4 * n, k, [0] * (4 * n * k))
    for j in range(k):
        for i in range(4 * n): big[i, j] = cols[j][i, 0]
    X = WA.solve(big)
    out = [[None] * k for _ in range(n)]
    for i in range(n):
        for j in range(k):
            out[i][j] = tuple(mpq(int(X[4 * i + s, j].p), int(X[4 * i + s, j].q)) for s in range(4))
    return out
def kv_red(v, idx, n):
    z = [0] * n
    for kk, a in v.items(): z[idx[kk]] = kred(a)
    return z
def left_kernel(rows):
    """rows: list of K-vectors (dict key -> K).  Exact basis of {c : sum c_i rows_i = 0} (K-coefficients) such that
    its reduction modulo p is a basis of the kernel modulo p.  Returns (kernel list of dicts i -> K, rank r,
    status) ; status 'ok' means the kernel dimension over K equals the kernel dimension modulo p (verified)."""
    h = len(rows)
    keys = sorted({k for v in rows for k in v})
    idx = {k: i for i, k in enumerate(keys)}; N = len(keys)
    if h == 0: return [], 0, 'ok'
    if N == 0: return [{i: K.ONE} for i in range(h)], 0, 'ok'
    A = nmod_mat(h, N, [x for v in rows for x in kv_red(v, idx, N)], P)
    R, r = A.rref()
    J = []; j = 0
    for i in range(r):
        while int(R[i, j]) == 0: j += 1
        J.append(j); j += 1
    # pivot rows I among the h rows for the columns J
    AJt = nmod_mat(r, h, [int(A[i, c]) for c in J for i in range(h)], P)
    R2, r2 = AJt.rref(); assert r2 == r
    I = []; j = 0
    for i in range(r):
        while int(R2[i, j]) == 0: j += 1
        I.append(j); j += 1
    others = [i for i in range(h) if i not in set(I)]
    if not others: return [], r, 'ok'
    MIJ = [[rows[i].get(keys[c], K.Z) for c in J] for i in I]           # r x r
    # x_i M[I,J] = M[i,J]  <=>  M[I,J]^T x_i^T = M[i,J]^T
    At = [[MIJ[a][b] for a in range(r)] for b in range(r)]
    Bt = [[rows[i].get(keys[c], K.Z) for i in others] for c in J]        # r x len(others)
    X = solve_K(At, Bt)                                                  # r x len(others)
    ker = []
    for t, i in enumerate(others):
        c = {i: K.ONE}
        for a, ia in enumerate(I):
            x = X[a][t]
            if not K.iszero(x): c[ia] = K.scal(mpq(-1), x)
        ker.append(c)
    # exact verification on all coordinates
    for c in ker:
        acc = {}
        for i, a in c.items():
            for kk, b in rows[i].items():
                acc[kk] = K.add(acc.get(kk, K.Z), K.mul(a, b))
        if any(not K.iszero(x) for x in acc.values()): return ker, r, 'FAILED'
    return ker, r, 'ok'
