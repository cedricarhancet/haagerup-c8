"""Exact linear algebra helpers (python-flint)."""
from flint import fmpq_mat, fmpq
from gmpy2 import mpq
import lattice_voa as V
def to_fq(x): return fmpq(int(x.numerator), int(x.denominator))
def to_mq(x): return mpq(int(x.p), int(x.q))
def matrix_from_cols(cols, keys=None):
    if keys is None: keys = sorted(set().union(*[set(c) for c in cols])) if cols else []
    kidx = {k: i for i, k in enumerate(keys)}
    M = fmpq_mat(max(len(keys), 1), len(cols))
    for j, c in enumerate(cols):
        for k, x in c.items():
            M[kidx[k], j] = to_fq(x)
    return M, keys
def nullspace(M):
    """basis of {x : M x = 0} as lists of fmpq"""
    Rr, rank = M.rref()
    n = M.ncols(); piv = []
    for i in range(rank):
        j = next(j for j in range(n) if Rr[i, j] != 0); piv.append(j)
    free = [j for j in range(n) if j not in piv]
    out = []
    for fj in free:
        x = [fmpq(0)]*n; x[fj] = fmpq(1)
        for i, pj in enumerate(piv):
            x[pj] = -Rr[i, fj]
        out.append(x)
    return out
def rank(M): return M.rref()[1]
def combine(vecs, coef):
    out = {}
    for v, c in zip(vecs, coef):
        if c != 0: V.vadd(out, v, to_mq(c) if not isinstance(c, type(mpq(1))) else c)
    return V.clean(out)
def coords(vecs, target):
    """exact coordinates of target in span(vecs), or None"""
    M, keys = matrix_from_cols(vecs + [target])
    kidx = {k: i for i, k in enumerate(keys)}
    for k in target:
        if k not in kidx: return None
    A = fmpq_mat(M.nrows(), len(vecs)); b = fmpq_mat(M.nrows(), 1)
    for i in range(M.nrows()):
        for j in range(len(vecs)): A[i, j] = M[i, j]
        b[i, 0] = M[i, len(vecs)]
    aug = fmpq_mat(M.nrows(), len(vecs) + 1)
    for i in range(M.nrows()):
        for j in range(len(vecs) + 1): aug[i, j] = M[i, j]
    Rr, r = aug.rref()
    sol = [fmpq(0)]*len(vecs)
    for i in range(r):
        j = next(j for j in range(len(vecs) + 1) if Rr[i, j] != 0)
        if j == len(vecs): return None
        sol[j] = Rr[i, len(vecs)]
    return sol
