"""Exact Q(sqrt(-6))-span with membership test (pivot rows + full verification)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import common as CM
import qs6 as Q
def smul(v): return ({k: -6*c for k, c in v[1].items()}, dict(v[0]))
class Span:
    """Q(s)-span of a list of Q(s)-vectors with exact membership: pivot rows modulo nothing (exact rref on a
    selection of rows), then full verification."""
    def __init__(self, vecs):
        self.B = []
        cols = []
        for v in vecs:
            cols += [v, smul(v)]
        keys = sorted(set(k for c in cols for part in c for k in part), key=repr)
        self.keys = keys; kidx = {k: i for i, k in enumerate(keys)}
        M = fmpq_mat(2*len(keys), len(cols))
        for j, c in enumerate(cols):
            for part in (0, 1):
                for k, x in c[part].items(): M[2*kidx[k] + part, j] = CM.to_fq(x)
        # independent columns
        R, rk = M.rref()
        pivc = [next(j for j in range(len(cols)) if R[i, j] != 0) for i in range(rk)]
        self.rank = rk//2
        keepv = sorted(set(j//2 for j in pivc))
        self.B = [vecs[j] for j in keepv]
        cols = []
        for v in self.B: cols += [v, smul(v)]
        self.cols = cols
        MT = fmpq_mat(len(cols), 2*len(keys))
        for j, c in enumerate(cols):
            for part in (0, 1):
                for k, x in c[part].items(): MT[j, 2*kidx[k] + part] = CM.to_fq(x)
        R2_, rk2 = MT.rref()
        assert rk2 == len(cols)
        self.pivr = [next(j for j in range(2*len(keys)) if R2_[i, j] != 0) for i in range(rk2)]
        A = fmpq_mat(rk2, rk2)
        for i, r in enumerate(self.pivr):
            for j in range(rk2): A[i, j] = MT[j, r]
        self.Ainv = A.inv()
        self.kidx = kidx
    def coords(self, w):
        w = Q.vclean(w)
        for part in w:
            for k in part:
                if k not in self.kidx: return None
        b = fmpq_mat(len(self.pivr), 1)
        for i, r in enumerate(self.pivr):
            k = self.keys[r//2]; part = r % 2
            b[i, 0] = CM.to_fq(w[part].get(k, mpq(0)))
        x = self.Ainv*b
        # verify
        rec = Q.vzero()
        for j, c in enumerate(self.cols):
            if x[j, 0] != 0: Q.vadd(rec, c, (CM.to_mq(x[j, 0]), mpq(0)))
        if not Q.eq(rec, w): return None
        return [(CM.to_mq(x[2*j, 0]), CM.to_mq(x[2*j + 1, 0])) for j in range(len(self.B))]
