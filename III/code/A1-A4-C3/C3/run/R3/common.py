"""Common helpers for R3a / R2 (exact over Q): weight spaces, kernels of mode conditions, zero modes on W_1."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, itertools
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import e8e8_voa as V
import weights_W as WW
HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.path.join(HERE, 'scratch-R3')   # unused (never read nor written by any script)
SL = WW.SL
NAMES = WW.NAMES
CUR = [SL[nm][i] for nm in NAMES for i in range(3)]          # 12 currents (E, H, F) x 4
E_OPS = [SL[nm][0] for nm in NAMES]
F_OPS = [SL[nm][2] for nm in NAMES]
H_OPS = [SL[nm][1] for nm in NAMES]
LEV = {'g1p': 3, 'g1pp': 13, 'g2p': 3, 'g2pp': 13}
_BYW = {}
def by_weight(wt):
    if wt not in _BYW:
        d = defaultdict(list)
        for k in V.basis(wt): d[WW.weight_of(k)].append(k)
        _BYW[wt] = d
    return _BYW[wt]
def wspace(wt, mu): return by_weight(wt)[tuple(mu)]
def to_fq(x): return fmpq(int(x.numerator), int(x.denominator))
def to_mq(x): return mpq(int(x.p), int(x.q))
def kernel(cols_images, ncols):
    """cols_images[j] = dict (target key -> mpq) image of column j; returns basis of kernel as lists of fmpq"""
    rows = sorted(set(k for img in cols_images for k in img), key=repr)
    ridx = {k: i for i, k in enumerate(rows)}
    if not rows:
        return [[fmpq(int(i == j)) for i in range(ncols)] for j in range(ncols)]
    M = fmpq_mat(len(rows), ncols)
    for j, img in enumerate(cols_images):
        for k, x in img.items(): M[ridx[k], j] = to_fq(x)
    Rr, rank = M.rref()
    piv = []
    for i in range(rank):
        j = next(j for j in range(ncols) if Rr[i, j] != 0); piv.append(j)
    free = [j for j in range(ncols) if j not in piv]
    out = []
    for fj in free:
        x = [fmpq(0)]*ncols; x[fj] = fmpq(1)
        for i, pj in enumerate(piv): x[pj] = -Rr[i, fj]
        out.append(x)
    return out
def combine(keys, coefs):
    v = {}
    for k, c in zip(keys, coefs):
        if c != 0: v[k] = to_mq(c)
    return v
def cond_kernel(wt, mu, ops, extra=None):
    """vectors of W_wt of sl2^4-weight mu killed by x_(n) for (x, n) in ops (and by extra linear maps)"""
    keys = wspace(wt, mu)
    imgs = []
    for k in keys:
        img = {}
        for t, (x, n) in enumerate(ops):
            r = V.mode_vec(x, n, {k: mpq(1)})
            for kk, c in r.items(): img[(t, kk)] = c
        if extra:
            for t, f in enumerate(extra):
                for kk, c in f(k).items(): img[('x', t, kk)] = c
        imgs.append(img)
    ker = kernel(imgs, len(keys))
    return [combine(keys, c) for c in ker]
def primary_ops(hw_for=(0, 1), inv_for=(2, 3)):
    """conditions: hw for factors hw_for (E(0) = 0), invariant for inv_for (E(0) = 0 with weight 0), all 12 currents
    in modes 1 and 2 (A x A primary)"""
    ops = [(E_OPS[i], 0) for i in range(4)]
    ops += [(x, n) for x in CUR for n in (1, 2)]
    return ops
def hA(a, b): return Fr(a*(a + 2), 20) + Fr(b*(b + 2), 60)
def span_coords(basis_vecs, target):
    """exact coordinates of target in span(basis_vecs) or None"""
    keys = sorted(set(k for v in basis_vecs for k in v) | set(target), key=repr)
    kidx = {k: i for i, k in enumerate(keys)}
    n = len(basis_vecs)
    M = fmpq_mat(len(keys), n + 1)
    for j, v in enumerate(basis_vecs):
        for k, x in v.items(): M[kidx[k], j] = to_fq(x)
    for k, x in target.items(): M[kidx[k], n] = to_fq(x)
    Rr, r = M.rref()
    sol = [fmpq(0)]*n
    for i in range(r):
        j = next(j for j in range(n + 1) if Rr[i, j] != 0)
        if j == n: return None
        sol[j] = Rr[i, n]
    return sol
# ---------------------------------------------------------------- W_1 isotypic data and sectors
SECTOR_R = {
    'S1': None,
    'c1': {(0, 2): 1}, 'c2': {(0, 0): 1, (0, 4): 1, (1, 1): 1, (1, 3): 1, (2, 2): 1},
    'd1': {(0, 2): 1, (1, 1): 1}, 'd5': {(0, 4): 1, (1, 3): 1},
    'd4': {(0, 0): 1, (1, 1): 1}, 'd6': {(0, 2): 1, (1, 3): 1, (2, 0): 1}}
HZ = {'0': Fr(0), 'c1': Fr(1, 3), 'c2': Fr(2, 3), 'd1': Fr(6, 13), 'd5': Fr(7, 13), 'd4': Fr(5, 13), 'd6': Fr(8, 13)}
PI = {'c1': 'c2', 'c2': 'c1', 'd1': 'd5', 'd5': 'd1', 'd4': 'd6', 'd6': 'd4'}
def w1_sectors():
    """component (a,b,c,d) -> list of sectors ('S1x1' / '1xS1' / z)"""
    comp = defaultdict(list)
    comp[(2, 0, 0, 0)].append('S1x1'); comp[(0, 2, 0, 0)].append('S1x1')
    comp[(0, 0, 2, 0)].append('1xS1'); comp[(0, 0, 0, 2)].append('1xS1')
    for z, pz in PI.items():
        for (a, b) in SECTOR_R[z]:
            for (c, d) in SECTOR_R[pz]:
                comp[(a, b, c, d)].append(z)
    return comp
def omega_C1_eig(comp, sector):
    a, b, c, d = comp
    if sector in ('S1x1', '1xS1'): return Fr(0)
    return HZ[sector] - hA(a, b)
def omega_1C_eig(comp, sector):
    a, b, c, d = comp
    if sector in ('S1x1', '1xS1'): return Fr(0)
    return HZ[PI[sector]] - hA(c, d)
_HWW1 = {}
def hw_W1(comp):
    """basis of highest weight vectors of W_1 of weight comp"""
    if comp not in _HWW1:
        _HWW1[comp] = cond_kernel(1, comp, [(E_OPS[i], 0) for i in range(4)])
    return _HWW1[comp]
def zero_mode_matrix(u, comp):
    """matrix of o(u) = u_(1) (u in W_2 commuting with sl2^4) on the hw space of W_1 of weight comp"""
    B = hw_W1(comp)
    M = []
    for w in B:
        img = V.clean(V.mode_vec(u, 1, w))
        c = span_coords(B, img)
        assert c is not None, 'o(u) does not preserve the hw space'
        M.append(c)
    # M[j] = coordinates of o(u) B_j ; return as column matrix (list of rows)
    n = len(B)
    return [[M[j][i] for j in range(n)] for i in range(n)]
