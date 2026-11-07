"""Common data for the Zhu computation (case II, subalgebra 9, choice 0): exact highest weight vectors of U_2 and
the lowering basis z_{a;i,j} = f1(0)^i f2(0)^j z_a (exact over K), mod-p reductions (p3 = 8388571, coherent
embedding), Cartan weights."""
import sys, os
GF = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'generation-forte'))
sys.path.insert(0, GF)
ZD = os.path.abspath(os.path.dirname(__file__))
os.chdir(GF)
import numpy as np
from gmpy2 import mpq
import lattice_voa as V, exact_engine as E, Kfield as K
if os.environ.get('HAAG_PSTART'):
    import modp
    E.F = modp.Fp(int(os.environ['HAAG_PSTART'])); E.P = E.F.p
    ZD = os.path.join(ZD, 'p%d' % E.P); os.makedirs(ZD, exist_ok=True)
P = E.P
def kv(v): return E.kv_from_terms([(K.ONE, v)])
def prod(u, n, v):               # exact product of two K-vectors
    return E.kprod(E.kv_terms(u), n, v)
def lin(*pairs):                 # sum of (rational scalar, K-vector)
    out = {}
    for c, v in pairs: E.kv_add(out, v, K.scal(mpq(c), K.ONE))
    return E.kv_clean(out)
def highest_vectors(S):
    T = S['T']; vac = {V.VAC: K.ONE}
    cur = {n: kv(T[n]) for n in ('H1', 'e1', 'f1', 'H2', 'e2', 'f2')}
    hv = []
    # affine part A_2
    hv.append(('A:e1(-2)1', (2, 0), prod(cur['e1'], -2, vac)))
    hv.append(('A:e1(-1)e1', (4, 0), prod(cur['e1'], -1, cur['e1'])))
    for i, lab in (('1', 'A:Cas1'), ('2', 'A:Cas2')):
        e, f, h = cur['e' + i], cur['f' + i], cur['H' + i]
        hv.append((lab, (0, 0), lin((1, prod(e, -1, f)), (1, prod(f, -1, e)), (mpq(1, 2), prod(h, -1, h)))))
    hv.append(('A:e2(-2)1', (0, 2), prod(cur['e2'], -2, vac)))
    hv.append(('A:e2(-1)e2', (0, 4), prod(cur['e2'], -1, cur['e2'])))
    hv.append(('A:e1(-1)e2', (2, 2), prod(cur['e1'], -1, cur['e2'])))
    # commutant and primary part
    for name, terms in S['hw'].items():
        v = E.kv_from_terms(terms)
        w = E.Blocks(T, 2).wt(v)
        hv.append(('P:' + name, w, v))
    return hv, cur
def lowering_basis(hv, cur):
    basis = []                   # (a, i, j, K-vector)
    for a, (lab, (r, l), z) in enumerate(hv):
        zi = z
        for i in range(r + 1):
            if i: zi = prod(cur['f1'], 0, zi)
            zij = zi
            for j in range(l + 1):
                if j: zij = prod(cur['f2'], 0, zij)
                basis.append((a, i, j, zij))
    return basis
def modp_sparse(kvec, pos):
    idx = np.fromiter((pos[k] for k in kvec), dtype=np.int64, count=len(kvec))
    val = np.fromiter((E.kred(a) for a in kvec.values()), dtype=np.int64, count=len(kvec))
    o = np.argsort(idx); return idx[o], val[o]
