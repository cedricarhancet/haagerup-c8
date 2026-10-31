"""Axis A1, EXACT over K = Q(sqrt13, sqrt-3): the weight-four space C_4 of the commutant C = Com_S(A) (case II,
subalgebra 9, choice 0) and the bi-weights of the two Virasoro vectors omega_1 = 2 f_1 (c = 32/19), omega_2 = 2 f_2
(c = 182/95), f_1 + f_2 = e the orthogonal idempotents of span(e, x) (exact, Q(sqrt13)).

Input: the exact K-basis of S_4 (pivots of N'_4 = S_4, file piv_w4/blk_XXX.pkl of certify_w4b.py, exact dimension
6010 by phase_d_w4.py), restricted to the block of g' + g'' weight (0, 0) -- C_4 lies in that block.
Steps.
 1. M = matrix over K of the maps y_(k), y in U_1 (rational), k = 0..4, on the 280 pivots of the block (0,0).
    Reduction modulo p3 (coherent embedding): rank, pivot columns and pivot rows.  dim_K ker <= dim_Fp ker.
 2. Exact kernel: for each free column, the square K-system on the pivot rows (Weil restriction, fmpz_mat.solve);
    the resulting vectors w are then checked EXACTLY: y_(k) w = 0 for all y in U_1, k = 0..4 (all rows).
    K-independence: independence of the reductions modulo p.  Hence dim_K C_4 = number of vectors found.
 3. L_i(0) = (omega_i)_(1) on C_4: exact images, exact coordinates in the basis w (square system on 4 dim C_4 key
    coordinates, then exact check of u = sum c_f w_f on all keys).  Exact identities: L_1 + L_2 = 4,
    [L_1, L_2] = 0, prod_lambda (L_1 - lambda) = 0; eigenspace dimensions; exact joint eigenvectors.
 4. Exact primary test for the eigenvectors of non-integral bi-weight: (omega_i)_(2), (omega_i)_(3) kill them."""
import sys, os, json, pickle, time, math
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../generation-forte')
sys.path.insert(0, GF); os.chdir(GF)
import numpy as np, scipy.sparse as sps, sympy as sp
from flint import nmod_mat, fmpz_mat, fmpq_mat, fmpq
from gmpy2 import mpq
import lattice_voa as V, exact_engine as E, Kfield as K
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']; P = E.P; F = E.F
B4 = E.Blocks(T, 4); mu0 = (0, 0); keys = B4.blk[mu0]; pos = B4.pos[mu0]; nk = len(keys)
bidx = sorted(B4.blk).index(mu0)
piv = []
with open('piv_w4/blk_%03d.pkl' % bidx, 'rb') as fh:
    while True:
        try: piv.append(pickle.load(fh))
        except EOFError: break
assert all(B4.wt(v) == mu0 for v in piv)
npv = len(piv); log('exact pivots of S_4 in the block (0,0):', npv, ' keys:', nk)
U1 = S['U1']; MODES = range(0, 5)
# ---------- 1. rational operators y_(k) on the block keys (sparse, by input key) and their reductions
ops = {}                      # (a, k) -> {input key: {output key: mpq}}
outkeys = {}
for a, y in enumerate(U1):
    for k in MODES:
        d = {}
        for kk in keys:
            img = {k2: c for k2, c in V.mode_vec(y, k, {kk: mpq(1)}).items() if c}
            if img:
                d[kk] = img
                for k2 in img: outkeys.setdefault((a, k, k2), len(outkeys))
        ops[(a, k)] = d
V._MODE.clear(); V._EMODE.clear()
nout = len(outkeys); log('rational operators built; output coordinates:', nout)
R_, C_, D_ = [], [], []
for (a, k), d in ops.items():
    for kk, img in d.items():
        for k2, c in img.items():
            R_.append(outkeys[(a, k, k2)]); C_.append(pos[kk]); D_.append(F.q(c))
Mop = sps.csr_matrix((np.array(D_, dtype=np.int64), (np.array(R_), np.array(C_))), shape=(nout, nk))
rows4 = np.array([E.modp_vec(v, pos, nk) for v in piv], dtype=np.int64)          # npv x nk  (reductions)
Z = np.asarray(Mop @ rows4.T) % P                                                  # nout x npv
nzr = np.nonzero(np.any(Z, axis=1))[0]
Zm = nmod_mat(len(nzr), npv, [int(x) for x in Z[nzr].flatten()], P)
Rr, rank = Zm.rref()
pivc = []
for i in range(rank):
    pivc.append(next(j for j in range(npv) if int(Rr[i, j]) != 0))
free = [j for j in range(npv) if j not in pivc]
log('rank mod p:', rank, ' kernel dim mod p (upper bound for dim_K C_4):', len(free))
# pivot rows (greedy echelon on the restriction to the pivot columns)
ech = E.Ech(rank); prow = []
for i in nzr:
    if ech.add(Z[i, pivc].copy()): prow.append(int(i))
    if len(prow) == rank: break
assert len(prow) == rank
inv_out = {v: k for k, v in outkeys.items()}
# ---------- 2. exact entries of M on the pivot rows: M[o, j] = sum_kk O_o,kk v_j[kk]   (K-valued)
need = {}
for o in prow: need.setdefault(inv_out[o][:2], {})[inv_out[o][2]] = o
colsof = {}                   # o -> list of (input key, coefficient)
for ak, outs in need.items():
    for kk, img in ops[ak].items():
        for k2, c in img.items():
            if k2 in outs: colsof.setdefault(outs[k2], []).append((kk, c))
def Mentry(o, v):
    acc = K.Z
    for kk, c in colsof.get(o, []):
        a = v.get(kk)
        if a is not None: acc = K.add(acc, K.scal(c, a))
    return acc
Mex = [[Mentry(o, piv[j]) for j in range(npv)] for o in prow]
log('exact entries on the pivot rows computed')
BASK = [K.ONE, K.R, K.Q, K.RQ]
def mulmat(a):                # 4x4 rational matrix of multiplication by a in the basis 1, r, q, rq
    return [[K.mul(a, BASK[j])[i] for j in range(4)] for i in range(4)]
def weil_system(rowsK, colsidx):
    n = len(rowsK); m = len(colsidx)
    A = [[mpq(0)]*(4*m) for _ in range(4*n)]
    for i in range(n):
        for jj, j in enumerate(colsidx):
            a = rowsK[i][j]
            if K.iszero(a): continue
            mm = mulmat(a)
            for s in range(4):
                for t in range(4): A[4*i + s][4*jj + t] = mm[s][t]
    return A
def to_fmpz_rows(A):
    out = []
    for r in A:
        den = 1
        for x in r:
            if x: den = den*int(x.denominator)//math.gcd(den, int(x.denominator))
        out.append([int(x*den) for x in r] + [den])
    return out
A = weil_system(Mex, pivc)                                 # (4 rank) x (4 rank)
Bcols = []                                                 # right-hand sides: -M[:, f] * 1  (lambda_f = 1)
for f in free:
    col = []
    for i in range(rank):
        a = Mex[i][f]; col += [-a[0], -a[1], -a[2], -a[3]]
    Bcols.append(col)
# clear denominators row by row (same scaling on the RHS)
N4 = 4*rank
Az = []; Bz = []
for i in range(N4):
    row = A[i] + [Bcols[t][i] for t in range(len(free))]
    den = 1
    for x in row:
        if x: den = den*int(x.denominator)//math.gcd(den, int(x.denominator))
    Az.append([int(x*den) for x in A[i]]); Bz.append([int(Bcols[t][i]*den) for t in range(len(free))])
log('Weil system %d x %d, %d right-hand sides; solving exactly' % (N4, N4, len(free)))
Am = fmpz_mat(N4, N4, [x for r in Az for x in r]); Bm = fmpz_mat(N4, len(free), [x for r in Bz for x in r])
Xsol = Am.solve(Bm)                                        # fmpq_mat exact
log('solved')
def fq(x):
    return mpq(int(x.p), int(x.q))
W = []
for t, f in enumerate(free):
    w = dict()
    E.kv_add(w, piv[f])
    for jj, j in enumerate(pivc):
        lam = (fq(Xsol[4*jj, t]), fq(Xsol[4*jj + 1, t]), fq(Xsol[4*jj + 2, t]), fq(Xsol[4*jj + 3, t]))
        if not K.iszero(lam): E.kv_add(w, piv[j], lam)
    W.append(E.kv_clean(w))
log('candidate kernel vectors:', len(W), ' sizes', [len(w) for w in W])
# exact check on all rows: y_(k) w = 0
fails = 0
for t, w in enumerate(W):
    for a, y in enumerate(U1):
        for k in MODES:
            r = E.kprod([(K.ONE, y)], k, w)
            if r: fails += 1
    V._MODE.clear(); V._EMODE.clear()
log('EXACT check y_(k) w = 0 (y in U_1, k = 0..4): failures =', fails)
assert fails == 0
Wp = np.array([E.modp_vec(w, pos, nk) for w in W], dtype=np.int64)
rkW = nmod_mat(len(W), nk, [int(x) for x in Wp.flatten()], P).rank()
log('rank of the reductions of the kernel vectors:', rkW, ' => dim_K C_4 =', rkW, '(exact)')
assert rkW == len(W) == len(free)
dimC4 = rkW
# ---------- 3. idempotents (exact, Q(sqrt13)) and L_i(0) on C_4
gb = json.load(open('griess_B_II.json')); n = len(gb['gram'])
prod = {tuple(map(int, kk.split(','))): [sp.sympify(c) for c in v] for kk, v in gb['prod'].items()}
gram = sp.Matrix([[sp.sympify(c) for c in r] for r in gb['gram']])
G9 = json.load(open('griess_subalgebras.json'))['II']['subalgebras'][9]
xv = sp.Matrix([0] + [sp.sympify(c) for c in G9['x']]); ev = sp.Matrix([1] + [0]*(n - 1))
def nf(t): return sp.expand(sp.radsimp(sp.expand(t)))
def mult(a, b):
    out = sp.zeros(n, 1)
    for i in range(n):
        for j in range(n):
            if a[i] != 0 and b[j] != 0: out += a[i]*b[j]*sp.Matrix(prod[(i, j)])
    return out.applyfunc(nf)
r13 = sp.sqrt(13)
xx = mult(xv, xv); i_ = next(i for i in range(1, n) if xv[i] != 0)
beta = nf(xx[i_]/xv[i_]); alpha = nf(xx[0] - beta*xv[0])
assert (mult(ev, ev) - ev).applyfunc(nf) == sp.zeros(n, 1) and (mult(ev, xv) - xv).applyfunc(nf) == sp.zeros(n, 1)
assert (xx - alpha*ev - beta*xv).applyfunc(nf) == sp.zeros(n, 1)
Dq = nf(sp.Rational(1, 4)/(beta**2/4 + alpha)); uu, vv = sp.symbols('uu vv')
sol = [s_ for s_ in sp.solve([uu**2 + 13*vv**2 - Dq.coeff(r13, 0), 2*uu*vv - Dq.coeff(r13, 1)], [uu, vv], dict=True)
       if s_[uu].is_Rational and s_[vv].is_Rational]
b0 = sol[0][uu] + sol[0][vv]*r13
idem = []
for b in (b0, -b0):
    v = ((1 - b*beta)/2*ev + b*xv).applyfunc(nf); assert (mult(v, v) - v).applyfunc(nf) == sp.zeros(n, 1); idem.append(v)
cc = [nf(8*(v.T*gram*v)[0]) for v in idem]
order = sorted(range(2), key=lambda i: float(cc[i])); idem = [idem[i] for i in order]; cc = [cc[i] for i in order]
log('idempotents f_1, f_2:', [list(v) for v in idem], ' central charges', cc)
Bb = pickle.load(open('griess_B_II_basis.pkl', 'rb'))
# exact images b_m(1) w and omega_i(n) w
_BOP = {}
def bop(m, nmode):
    if (m, nmode) not in _BOP:
        d = {}
        for kk in keys:
            img = {k2: c for k2, c in V.mode_vec(Bb[m], nmode, {kk: mpq(1)}).items() if c}
            if img: d[kk] = img
        V._MODE.clear(); V._EMODE.clear(); _BOP[(m, nmode)] = d
    return _BOP[(m, nmode)]
def omega_mode(i, nmode, w):
    out = {}
    for m in range(n):
        c = 2*idem[i][m]
        if c == 0: continue
        cK = K.from_sympy(c); O = bop(m, nmode)
        for kk, a in w.items():
            img = O.get(kk)
            if not img: continue
            ca = K.mul(cK, a)
            for k2, x in img.items(): out[k2] = K.add(out.get(k2, K.Z), K.scal(x, ca))
    return E.kv_clean(out)
# coordinates in the basis W: 4*dim pivot key-coordinates (Weil), exact solve, exact check on all keys
Wcols = []; Wden = []
for w in W:
    for z in BASK:
        col, dd = E.weil_column(E.kv_scale(z, w), pos, nk); Wcols.append(col); Wden.append(dd)
mW = len(Wcols)
Mt = nmod_mat(mW, 4*nk, [c % P for col in Wcols for c in col], P); Rw, rkw = Mt.rref(); assert rkw == mW
prw = []; j = 0
for i in range(rkw):
    while int(Rw[i, j]) == 0: j += 1
    prw.append(j); j += 1
BP = fmpz_mat(mW, mW, [Wcols[c][r] for r in prw for c in range(mW)])
def coords(u):
    """K-coordinates of u in the basis W (exact, verified on all keys)"""
    uc, den = E.weil_column(u, pos, nk) if u else ([0]*(4*nk), 1)
    X = BP.solve(fmpz_mat(mW, 1, [uc[r] for r in prw]))
    lam = [fq(X[i, 0])*Wden[i]/den for i in range(mW)]
    cK = [(lam[4*t], lam[4*t + 1], lam[4*t + 2], lam[4*t + 3]) for t in range(len(W))]
    chk = {}
    for t, w in enumerate(W): E.kv_add(chk, w, cK[t])
    E.kv_add(chk, u, K.scal(mpq(-1), K.ONE))
    assert not E.kv_clean(chk), 'u not in span_K(W)'
    return cK
L = []
for i in range(2):
    Li = [coords(omega_mode(i, 1, w)) for w in W]       # Li[t] = coordinates of L_i(0) w_t
    L.append([[Li[t][s] for t in range(dimC4)] for s in range(dimC4)])   # matrix: column t
    log('L_%d(0) on C_4: exact matrix computed' % (i + 1))
def kmat_mul(X, Y):
    d = len(X); return [[__import__('functools').reduce(K.add, [K.mul(X[a][c], Y[c][b]) for c in range(d)], K.Z) for b in range(d)] for a in range(d)]
def kmat_sub_scalar(X, lam):
    d = len(X); return [[K.add(X[a][b], K.scal(-lam, K.ONE)) if a == b else X[a][b] for b in range(d)] for a in range(d)]
def kmat_zero(X): return all(K.iszero(x) for r in X for x in r)
def kmat_add(X, Y): return [[K.add(x, y) for x, y in zip(rx, ry)] for rx, ry in zip(X, Y)]
res = {'dim_C4_exact': dimC4, 'npivots_block00': npv, 'rank_modp': rank, 'central_charges': [str(c) for c in cc]}
Ssum = kmat_add(L[0], L[1])
res['L1+L2 == 4 (exact)'] = kmat_zero(kmat_sub_scalar(Ssum, mpq(4)))
res['[L1,L2] == 0 (exact)'] = kmat_zero(kmat_add(kmat_mul(L[0], L[1]), [[K.scal(mpq(-1), x) for x in r] for r in kmat_mul(L[1], L[0])]))
# eigenvalues: candidates from the modular computation, then exact annihilating polynomial
cands = [[mpq(0), mpq(2), mpq(4), mpq(25, 19), mpq(50, 19)], [mpq(0), mpq(2), mpq(4), mpq(26, 19), mpq(51, 19)]]
for i in range(2):
    Pm = None
    for lam in cands[i]:
        X = kmat_sub_scalar(L[i], lam); Pm = X if Pm is None else kmat_mul(Pm, X)
    res['prod_lambda (L%d - lambda) == 0 (exact), lambda in %s' % (i + 1, [str(x) for x in cands[i]])] = kmat_zero(Pm)
    mult_ = {}
    for lam in cands[i]:
        X = kmat_sub_scalar(L[i], lam)
        Xp = [[E.kred(x) for x in r] for r in X]
        mult_[str(lam)] = dimC4 - nmod_mat(dimC4, dimC4, [x for r in Xp for x in r], P).rank()
    res['eigenspace dims of L%d (mod p upper bounds, sum = dim when diagonalisable)' % (i + 1)] = mult_
log(res)
# joint eigenvectors of non-integral bi-weight (exact): kernel of L1 - lambda, as K-vector in V'_4
special = {}
for lam1, lam2 in ((mpq(25, 19), mpq(51, 19)), (mpq(50, 19), mpq(26, 19))):
    X = kmat_sub_scalar(L[0], lam1)
    Xp = nmod_mat(dimC4, dimC4, [E.kred(x) for r in X for x in r], P)
    nul, nd = Xp.nullspace()
    assert nd == 1
    # exact kernel vector: solve with the pivot structure of the modular nullspace (Weil restriction, 7x7 over K)
    Rx, rkx = Xp.rref(); pc = []
    for ii in range(rkx): pc.append(next(jj for jj in range(dimC4) if int(Rx[ii, jj]) != 0))
    fr = [jj for jj in range(dimC4) if jj not in pc][0]
    rr = []; ech2 = E.Ech(rkx)
    Xpn = np.array([[E.kred(x) for x in r] for r in X], dtype=np.int64)
    for ii in range(dimC4):
        if ech2.add(Xpn[ii, pc].copy()): rr.append(ii)
    Aw = weil_system([X[ii] for ii in rr], pc)
    rhs = []
    for ii in rr: a = X[ii][fr]; rhs += [-a[0], -a[1], -a[2], -a[3]]
    Aq = fmpq_mat(len(Aw), len(Aw), [fmpq(int(x.numerator), int(x.denominator)) for r in Aw for x in r])
    bq = fmpq_mat(len(rhs), 1, [fmpq(int(x.numerator), int(x.denominator)) for x in rhs])
    sq = Aq.solve(bq)
    cvec = [K.Z]*dimC4; cvec[fr] = K.ONE
    for jj, j in enumerate(pc): cvec[j] = tuple(fq(sq[4*jj + s, 0]) for s in range(4))
    # exact check X c = 0
    assert all(K.iszero(__import__('functools').reduce(K.add, [K.mul(X[a][b], cvec[b]) for b in range(dimC4)], K.Z)) for a in range(dimC4))
    v = {}
    for t, w in enumerate(W): E.kv_add(v, w, cvec[t])
    v = E.kv_clean(v)
    # exact eigen-identities in V'_4 and primary conditions
    chk = {}
    ok = {}
    for i, lam in ((0, lam1), (1, lam2)):
        u = omega_mode(i, 1, v); d = dict(u); E.kv_add(d, v, K.scal(-lam, K.ONE))
        ok['L%d(0) v = %s v' % (i + 1, lam)] = not E.kv_clean(d)
        for nm in (2, 3):
            ok['L%d(%d) v = 0' % (i + 1, nm - 1)] = not omega_mode(i, nm, v)
    special['(%s, %s)' % (lam1, lam2)] = {'checks': ok, 'n_terms': len(v)}
    pickle.dump(v, open(HERE + '/c2_primary_%d_%d.pkl' % (lam1.numerator, lam2.numerator), 'wb'))
    log('bi-weight (%s, %s):' % (lam1, lam2), ok)
res['primaries_of_nonintegral_biweight'] = special
json.dump(res, open(HERE + '/c2_exact_C4.json', 'w'), indent=1)
pickle.dump({'W': W, 'L1': L[0], 'L2': L[1]}, open(HERE + '/c2_exact_C4_basis.pkl', 'wb'))
log('done')
