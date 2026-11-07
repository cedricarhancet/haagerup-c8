"""R2: signature test of the invariant Hermitian form on the lowest weight spaces R_z of the S-modules M_z,
realized inside W_1 = e8 + e8 as slices R_z x w0 (w0 = fixed highest weight vector of one component of R_{pi z}),
exact over Q(sqrt(-6)).
  Slice: Y = {x in W_1 : (g'_2, g''_2)-weight = (c0, d0), E'_2(0) x = E''_2(0) x = 0};  X = eigenspace of
  D = o(omega_{S x 1}) = o(omega_{A x 1}) + o(omega_{C x 1}) on Y for the eigenvalue h_z.
  Unknown: Hermitian H on X (<x, y> = x^dagger H y).  Conditions (CKLW invariance, <a_n b, c> = (-1)^d <b, (theta a)_{-n} c>
  for quasi-primary a of weight d, zero modes n = 0):
     weight one (d = 1, theta E = -F, theta H = -H):  O_E^dag H = H O_F,  O_F^dag H = H O_E,  O_H^dag H = H O_H ;
     weight two primaries p_tau (hw, d = 2):  O_p^dag H = c_tau H O_{p_lw},  p_lw = F'^a F''^b p,
        theta p = c_tau p_lw with c_tau = sign(N_tau)/(a! b!), N_tau = (p_lw, p)  (forced by positivity of S and theta^2 = 1);
     omega_{C x 1} (theta-fixed):  O^dag H = H O.
  The solution space (real dimension) and the signature of its generator (exact Hermitian LDL over Q(sqrt(-6))).
Usage: python3 -I r2_signature.py CI   (CI = candidate omega_{C x 1}, 0..3).  Output: r2_CI.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time
from math import factorial
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
S2 = pickle.load(open(os.path.join(CM.HERE, 'S2_%d.pkl' % CI), 'rb'))
u, wC1 = OC['cands'][CI]
omA1 = V.vadd(dict(DG['SUG']['g1p']), DG['SUG']['g1pp'])
Ep, Hp, Fp = CM.SL['g1p']; Epp, Hpp, Fpp = CM.SL['g1pp']
E2p, E2pp = CM.SL['g2p'][0], CM.SL['g2pp'][0]
def qs_neg(v): return Q.vscale(v, (mpq(-1), mpq(0)))
def smul(v): return ({k: -6*c for k, c in v[1].items()}, dict(v[0]))     # s * v
def qcoords(B, w):
    """coordinates (Q(s)) of the Q(s)-vector w in the Q(s)-span of the list B, or None"""
    cols = []
    for b in B: cols += [b, smul(b)]
    keys = sorted(set(k for c in cols + [w] for part in c for k in part), key=repr)
    kidx = {k: i for i, k in enumerate(keys)}
    n = len(cols)
    M = fmpq_mat(2*len(keys), n + 1)
    for j, c in enumerate(cols + [w]):
        for part in (0, 1):
            for k, x in c[part].items(): M[2*kidx[k] + part, j] = CM.to_fq(x)
    Rr, r = M.rref()
    sol = [fmpq(0)]*n
    for i in range(r):
        j = next(j for j in range(n + 1) if Rr[i, j] != 0)
        if j == n: return None
        sol[j] = Rr[i, n]
    return [(CM.to_mq(sol[2*i]), CM.to_mq(sol[2*i + 1])) for i in range(len(B))]
def op_matrix(B, f):
    """matrix (list of rows) of the operator f on span(B) (must be stable)"""
    cols = []
    for b in B:
        c = qcoords(B, Q.vclean(f(b)))
        assert c is not None, 'subspace not stable'
        cols.append(c)
    m = len(B)
    return [[cols[j][i] for j in range(m)] for i in range(m)]
def build_slice(c0, d0, hz):
    Y = []
    for mu in sorted(set(CM.WW.weight_of(k)[:2] for k in V.basis(1) if CM.WW.weight_of(k)[2:] == (c0, d0))):
        Y += CM.cond_kernel(1, mu + (c0, d0), [(E2p, 0), (E2pp, 0)])
    Yq = [Q.vreal(y) for y in Y]
    def Dop(x): return Q.vadd(Q.mode_rq(omA1, 1, x), Q.mode(wC1, 1, x))
    Dm = op_matrix(Yq, Dop)
    # eigenspace for hz
    lam = (mpq(hz.numerator, hz.denominator), mpq(0))
    M = [[Q.sub(Dm[i][j], lam) if i == j else Dm[i][j] for j in range(len(Yq))] for i in range(len(Yq))]
    ker = qkernel(M)
    X = []
    for c in ker:
        v = Q.vzero()
        for y, x in zip(Yq, c):
            if not Q.iszero(x): Q.vadd(v, y, x)
        X.append(Q.vclean(v))
    return Yq, Dm, X
def qkernel(M):
    r, c = len(M), len(M[0])
    A = fmpq_mat(2*r, 2*c)
    for i in range(r):
        for j in range(c):
            a, b = M[i][j]
            A[2*i, 2*j] = CM.to_fq(a); A[2*i, 2*j + 1] = CM.to_fq(-6*b)
            A[2*i + 1, 2*j] = CM.to_fq(b); A[2*i + 1, 2*j + 1] = CM.to_fq(a)
    Rr, rk = A.rref()
    n = 2*c
    piv = [next(j for j in range(n) if Rr[i, j] != 0) for i in range(rk)]
    free = [j for j in range(n) if j not in piv]
    out = []
    for fj in free:
        x = [fmpq(0)]*n; x[fj] = fmpq(1)
        for i, pj in enumerate(piv): x[pj] = -Rr[i, fj]
        out.append([(CM.to_mq(x[2*j]), CM.to_mq(x[2*j + 1])) for j in range(c)])
    # Q(s)-basis
    basis = []
    for v in out:
        if qrank(basis + [v]) > len(basis): basis.append(v)
    return basis
def qrank(vecs):
    if not vecs: return 0
    c = len(vecs[0])
    A = fmpq_mat(2*len(vecs), 2*c)
    for i, v in enumerate(vecs):
        for j, (a, b) in enumerate(v):
            A[2*i, 2*j] = CM.to_fq(a); A[2*i, 2*j + 1] = CM.to_fq(b)
            A[2*i + 1, 2*j] = CM.to_fq(-6*b); A[2*i + 1, 2*j + 1] = CM.to_fq(a)
    return A.rref()[1]//2
def qrank_vecs(vecs):
    """rank over Q(s) of a list of Q(s)-vectors (W-keys)"""
    if not vecs: return 0
    keys = sorted(set(k for v in vecs for part in v for k in part), key=repr)
    kidx = {k: i for i, k in enumerate(keys)}
    A = fmpq_mat(2*len(vecs), 2*len(keys))
    for i, v in enumerate(vecs):
        sv = smul(v)
        for part in (0, 1):
            for k, x in v[part].items(): A[2*i, 2*kidx[k] + part] = CM.to_fq(x)
            for k, x in sv[part].items(): A[2*i + 1, 2*kidx[k] + part] = CM.to_fq(x)
    return A.rref()[1]//2
# ---------------------------------------------------------------- generators and theta
def lower(p, a, b):
    for _ in range(a): p = Q.mode_rq(Fp, 0, p)
    for _ in range(b): p = Q.mode_rq(Fpp, 0, p)
    return p
GENS2 = []
for tau, vecs in sorted(S2['prims'].items()):
    if tau == (0, 0) or not vecs: continue
    p = vecs[0]; a, b = tau
    plw = lower(p, a, b)
    N = Q.vac(Q.mode(plw, 3, p)); assert N[1] == 0 and N[0] != 0
    c = mpq(1 if N[0] > 0 else -1, factorial(a)*factorial(b))
    if os.environ.get('FLIP') == '%d,%d' % tau: c = -c          # NEGATIVE CONTROL: wrong sign of theta on this type
    GENS2.append((tau, p, plw, c, N[0]))
def hermitian_solutions(X, verbose=True):
    m = len(X)
    W1ops = {'E1': (Ep, 0), 'F1': (Fp, 0), 'H1': (Hp, 0), 'E2': (Epp, 0), 'F2': (Fpp, 0), 'H2': (Hpp, 0)}
    O = {nm: op_matrix(X, lambda x, z=z: Q.mode_rq(z[0], 0, x)) for nm, z in W1ops.items()}
    conds = [(O['E1'], (mpq(1), mpq(0)), O['F1']), (O['F1'], (mpq(1), mpq(0)), O['E1']), (O['H1'], (mpq(1), mpq(0)), O['H1']),
             (O['E2'], (mpq(1), mpq(0)), O['F2']), (O['F2'], (mpq(1), mpq(0)), O['E2']), (O['H2'], (mpq(1), mpq(0)), O['H2'])]
    for tau, p, plw, c, N in GENS2:
        Op = op_matrix(X, lambda x, p=p: Q.mode(p, 1, x))
        Olw = op_matrix(X, lambda x, plw=plw: Q.mode(plw, 1, x))
        conds.append((Op, (c, mpq(0)), Olw))
        # theta(p_lw) = c (a! b!)^2 p : equivalent to the adjoint of the previous condition, included as a check
        conds.append((Olw, (c*(factorial(tau[0])*factorial(tau[1]))**2, mpq(0)), Op))
    OC1 = op_matrix(X, lambda x: Q.mode(wC1, 1, x))
    conds.append((OC1, (mpq(1), mpq(0)), OC1))
    # unknowns: H_jj = r_j (real) ; H_jk = a_jk + s b_jk (j < k)
    idx = {}
    for j in range(m): idx[('r', j)] = len(idx)
    for j in range(m):
        for k in range(j + 1, m):
            idx[('a', j, k)] = len(idx); idx[('b', j, k)] = len(idx)
    nv = len(idx)
    def Hent(j, k):
        """H_jk as dict var -> Q(s) coefficient"""
        if j == k: return {idx[('r', j)]: Q.ONE}
        if j < k: return {idx[('a', j, k)]: Q.ONE, idx[('b', j, k)]: Q.S}
        return {idx[('a', k, j)]: Q.ONE, idx[('b', k, j)]: Q.scal(mpq(-1), Q.S)}
    rows = []
    for (A, c, B) in conds:
        # (A^dag H)_{jk} - c (H B)_{jk} = 0 ;  (A^dag)_{jl} = conj(A_{lj})
        for j in range(m):
            for k in range(m):
                eq = {}
                for l in range(m):
                    co = Q.conj(A[l][j])
                    if not Q.iszero(co):
                        for var, x in Hent(l, k).items(): eq[var] = Q.add(eq.get(var, Q.Z), Q.mul(co, x))
                    co2 = B[l][k]
                    if not Q.iszero(co2):
                        for var, x in Hent(j, l).items(): eq[var] = Q.sub(eq.get(var, Q.Z), Q.mul(Q.mul(c, x), co2))
                if any(not Q.iszero(x) for x in eq.values()):
                    rows.append({v: x[0] for v, x in eq.items() if x[0]})
                    rows.append({v: x[1] for v, x in eq.items() if x[1]})
    rows = [r for r in rows if r]
    M = fmpq_mat(max(len(rows), 1), nv)
    for i, r in enumerate(rows):
        for v, x in r.items(): M[i, v] = CM.to_fq(x)
    Rr, rk = M.rref()
    piv = [next(j for j in range(nv) if Rr[i, j] != 0) for i in range(rk)]
    free = [j for j in range(nv) if j not in piv]
    sols = []
    for fj in free:
        x = [fmpq(0)]*nv; x[fj] = fmpq(1)
        for i, pj in enumerate(piv): x[pj] = -Rr[i, fj]
        Hm = [[Q.Z]*m for _ in range(m)]
        for j in range(m):
            for k in range(m):
                t = Q.Z
                for var, co in Hent(j, k).items():
                    if x[var] != 0: t = Q.add(t, Q.scal(CM.to_mq(x[var]), co))
                Hm[j][k] = t
        sols.append(Hm)
    return sols, len(rows)
def inertia(Hm):
    """exact inertia (n+, n-, n0) of a Hermitian matrix over Q(s): congruence H -> P^dag H P (Sylvester)"""
    A = [row[:] for row in Hm]; n = len(A)
    pos = neg = zero = 0
    active = list(range(n))
    while active:
        p = next((i for i in active if not Q.iszero(A[i][i])), None)
        if p is None:
            pair = next(((i, j) for i in active for j in active if i != j and not Q.iszero(A[i][j])), None)
            if pair is None:
                zero += len(active); break
            i, j = pair
            t = Q.conj(A[i][j])            # e_i -> e_i + t e_j : new A_ii = 2 |A_ij|^2 > 0
            for k in range(n): A[k][i] = Q.add(A[k][i], Q.mul(A[k][j], t))
            for k in range(n): A[i][k] = Q.add(A[i][k], Q.mul(Q.conj(t), A[j][k]))
            p = i
        d = A[p][p]; assert d[1] == 0, 'non-real diagonal'
        if d[0] > 0: pos += 1
        else: neg += 1
        active.remove(p)
        for i in active:
            if Q.iszero(A[i][p]): continue
            f = Q.div(A[i][p], d)
            for k in active: A[i][k] = Q.sub(A[i][k], Q.mul(f, A[p][k]))
    return pos, neg, zero
CASES = [  # (label, (c0, d0) of w0 in R_{pi z}, h_z, expected dim R_z)
    ('d4', (2, 0), Fr(5, 13), 5), ('d1', (0, 4), Fr(6, 13), 7), ('d5', (0, 2), Fr(7, 13), 13), ('d6', (0, 0), Fr(8, 13), 14),
    ('c1', (0, 4), Fr(1, 3), 3), ('c2', (0, 2), Fr(2, 3), 27)]
if __name__ == '__main__':
    print('candidate', CI, '; primaries: ', [(t, str(c), str(N)) for t, p, plw, c, N in GENS2], flush=True)
    out = {}
    only = sys.argv[2].split(',') if len(sys.argv) > 2 else None
    for lab, (c0, d0), hz, dimR in CASES:
        if only and lab not in only: continue
        Yq, Dm, X = build_slice(c0, d0, hz)
        print('%s: slice Y dim %d, X = eigenspace of D for h = %s : dim %d (expected %d) (%.0fs)' % (lab, len(Yq), hz, len(X), dimR, time.time() - t0), flush=True)
        sols, nrows = hermitian_solutions(X)
        res = {'dimX': len(X), 'dim_solutions': len(sols), 'equations': nrows}
        if len(sols) == 1:
            ine = inertia(sols[0])
            res['inertia'] = ine
            print('   invariant Hermitian forms: real dimension %d ; inertia (n+, n-, n0) = %s  -> %s' %
                  (len(sols), ine, 'DEFINITE' if (ine[2] == 0 and (ine[0] == 0 or ine[1] == 0)) else 'INDEFINITE/DEGENERATE'), flush=True)
            res['H_diag'] = [Q.s(sols[0][i][i]) for i in range(len(X))]
        else:
            print('   invariant Hermitian forms: real dimension', len(sols), flush=True)
            res['inertias'] = [inertia(Hm) for Hm in sols]
            print('   inertias of basis solutions:', res['inertias'])
        out[lab] = res
    fl = os.environ.get('FLIP')
    json.dump(out, open(os.path.join(CM.HERE, 'r2_%d%s%s.json' % (CI, '' if not only else '_' + '_'.join(only), '' if not fl else '_FLIP' + fl.replace(',', ''))), 'w'), indent=1)
    print('time %.0fs' % (time.time() - t0))
