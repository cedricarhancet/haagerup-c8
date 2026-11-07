"""R2 (general version, multiplicities allowed, top levels in W_1 or W_2): invariant Hermitian forms on R_z.
Slice X = R_z x w0 = eigenspace of D = o(omega_{S x 1}) for h_z in Y = {x in W_k : (g'_2,g''_2)-weight (c0,d0),
E'_2(0) x = E''_2(0) x = 0}, computed weight space by weight space for (g'_1, g''_1).
By S_1-invariance (weight-one conditions E^dag = F, H^dag = H), H = sum_j H_std(j) (x) K_j, K_j Hermitian on the
multiplicity space of the type j (basis: highest weight vectors h_{j,alpha} of X, H(h_{j,g}, h_{j,b}) = K_j[g][b]).
Weight-two conditions (CKLW, d = 2): for each S_2 primary p (type tau, theta p = c_tau p_lw), each pair of types
(i, j) with mu = lambda_j - lambda_i a weight of V(tau), p^(mu) = F'^k1 F''^k2 p, theta p^(mu) = (-1)^{k1+k2} c_tau E'^k1 E''^k2 p_lw:
   H(o(p^(mu)) h_{i,a}, h_{j,b}) = H(h_{i,a}, o(theta p^(mu)) h_{j,b})     (hw-to-hw matrix elements suffice, CG argument),
plus omega_{C x 1} (theta-fixed).  Projections onto hw spaces along F'X + F''X (exact over Q(sqrt(-6))).
Usage: python3 -I r2_reduced.py CI label[,label...]   labels: d1..d6, c1, c2, b.  Output: r2red_CI_labels.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time, itertools
from math import factorial
from fractions import Fraction as Fr
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import e8e8_voa as V
import common as CM
import qs6 as Q
import r2_signature as R2      # reuses qcoords, qkernel, qrank, inertia, GENS2, wC1, omA1 (same candidate CI)
t0 = time.time()
CI = R2.CI
T12 = json.load(open(os.path.join(CM.HERE, '..', 'data', 't12_table.json')))
Ep, Hp, Fp = CM.SL['g1p']; Epp, Hpp, Fpp = CM.SL['g1pp']
E2p, E2pp = CM.SL['g2p'][0], CM.SL['g2pp'][0]
CASES = {  # label: (weight k of the top level R_z x R_{pi z} in W, (c0, d0) type of w0 in R_{pi z}, h_z)
    'd4': (1, (2, 0), Fr(5, 13)), 'd1': (1, (0, 4), Fr(6, 13)), 'd5': (1, (0, 2), Fr(7, 13)), 'd6': (1, (0, 0), Fr(8, 13)),
    'c1': (1, (0, 4), Fr(1, 3)), 'c2': (1, (0, 2), Fr(2, 3)),
    'd2': (2, (3, 3), Fr(11, 13)), 'd3': (2, (1, 5), Fr(15, 13)), 'b': (2, (3, 1), Fr(1))}
def wmode(g, x):
    """zero mode of a W-vector g (rational dict or Q(s)-pair) of weight 1 or 2 on the Q(s)-vector x"""
    if isinstance(g, tuple): return Q.mode(g, 1, x)
    wt = V.weight(next(iter(g)))
    return Q.mode_rq(g, wt - 1, x)
class Slice:
    def __init__(self, k, c0, d0, hz):
        self.k, self.c0, self.d0, self.hz = k, c0, d0, hz
        self.X = {}
    def Xw(self, lam):
        lam = tuple(lam)
        if lam in self.X: return self.X[lam]
        mu = lam + (self.c0, self.d0)
        if not CM.wspace(self.k, mu):
            self.X[lam] = []; return []
        Y = CM.cond_kernel(self.k, mu, [(E2p, 0), (E2pp, 0)])
        if not Y:
            self.X[lam] = []; return []
        Yq = [Q.vreal(y) for y in Y]
        def Dop(x): return Q.vadd(Q.mode_rq(R2.omA1, 1, x), Q.mode(R2.wC1, 1, x))
        Dm = R2.op_matrix(Yq, Dop)
        L = (mpq(self.hz.numerator, self.hz.denominator), mpq(0))
        M = [[Q.sub(Dm[i][j], L) if i == j else Dm[i][j] for j in range(len(Yq))] for i in range(len(Yq))]
        ker = R2.qkernel(M)
        X = []
        for c in ker:
            v = Q.vzero()
            for y, x in zip(Yq, c):
                if not Q.iszero(x): Q.vadd(v, y, x)
            X.append(Q.vclean(v))
        self.X[lam] = X
        return X
def lin_indep(vecs):
    out = []
    for v in vecs:
        if Q.iszero_v(v): continue
        if R2.qrank_vecs(out + [v]) > len(out): out.append(v)
    return out
def qrank_vecs(vecs):
    if not vecs: return 0
    keys = sorted(set(k for v in vecs for part in v for k in part), key=repr)
    kidx = {k: i for i, k in enumerate(keys)}
    A = fmpq_mat(2*len(vecs), 2*len(keys))
    for i, v in enumerate(vecs):
        sv = R2.smul(v)
        for part in (0, 1):
            for k, x in v[part].items(): A[2*i, 2*kidx[k] + part] = CM.to_fq(x)
            for k, x in sv[part].items(): A[2*i + 1, 2*kidx[k] + part] = CM.to_fq(x)
    return A.rref()[1]//2
R2.qrank_vecs = qrank_vecs
def run(label):
    k, (c0, d0), hz = CASES[label]
    Rtypes = {tuple(map(int, t.split(','))): m for t, m in T12[label]['R'].items()}
    SL_ = Slice(k, c0, d0, hz)
    dimX = sum(len(SL_.Xw((a1, a2))) for a1 in range(-6, 7) for a2 in range(-14, 15))
    hw = {}; proj = {}
    for (a, b), mult in sorted(Rtypes.items()):
        Xl = SL_.Xw((a, b))
        # hw vectors: kernel of E', E'' on X_lambda
        imgsE = [(Q.mode_rq(Ep, 0, x), Q.mode_rq(Epp, 0, x)) for x in Xl]
        # coordinates of images in X_{lambda+(2,0)}, X_{lambda+(0,2)}
        Xa, Xb = SL_.Xw((a + 2, b)), SL_.Xw((a, b + 2))
        rows = []
        for (ea, eb) in imgsE:
            ca = R2.qcoords(Xa, Q.vclean(ea)) if Xa else []
            cb = R2.qcoords(Xb, Q.vclean(eb)) if Xb else []
            assert (Xa and ca is not None) or (not Xa and Q.iszero_v(ea)), 'E image outside X'
            assert (Xb and cb is not None) or (not Xb and Q.iszero_v(eb)), 'E image outside X'
            rows.append((ca or []) + (cb or []))
        if rows and rows[0]:
            M = [[rows[j][i] for j in range(len(Xl))] for i in range(len(rows[0]))]
            ker = R2.qkernel(M)
        else:
            ker = [[Q.ONE if i == j else Q.Z for i in range(len(Xl))] for j in range(len(Xl))]
        H = []
        for c in ker:
            v = Q.vzero()
            for x, cc in zip(Xl, c):
                if not Q.iszero(cc): Q.vadd(v, x, cc)
            H.append(Q.vclean(v))
        assert len(H) == mult, (label, (a, b), len(H), mult)
        hw[(a, b)] = H
        # adapted basis of X_lambda: hw + F' X_{+(2,0)} + F'' X_{+(0,2)}
        desc = [Q.vclean(Q.mode_rq(Fp, 0, x)) for x in Xa] + [Q.vclean(Q.mode_rq(Fpp, 0, x)) for x in Xb]
        desc = lin_indep(desc)
        assert len(desc) + len(H) == len(Xl), 'adapted basis has wrong size'
        proj[(a, b)] = H + desc
    print('%s: X = R_%s x w0, dim %d (expected %d) ; types %s (%.0fs)' % (label, label, dimX, T12[label]['dimR'], dict(Rtypes), time.time() - t0), flush=True)
    # unknowns
    idx = {}
    for t, H in hw.items():
        for g in range(len(H)):
            idx[(t, 'r', g)] = len(idx)
            for bb in range(g + 1, len(H)):
                idx[(t, 'a', g, bb)] = len(idx); idx[(t, 'b', g, bb)] = len(idx)
    nv = len(idx)
    def Kent(t, g, bb):
        if g == bb: return {idx[(t, 'r', g)]: Q.ONE}
        if g < bb: return {idx[(t, 'a', g, bb)]: Q.ONE, idx[(t, 'b', g, bb)]: Q.S}
        return {idx[(t, 'a', bb, g)]: Q.ONE, idx[(t, 'b', bb, g)]: Q.scal(mpq(-1), Q.S)}
    def hwcoords(t, v):
        c = R2.qcoords(proj[t], Q.vclean(v))
        assert c is not None, 'image outside X'
        return c[:len(hw[t])]
    gens = [(tau, p, plw, c) for tau, p, plw, c, N in R2.GENS2]
    gens.append(((0, 0), R2.wC1, R2.wC1, mpq(1)))
    rows = []
    def lowerw(v, k1, k2, ops):
        for _ in range(k1): v = Q.mode_rq(ops[0], 0, v)
        for _ in range(k2): v = Q.mode_rq(ops[1], 0, v)
        return v
    for tau, p, plw, c in gens:
        for ti, tj in itertools.product(hw, hw):
            mu = (tj[0] - ti[0], tj[1] - ti[1])
            if abs(mu[0]) > tau[0] or abs(mu[1]) > tau[1] or (tau[0] - mu[0]) % 2 or (tau[1] - mu[1]) % 2: continue
            k1, k2 = (tau[0] - mu[0])//2, (tau[1] - mu[1])//2
            pmu = lowerw(p, k1, k2, (Fp, Fpp))
            tpmu = Q.vscale(lowerw(plw, k1, k2, (Ep, Epp)), (c*(-1)**(k1 + k2), mpq(0)))
            for al, hi in enumerate(hw[ti]):
                x = hwcoords(tj, wmode(pmu, hi))
                for be, hj in enumerate(hw[tj]):
                    y = hwcoords(ti, wmode(tpmu, hj))
                    eq = {}
                    for g, xg in enumerate(x):
                        co = Q.conj(xg)
                        if Q.iszero(co): continue
                        for var, z in Kent(tj, g, be).items(): eq[var] = Q.add(eq.get(var, Q.Z), Q.mul(co, z))
                    for g, yg in enumerate(y):
                        if Q.iszero(yg): continue
                        for var, z in Kent(ti, al, g).items(): eq[var] = Q.sub(eq.get(var, Q.Z), Q.mul(yg, z))
                    r0 = {v_: z[0] for v_, z in eq.items() if z[0]}; r1 = {v_: z[1] for v_, z in eq.items() if z[1]}
                    if r0: rows.append(r0)
                    if r1: rows.append(r1)
    M = fmpq_mat(max(len(rows), 1), nv)
    for i, r in enumerate(rows):
        for v_, x in r.items(): M[i, v_] = CM.to_fq(x)
    Rr, rk = M.rref()
    piv = [next(j for j in range(nv) if Rr[i, j] != 0) for i in range(rk)]
    free = [j for j in range(nv) if j not in piv]
    sols = []
    for fj in free:
        xx = [fmpq(0)]*nv; xx[fj] = fmpq(1)
        for i, pj in enumerate(piv): xx[pj] = -Rr[i, fj]
        Ks = {}
        for t, H in hw.items():
            m_ = len(H)
            Ks[t] = [[Q.Z]*m_ for _ in range(m_)]
            for g in range(m_):
                for bb in range(m_):
                    tot = Q.Z
                    for var, z in Kent(t, g, bb).items():
                        if xx[var] != 0: tot = Q.add(tot, Q.scal(CM.to_mq(xx[var]), z))
                    Ks[t][g][bb] = tot
        sols.append(Ks)
    res = {'dimX': dimX, 'unknowns': nv, 'equations': len(rows), 'dim_solutions': len(sols)}
    if len(sols) == 1:
        Ks = sols[0]
        tot = [0, 0, 0]; per = {}
        for t, K in Ks.items():
            ine = R2.inertia(K); d = (t[0] + 1)*(t[1] + 1)
            per[str(t)] = ine
            for i in range(3): tot[i] += d*ine[i]
        res['inertia'] = tot; res['per_type'] = per
        verdict = 'DEFINITE' if tot[2] == 0 and (tot[0] == 0 or tot[1] == 0) else 'INDEFINITE/DEGENERATE'
        print('   unknowns %d, equations %d ; solution space real dim 1 ; inertia on R_%s = %s -> %s ; per type %s' %
              (nv, len(rows), label, tuple(tot), verdict, per), flush=True)
    else:
        print('   unknowns %d, equations %d ; solution space real dim %d' % (nv, len(rows), len(sols)), flush=True)
    return res
if __name__ == '__main__':
    labels = sys.argv[2].split(',')
    out = {}
    for lab in labels:
        out[lab] = run(lab)
    fl = os.environ.get('FLIP')
    json.dump(out, open(os.path.join(CM.HERE, 'r2red_%d_%s%s.json' % (CI, '_'.join(labels), '' if not fl else '_FLIP' + fl.replace(',', ''))), 'w'), indent=1)
    print('time %.0fs' % (time.time() - t0))
