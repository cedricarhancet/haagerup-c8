"""Axis B3, exact version for candidates defined over Q(sqrt13) (same logic as exact_isotropy.py).
Linear algebra over Q(sqrt13) by the doubled rational representation (flint): z = p + q sqrt13 acts on
Q^2 by [[p, 13q], [q, p]].  Characteristic polynomial of x_1 from the doubled matrix (= chi_X * chi_Xbar),
factored over Q(sqrt13); eigenspaces and kernels computed over Q and reassembled into Q(sqrt13)-bases;
isotropy decided by Groebner bases over QQ<sqrt(13)>."""
import sys, json, time
import sympy as sp
from flint import fmpq_mat, fmpq, acb_mat, acb_poly, acb, arb, ctx
ctx.prec = 400
import griess_subalgebras as GS, primaries as PR, isotropy_test as IT, exact_isotropy as EI
S = sp.sqrt(13)
def split(z):
    z = sp.expand(sp.radsimp(sp.sympify(z)))
    q = z.coeff(S) if z.has(S) else sp.Integer(0)
    p = sp.expand(z - q*S)
    assert p.is_Rational and q.is_Rational, z
    return sp.Rational(p), sp.Rational(q)
def fq(r): r = sp.Rational(r); return fmpq(int(r.p), int(r.q))
def dbl(Mp, Mq):
    r, c = Mp.shape
    M = fmpq_mat(2*r, 2*c)
    for i in range(r):
        for j in range(c):
            a, b = Mp[i, j], Mq[i, j]
            if a: M[i, j] = fq(a); M[i + r, j + c] = fq(a)
            if b: M[i, j + c] = fq(13*b); M[i + r, j] = fq(b)
    return M
def nullspace(M):
    Rr, rank = M.rref(); n = M.ncols(); piv = []
    for i in range(rank): piv.append(next(j for j in range(n) if Rr[i, j] != 0))
    free = [j for j in range(n) if j not in piv]; out = []
    for fj in free:
        x = [fmpq(0)]*n; x[fj] = fmpq(1)
        for i, pj in enumerate(piv): x[pj] = -Rr[i, fj]
        out.append(x)
    return out
def qbasis(ns, m):
    """Q(sqrt13)-basis (sympy column vectors) of the Q-space spanned by ns (closed under sqrt13)"""
    vecs = [sp.Matrix([sp.Rational(int(v[i].p), int(v[i].q)) + S*sp.Rational(int(v[i + m].p), int(v[i + m].q)) for i in range(m)]) for v in ns]
    basis = []; cur = fmpq_mat(0, 2*m)
    def rk(rows):
        if not rows: return 0
        M = fmpq_mat(len(rows), 2*m, [x for r in rows for x in r]); return M.rref()[1]
    rows = []
    for v, raw in zip(vecs, ns):
        J = [raw[i + m]*13 for i in range(m)] + [raw[i] for i in range(m)]
        if rk(rows + [raw, J]) > rk(rows):
            rows += [raw, J]; basis.append(v)
        if 2*len(basis) == len(ns): break
    return basis
def lin_comb(mats, coefs, shape):
    P = sp.zeros(*shape); Q = sp.zeros(*shape)
    for Mx, (p, q) in zip(mats, coefs):
        if p: P += p*Mx
        if q: Q += q*Mx
    return P, Q
def run(case, ci):
    T, basis, prod, gram = GS.griess_data(case)
    D = PR.setup(case)
    cand = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    n = len(basis); gramM = sp.Matrix(gram)
    xs = [split(c) for c in cand['x']]
    xcoef = [(0, 0)] + xs
    xv = sp.Matrix([0] + [sp.sympify(c) for c in cand['x']]); ev = sp.Matrix([1] + [0]*(n - 1))
    Ycond = sp.Matrix.vstack(ev.T*gramM, xv.T*gramM).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
    Ys = Ycond.nullspace(simplify=True)
    lam = sp.sympify(cand['lambda']); Dd = sp.sympify(cand['D'])
    b = 1/sp.sqrt(Dd); a = (1 - b*lam)/2
    row = {'c': cand.get('c'), 'types': {}, 'excluded_by': []}
    t = sp.Symbol('t')
    for tau, d in D['types'].items():
        m = len(d['P'])
        Xp, Xq = lin_comb(d['X1'], xcoef, (m, m))
        r2 = d['X2'][0].rows
        Kp = []; Kq = []
        if r2:
            K2p, K2q = lin_comb(d['X2'], xcoef, (r2, m))
            Kp = [d['X2'][0], K2p]; Kq = [sp.zeros(r2, m), K2q]
        Yops = []
        for y in Ys:
            yc = [split(z) for z in y]
            Yops.append(lin_comb(d['X1'], yc, (m, m)))
        Pm = IT.pairing(D['T'], tau, d['P'])
        M2 = dbl(Xp, Xq)
        cp2 = M2.charpoly()
        coeffs = [sp.Rational(int(cp2[i].p), int(cp2[i].q)) for i in range(cp2.degree() + 1)]
        poly = sp.Poly(sum(c*t**i for i, c in enumerate(coeffs)), t)
        fac = sp.factor_list(poly, extension=[S])[1]
        roots = []
        for f, mult in fac:
            if f.degree() == 1:
                cs = f.all_coeffs(); roots.append(sp.expand(sp.radsimp(-cs[1]/cs[0])))
            else:
                roots.append(None)
        ent = []
        # simple roots of nonlinear factors: certified balls (Arb); eigenvector = a column of adj(X - xi)
        for (f, mult) in fac:
            if f.degree() == 1: continue
            if mult > 1:
                ent.append({'factor': str(f.as_expr()), 'verdict': 'exists', 'why': 'nonlinear factor with multiplicity > 1 (undecided)'}); continue
            cs = [IT.to_acb(cc) for cc in reversed(f.all_coeffs())]
            rts = acb_poly(cs).roots(maxprec=4000)
            Xa = IT.acbmat(Xp + S*Xq); Pa = IT.acbmat(Pm)
            Ya = [IT.acbmat(Yp + S*Yq) for (Yp, Yq) in Yops]
            Ka = [IT.acbmat(kp + S*kq) for kp, kq in zip(Kp, Kq)]
            for rt in rts:
                A = acb_mat(m, m, [Xa[i, j] - (rt if i == j else 0) for i in range(m) for j in range(m)])
                if not A.det().contains(0):
                    continue              # not an eigenvalue of X (root of the conjugate factor)
                u = None
                for col in range(m):
                    cand_u = []
                    for i in range(m):
                        minor = acb_mat(m - 1, m - 1, [A[r, c] for r in range(m) if r != col for c in range(m) if c != i])
                        cand_u.append(minor.det()*(-1)**(i + col))
                    if any(not z.contains(0) for z in cand_u):
                        u = acb_mat(m, 1, cand_u); break
                if u is None:
                    ent.append({'root': str(rt)[:30], 'verdict': 'exists', 'why': 'adjugate not certified nonzero (undecided)'}); continue
                ut = u.transpose(); reasons = []
                for yi, Yb in enumerate(Ya):
                    if not (ut*Pa*Yb*u)[0, 0].contains(0): reasons.append('(u,y%d u)' % yi)
                for ki, Kb in enumerate(Ka):
                    w = Kb*u
                    if any(not w[i, 0].contains(0) for i in range(w.nrows())): reasons.append('K%d' % ki)
                ent.append({'root': str(rt)[:30], 'dimE': 1, 'verdict': 'none' if reasons else 'exists',
                            'why': 'certified (Arb): ' + ', '.join(reasons) if reasons else 'not excluded numerically'})
        seen = set()
        for xi in roots:
            if xi is None or xi in seen: continue
            seen.add(xi)
            xp, xq = split(xi)
            Mp = Xp - xp*sp.eye(m); Mq = Xq - xq*sp.eye(m)
            blocks_p = [Mp] + Kp; blocks_q = [Mq] + Kq
            BigP = sp.Matrix.vstack(*blocks_p); BigQ = sp.Matrix.vstack(*blocks_q)
            ns = nullspace(dbl(BigP, BigQ))
            if not ns:
                continue          # xi is an eigenvalue of Xbar only, or E cap K = 0
            Eb = qbasis(ns, m); dE = len(Eb)
            Em = sp.Matrix.hstack(*Eb)
            Qs = []
            for (Yp, Yq) in Yops:
                Y = Yp + S*Yq
                Q = (Em.T*Pm*Y*Em).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
                Qs.append(((Q + Q.T)/2).applyfunc(sp.expand))
            Q0 = (Em.T*Pm*Em).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
            verdict, why = EI.common_isotropic(Qs, Q0, dE)
            h1 = sp.expand(sp.radsimp(a*(2 - EI.h_tau(case, tau)) + 2*b*xi)); h2 = sp.expand(2 - EI.h_tau(case, tau) - h1)
            ent.append({'xi': str(xi), 'dimE': dE, 'h1': str(h1), 'h2': str(h2), 'verdict': verdict, 'why': why})
        nposs = sum(1 for z in ent if z['verdict'] == 'exists')
        row['types'][str(tau)] = {'dimPi': m, 'eigen': ent, 'possible_eigenspaces': nposs}
        if nposs == 0: row['excluded_by'].append(str(tau))
    return row
if __name__ == '__main__':
    case = sys.argv[1]; ci = int(sys.argv[2])
    t0 = time.time()
    r = run(case, ci)
    tag = 'EXCLUDED by ' + ', '.join(r['excluded_by']) if r['excluded_by'] else 'survives'
    print('cand %2d c = %s : %s' % (ci, r['c'], tag))
    for tau, v in r['types'].items():
        print('    ', tau, [(z.get('h1', z.get('root', z.get('factor'))), z.get('h2', ''), z.get('dimE', '?'), z['verdict']) for z in v['eigen']] if 'eigen' in v else v)
    print('time %.1fs' % (time.time() - t0))
    json.dump(r, open(f'exact_isotropy_s13_{case}_{ci}.json', 'w'), indent=1)
