"""Axis B3 (first part): type-by-type necessary conditions on the line p_tau = U_2 cap Pi_tau.
For a candidate C_2 = span(e, x) (Axis A2), let Y = B cap C_2^perp = Com_{M_b}(A)_2 (B = C_2 + Y, U _|_ M_b).
For u in p_tau and y in Y, y_1 u lies in M_b, hence (u, y_1 u) = 0, where (u, v) = <u, theta v> and theta v is
the lowest weight vector f'_0^a f''_0^b v of the g'+g''-module of v.  Also (u, u) != 0, e_2 u = x_2 u = 0, and u
is an eigenvector of x_1.  The eigenvectors of x_1 on Pi_tau are computed with rigorous enclosures (Arb, Rump's
method); an eigenline is EXCLUDED when some condition is certified nonzero.  A candidate is excluded when, for
some type tau, fewer admissible lines remain than the multiplicity of tau in U_2."""
import sys, json, time, itertools
import sympy as sp
from flint import acb_mat, arb, acb, ctx
from gmpy2 import mpq
import lattice_voa as V, cases as C, griess_subalgebras as GS, primaries as PR, linalg as LA
ctx.prec = 320
S13 = arb(13).sqrt()
def to_acb(z):
    z = sp.expand(sp.radsimp(sp.sympify(z)))
    q = z.coeff(sp.sqrt(13)) if z.has(sp.sqrt(13)) else 0
    p = sp.expand(z - q*sp.sqrt(13))
    assert p.is_Rational and sp.sympify(q).is_Rational, ('not in Q(sqrt13)', z)
    p = sp.Rational(p); q = sp.Rational(q)
    return acb(arb(int(p.p))/int(p.q) + S13*arb(int(q.p))/int(q.q))
def acbmat(M):
    return acb_mat(M.rows, M.cols, [to_acb(M[i, j]) for i in range(M.rows) for j in range(M.cols)])
def theta(T, tau, v):
    for _ in range(tau[0]): v = V.clean(V.mode_vec(T['f1'], 0, v))
    for _ in range(tau[1]): v = V.clean(V.mode_vec(T['f2'], 0, v))
    return v
def pairing(T, tau, P):
    th = [theta(T, tau, p) for p in P]
    return sp.Matrix(len(P), len(P), lambda i, j: sp.Rational(int(V.form(P[i], th[j]).numerator), int(V.form(P[i], th[j]).denominator)))
def nonzero(z):          # certified nonzero acb
    return not z.contains(0)
def run(case):
    T, basis, prod, gram = GS.griess_data(case)
    D = PR.setup(case)
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras']
    n = len(basis); gramM = sp.Matrix(gram)
    pair = {tau: pairing(D['T'], tau, d['P']) for tau, d in D['types'].items()}
    out = []
    for ci, cand in enumerate(G):
        xs = [sp.sympify(c) for c in cand['x']]
        xv = sp.Matrix([0] + xs)                     # coordinates in basis [e] + perp
        ev = sp.Matrix([1] + [0]*(n - 1))
        # Y = vectors y with <y, e> = <y, x> = 0
        Ycond = sp.Matrix.vstack((ev.T*gramM), (xv.T*gramM)).applyfunc(lambda t: sp.expand(sp.radsimp(t)))
        Ys = Ycond.nullspace(simplify=True)
        assert len(Ys) == n - 2
        row = {'c': cand.get('c'), 'types': {}, 'excluded_by': []}
        for tau, d in D['types'].items():
            m = len(d['P'])
            X = sum((xv[i]*d['X1'][i] for i in range(n)), sp.zeros(m, m))
            Yops = [sum((y[i]*d['X1'][i] for i in range(n)), sp.zeros(m, m)) for y in Ys]
            r2 = d['X2'][0].rows
            Kops = []
            if r2:
                Kops = [d['X2'][0], sum((xv[i]*d['X2'][i] for i in range(n)), sp.zeros(r2, m))]
            Pm = pair[tau]
            Xa = acbmat(X); Pa = acbmat(Pm); Ya = [acbmat(Y) for Y in Yops]; Ka = [acbmat(K) for K in Kops]
            try:
                E, R = Xa.eig(right=True, algorithm='rump')
            except Exception as exc:
                row['types'][str(tau)] = {'dimPi': m, 'status': 'eig failed: ' + str(exc)[:60]}
                continue
            lines = []
            for j in range(m):
                u = acb_mat(m, 1, [R[i, j] for i in range(m)])
                ut = u.transpose()
                reasons = []
                for yi, Yb in enumerate(Ya):
                    q = (ut*Pa*Yb*u)[0, 0]
                    if nonzero(q): reasons.append('(u,y%d u)' % yi)
                for ki, Kb in enumerate(Ka):
                    w = Kb*u
                    if any(nonzero(w[i, 0]) for i in range(w.nrows())): reasons.append(['e_2 u', 'x_2 u'][ki])
                q0 = (ut*Pa*u)[0, 0]
                lines.append({'eigenvalue': str(E[j])[:40], 'excluded': bool(reasons), 'reasons': reasons,
                              'norm_certified_nonzero': nonzero(q0)})
            nadm = sum(1 for l in lines if not l['excluded'])
            row['types'][str(tau)] = {'dimPi': m, 'admissible_lines': nadm, 'lines': lines}
            if nadm < PR.PRIMARY_TYPES[case][tau]:
                row['excluded_by'].append(str(tau))
        out.append(row)
    return out
if __name__ == '__main__':
    res = {}
    for case in sys.argv[1:] or ['I', 'II']:
        t0 = time.time()
        res[case] = run(case)
        print('case', case)
        for i, r in enumerate(res[case]):
            adm = {t: v.get('admissible_lines', v.get('status')) for t, v in r['types'].items()}
            tag = 'EXCLUDED by ' + ', '.join(r['excluded_by']) if r['excluded_by'] else 'SURVIVES'
            print('  cand %2d c = %-46s %s' % (i, r['c'], tag))
            print('        admissible lines per type:', adm)
        print('  time %.1fs' % (time.time() - t0))
    json.dump(res, open('isotropy_test.json', 'w'), indent=1)
