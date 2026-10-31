"""Axis B3/B4: for a candidate C_2 over Q(sqrt13) (or Q), the admissible eigenspaces E_tau (exact) and, when
dim E_tau = 2 (or 1), the finite set of lines u in P(E_tau) with Q_y(u) = 0 for all y in Y and (u,u) != 0.
The lines are the common roots of the binary quadratic forms Q_y restricted to E_tau (gcd over Q(sqrt13))."""
import sys, json, time, pickle
import sympy as sp
from flint import fmpq_mat, fmpq
import lattice_voa as V, griess_subalgebras as GS, primaries as PR, isotropy_test as IT, exact_isotropy as EI
import exact_isotropy_s13 as ES
S = sp.sqrt(13)
def admissible_spaces(case, ci):
    T, basis, prod, gram = GS.griess_data(case)
    D = PR.setup(case)
    cand = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    n = len(basis); gramM = sp.Matrix(gram)
    xs = [ES.split(c) for c in cand['x']]; xcoef = [(0, 0)] + xs
    xv = sp.Matrix([0] + [sp.sympify(c) for c in cand['x']]); ev = sp.Matrix([1] + [0]*(n - 1))
    Ycond = sp.Matrix.vstack(ev.T*gramM, xv.T*gramM).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
    Ys = Ycond.nullspace(simplify=True)
    lam = sp.sympify(cand['lambda']); Dd = sp.sympify(cand['D']); b = 1/sp.sqrt(Dd); a = (1 - b*lam)/2
    t = sp.Symbol('t'); out = {}
    for tau, d in D['types'].items():
        m = len(d['P'])
        Xp, Xq = ES.lin_comb(d['X1'], xcoef, (m, m))
        r2 = d['X2'][0].rows; Kp = []; Kq = []
        if r2:
            K2p, K2q = ES.lin_comb(d['X2'], xcoef, (r2, m)); Kp = [d['X2'][0], K2p]; Kq = [sp.zeros(r2, m), K2q]
        Yops = [ES.lin_comb(d['X1'], [ES.split(z) for z in y], (m, m)) for y in Ys]
        Pm = IT.pairing(D['T'], tau, d['P'])
        cp2 = ES.dbl(Xp, Xq).charpoly()
        poly = sp.Poly(sum(sp.Rational(int(cp2[i].p), int(cp2[i].q))*t**i for i in range(cp2.degree() + 1)), t)
        fac = sp.factor_list(poly, extension=[S])[1]
        spaces = []
        for f, mult in fac:
            if f.degree() != 1: continue
            cs = f.all_coeffs(); xi = sp.expand(sp.radsimp(-cs[1]/cs[0]))
            xp, xq = ES.split(xi)
            BigP = sp.Matrix.vstack(Xp - xp*sp.eye(m), *Kp); BigQ = sp.Matrix.vstack(Xq - xq*sp.eye(m), *Kq)
            ns = ES.nullspace(ES.dbl(BigP, BigQ))
            if not ns: continue
            Eb = ES.qbasis(ns, m); Em = sp.Matrix.hstack(*Eb)
            Qs = []
            for (Yp, Yq) in Yops:
                Q = (Em.T*Pm*(Yp + S*Yq)*Em).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
                Qs.append(((Q + Q.T)/2).applyfunc(sp.expand))
            Q0 = (Em.T*Pm*Em).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
            verdict, why = EI.common_isotropic(Qs, Q0, len(Eb))
            if verdict != 'exists': continue
            h1 = sp.nsimplify(sp.N(a*(2 - EI.h_tau(case, tau)) + 2*b*xi, 60), rational=True, tolerance=1e-40)
            spaces.append({'xi': xi, 'E': Em, 'Qs': Qs, 'Q0': Q0, 'h1': h1, 'h2': 2 - EI.h_tau(case, tau) - h1})
        out[tau] = {'P': d['P'], 'spaces': spaces, 'pair': Pm}
    return D, out
def lines_in(space):
    """finite list of isotropic lines (coefficient vectors in the basis E) when dim E <= 2; else None"""
    E = space['E']; d = E.shape[1]
    if d == 1:
        return [sp.Matrix([1])] if all(Q[0, 0] == 0 for Q in space['Qs']) and space['Q0'][0, 0] != 0 else []
    if d != 2: return None
    s = sp.Symbol('s')
    polys = []
    for Q in space['Qs']:
        p = sp.Poly(sp.expand(Q[0, 0] + 2*Q[0, 1]*s + Q[1, 1]*s**2), s, extension=[S])
        if not p.is_zero: polys.append(p)
    res = []
    # line u = (1, s)
    if not polys:
        return 'all'
    g = polys[0]
    for p in polys[1:]: g = sp.gcd(g, p)
    if g.degree() >= 1:
        for r in sp.roots(g, multiple=True):
            u = sp.Matrix([1, r])
            if sp.simplify((u.T*space['Q0']*u)[0, 0]) != 0: res.append(u)
    # line u = (0, 1)
    u = sp.Matrix([0, 1])
    if all(sp.simplify((u.T*Q*u)[0, 0]) == 0 for Q in space['Qs']) and sp.simplify((u.T*space['Q0']*u)[0, 0]) != 0:
        res.append(u)
    return res
if __name__ == '__main__':
    case = sys.argv[1]; ci = int(sys.argv[2])
    t0 = time.time()
    D, out = admissible_spaces(case, ci)
    summary = {}
    for tau, v in out.items():
        info = []
        for sp_ in v['spaces']:
            L = lines_in(sp_)
            info.append({'dimE': sp_['E'].shape[1], 'h1': str(sp_['h1']), 'h2': str(sp_['h2']),
                         'lines': 'positive-dimensional' if L is None else (L if L == 'all' else [str(list(l)) for l in L])})
        summary[str(tau)] = info
        print(tau, [(z['dimE'], z['h1'], z['h2'], z['lines'] if isinstance(z['lines'], str) else len(z['lines'])) for z in info])
    json.dump(summary, open(f'isotropic_lines_{case}_{ci}.json', 'w'), indent=1)
    print('time %.1fs' % (time.time() - t0))
