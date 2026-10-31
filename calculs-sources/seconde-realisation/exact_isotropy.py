"""Axis B3 (exact version): for each candidate C_2, each primary type tau and each eigenvalue xi of x_1 on Pi_tau,
the space E = ker(x_1 - xi) cap ker(e_2) cap ker(x_2) (exact), the restrictions to E of the quadratic forms
Q_y(u) = (u, y_1 u), y in Y = B cap C_2^perp, and of the norm Q_0(u) = (u, u).  The line p_tau must lie in one
such E, be Q_y-isotropic for all y and satisfy Q_0(u) != 0.  We also record h1 = eigenvalue of (omega1)_1 and
h2 = 2 - h_tau - h1 (omega1 = 2f with f = a e + b x, b = 1/sqrt(D), a = (1 - b lambda)/2; both signs of b
correspond to swapping omega1 and omega2)."""
import sys, json, time, pickle
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, cases as C, griess_subalgebras as GS, primaries as PR, isotropy_test as IT
def h_tau(case, tau):
    k1, k2 = {'I': (1, 39), 'II': (3, 13)}[case]
    a, b = tau
    return sp.Rational(a*(a + 2), 4*(k1 + 2)) + sp.Rational(b*(b + 2), 4*(k2 + 2))
def red(z): return sp.expand(sp.radsimp(z))
def common_isotropic(Qs, Q0, d):
    """decide whether there is u != 0 in C^d with Q(u) = 0 for Q in Qs and Q0(u) != 0; returns
    ('none'|'exists'|'undecided', info).  Exact (sympy)."""
    if d == 0: return 'none', 'E = 0'
    Qs = [Q for Q in Qs if not Q.is_zero_matrix]
    if Q0.is_zero_matrix: return 'none', 'Q0 = 0 on E'
    if not Qs: return 'exists', 'all Q_y vanish on E'
    us = sp.symbols('u0:%d' % d)
    uv = sp.Matrix(us)
    eqs = [red((uv.T*Q*uv)[0, 0]) for Q in Qs]
    q0 = red((uv.T*Q0*uv)[0, 0])
    if d == 1:
        return ('exists', 'd=1, Q_y vanish') if all(e == 0 for e in eqs) else ('none', 'd=1, some Q_y != 0')
    # charts u_k = 1, u_j = 0 for j < k ; solve; check Q0 != 0
    sols_found = []
    for k in range(d):
        sub = {us[j]: 0 for j in range(k)}; sub[us[k]] = 1
        E2 = [red(e.subs(sub)) for e in eqs]
        E2 = [e for e in E2 if e != 0]
        vars_ = list(us[k + 1:])
        if any(e.is_number for e in E2): continue
        if not vars_:
            if not E2 and red(q0.subs(sub)) != 0: sols_found.append(sub)
            continue
        if not E2:
            return 'exists', 'positive-dimensional (no equation in chart %d)' % k
        Gb = sp.groebner(E2 + [sp.Symbol('w')*red(q0.subs(sub)) - 1], *(vars_ + [sp.Symbol('w')]), order='grevlex')
        if list(Gb) == [1]:
            continue
        return 'exists', 'Groebner basis != 1 in chart %d (dim>=0 solution with Q0 != 0)' % k
    return ('exists', 'isolated point') if sols_found else ('none', 'Groebner = 1 in all charts')
def run(case, only=None):
    T, basis, prod, gram = GS.griess_data(case)
    D = PR.setup(case)
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras']
    n = len(basis); gramM = sp.Matrix(gram)
    pair = {tau: IT.pairing(D['T'], tau, d['P']) for tau, d in D['types'].items()}
    out = []
    for ci, cand in enumerate(G):
        if only is not None and ci not in only: continue
        xs = [sp.sympify(c) for c in cand['x']]
        xv = sp.Matrix([0] + xs); ev = sp.Matrix([1] + [0]*(n - 1))
        Ycond = sp.Matrix.vstack(ev.T*gramM, xv.T*gramM).applyfunc(red)
        Ys = Ycond.nullspace(simplify=True)
        lam = sp.sympify(cand['lambda']); Dd = sp.sympify(cand['D'])
        b = 1/sp.sqrt(Dd); a = (1 - b*lam)/2
        row = {'c': cand.get('c'), 'types': {}, 'excluded_by': []}
        for tau, d in D['types'].items():
            m = len(d['P'])
            X = sum((xv[i]*d['X1'][i] for i in range(n)), sp.zeros(m, m)).applyfunc(red)
            r2 = d['X2'][0].rows
            Kops = [d['X2'][0], sum((xv[i]*d['X2'][i] for i in range(n)), sp.zeros(r2, m)).applyfunc(red)] if r2 else []
            Yops = [sum((y[i]*d['X1'][i] for i in range(n)), sp.zeros(m, m)).applyfunc(red) for y in Ys]
            Pm = pair[tau]
            t = sp.Symbol('t')
            cp = X.charpoly(t).as_expr()
            fac = sp.factor_list(sp.Poly(cp, t, extension=[sp.sqrt(13)] if any(sp.sympify(c).has(sp.sqrt(13)) for c in cand['x']) else None))[1]
            if any(f.degree() > 1 for f, _ in fac):
                row['types'][str(tau)] = {'status': 'non-split charpoly (handled by the certified Arb test)', 'charpoly': str(sp.factor(cp))}
                continue
            roots = sp.roots(sp.Poly(cp, t), multiple=False)
            if sum(roots.values()) != m:
                row['types'][str(tau)] = {'status': 'non-split charpoly', 'charpoly': str(sp.factor(cp))}
                continue
            ent = []
            for xi in roots:
                M = (X - xi*sp.eye(m))
                for Kb in Kops: M = M.col_join(Kb)
                Eb = M.applyfunc(red).nullspace(simplify=True)
                dE = len(Eb)
                h1 = red(a*(2 - h_tau(case, tau)) + 2*b*xi); h2 = red(2 - h_tau(case, tau) - h1)
                if dE == 0:
                    ent.append({'xi': str(xi), 'dimE': 0, 'h1': str(h1), 'h2': str(h2), 'verdict': 'none', 'why': 'E = 0'}); continue
                Em = sp.Matrix.hstack(*Eb)
                Qs = [(Em.T*Pm*Y*Em).applyfunc(red) for Y in Yops]
                Qs = [((Q + Q.T)/2).applyfunc(red) for Q in Qs]
                Q0 = (Em.T*Pm*Em).applyfunc(red)
                verdict, why = common_isotropic(Qs, Q0, dE)
                ent.append({'xi': str(xi), 'dimE': dE, 'h1': str(h1), 'h2': str(h2), 'verdict': verdict, 'why': why,
                            'rankQ': [Q.rank() for Q in Qs]})
            nposs = sum(1 for z in ent if z['verdict'] == 'exists')
            row['types'][str(tau)] = {'dimPi': m, 'eigen': ent, 'possible_eigenspaces': nposs}
            if nposs == 0 or (PR.PRIMARY_TYPES[case][tau] == 2 and False):
                row['excluded_by'].append(str(tau))
        out.append((ci, row))
    return out
if __name__ == '__main__':
    case = sys.argv[1]; only = [int(a) for a in sys.argv[2:]] or None
    t0 = time.time()
    res = run(case, only)
    for ci, r in res:
        tag = 'EXCLUDED by ' + ', '.join(r['excluded_by']) if r['excluded_by'] else 'survives'
        print('cand %2d c = %s : %s' % (ci, r['c'], tag))
        for tau, v in r['types'].items():
            if 'eigen' in v:
                print('    ', tau, [(z['h1'], z['h2'], z['dimE'], z['verdict']) for z in v['eigen']])
            else:
                print('    ', tau, v)
    print('time %.1fs' % (time.time() - t0))
    json.dump([(ci, r) for ci, r in res], open(f'exact_isotropy_{case}_{"_".join(map(str, only)) if only else "all"}.json', 'w'), indent=1)
