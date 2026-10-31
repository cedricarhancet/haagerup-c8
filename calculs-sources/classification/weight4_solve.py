"""Insert the explicit data into the exact weight-four relations (Axis B2), gamma = 0 and gamma = 1.

W4(h) = Tr h(0)^4 q^{L0-1/3} - (E2/2) Tr h(0)^2 q^{L0-1/3} + (E2^2/48) ch  is in M_4(rho, mu) (<h,h> = 1).
Leading coefficients:  w0 = 1/48 ; w0p = T4 - T/2 + (dim U_1 - 48)/48 ; w_z = Q_z - tau_z/2 + dim R_z/48,
with Q_z = Tr_{R_z} h^4, tau_z = Tr_{R_z} h^2, T = Tr_s (ad h)^2, T4 = Tr_s (ad h)^4, u = <h1,h1>.
Known quartic traces (e6_tools.py and sl3):  Tr_27 h1^4 = 3u^2, Tr_e6 (ad h1)^4 = 18u^2,
Tr_3 h2^4 = (1-u)^2/2, Tr_sl3 (ad h2)^4 = 9(1-u)^2.
"""
import json, sympy as sp
import obstruction as O, vvmf_exact as V

u, T, T4 = sp.symbols('u T T4', real=True)
tau = sp.symbols('tau1:7', real=True)
Q = sp.symbols('Q1:7', real=True)
names = V.names; d = V.d
labels = ['w0', 'w0p'] + ['w_' + names[z] for z in range(1, d)]

def weight2_tau(gamma):
    R = json.load(open('weight2_relations_exact_gamma%d.json' % gamma))
    xd = [7 - gamma, 42 + 16*gamma, 1 if gamma == 1 else 119, 5 - 3*gamma, 13 - 5*gamma, 14 + 7*gamma]
    w = {'w0': sp.Rational(-1, 12), 'w0p': T - sp.Rational(1, 12)*((6 + 13*gamma) - 24),
         'w_b': 6 + 18*u - T - sp.Rational(80 - 13*gamma, 12), 'w_a': 27 - 9*u - sp.Rational(81, 12),
         'w_c0': 27 - 9*u - sp.Rational(81, 12), 'w_c1': 1 - u - sp.Rational(3, 12), 'w_c2': 6*u - sp.Rational(27, 12)}
    for j in range(6):
        w['w_d%d' % (j + 1)] = tau[j] - sp.Rational(xd[j], 12)
    eqs = [sum(sp.Rational(c)*w[l] for c, l in zip(row, R['labels'])) for row in R['relations']]
    free = [tau[4], tau[5]] if gamma == 1 else [tau[5]]
    sol = sp.solve(eqs, [t for t in tau if t not in free], dict=True)[0]
    return {t: sp.expand(sol.get(t, t)) for t in tau}, xd

def weight4_system(gamma):
    P = [0, names.index('d3')] if gamma == 1 else [0]
    A = [(0, -1), (0, 0)] + [(z, -1 if z in P else 0) for z in range(1, d)]
    R4 = O.relations(4, P, A, 1)
    tsol, xd = weight2_tau(gamma)
    taub = 6 + 18*u - T; taua = 27 - 9*u; tc1 = 1 - u; tc2 = 6*u
    Qb = 18*u**2 + 9*(1 - u)**2 - T4
    Qa = 9*u**2 + 36*u*(1 - u) + sp.Rational(27, 2)*(1 - u)**2
    Qc1 = (1 - u)**2/2; Qc2 = 3*u**2
    dimb = 80 - 13*gamma
    w = {'w0': sp.Rational(1, 48), 'w0p': T4 - T/2 + sp.Rational(6 + 13*gamma - 48, 48),
         'w_b': Qb - taub/2 + sp.Rational(dimb, 48), 'w_a': Qa - taua/2 + sp.Rational(81, 48),
         'w_c0': Qa - taua/2 + sp.Rational(81, 48), 'w_c1': Qc1 - tc1/2 + sp.Rational(3, 48),
         'w_c2': Qc2 - tc2/2 + sp.Rational(27, 48)}
    for j in range(6):
        w['w_d%d' % (j + 1)] = Q[j] - tsol[tau[j]]/2 + sp.Rational(xd[j], 48)
    eqs = [sp.expand(sum(c*w[l] for c, l in zip(row, labels))) for row in R4]
    return eqs, tsol, R4

if __name__ == '__main__':
    out = {}
    for gamma in (0, 1):
        eqs, tsol, R4 = weight4_system(gamma)
        unknowns = list(Q) + [T4]
        M, b = sp.linear_eq_to_matrix(eqs, unknowns)
        print('=== gamma =', gamma, ': %d weight-four relations, rank in the unknowns (Q1..Q6, T4) = %d' % (len(eqs), M.rank()))
        sol = sp.solve(eqs, unknowns, dict=True)
        if sol:
            sol = sol[0]
            for k in unknowns:
                print('   ', k, '=', sp.factor(sol.get(k, k)) if k in sol else '(free)')
        # consistency conditions: eliminate the unknowns
        aug = M.row_join(b)
        # left kernel of M gives conditions on b
        ker = M.T.nullspace()
        conds = [sp.factor(sp.expand((k.T*b)[0])) for k in ker]
        print('   consistency conditions (must vanish):', conds)
        out[gamma] = {'solution': {str(k): str(v) for k, v in sol.items()} if sol else None,
                      'conditions': [str(c) for c in conds]}
    json.dump(out, open('weight4_solve_output.json', 'w'), indent=1)
