# Step 3b: polynomial system for u = omega_{C x 1} in G
import sys, pickle, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import sympy as sp
t0 = time.time()
G = pickle.load(open('G.pkl', 'rb'))
D = pickle.load(open('griess.pkl', 'rb'))
gram, P, Mt = D['gram'], D['P'], D['Mt']
n = 16
pivots, Gb, omCC = G['pivots'], G['Gb'], G['omCC']
wCC = [omCC.get(p, Fr(0)) for p in pivots]
Q = lambda f: sp.Rational(f.numerator, f.denominator)
def hA(a, b):  # first factor: levels 3 and 13
    return Fr(a*(a+2), 20) + Fr(b*(b+2), 60)
# check o(omega_CC) = (1 - h_AA) Id on each HW_tau
ok = True
for tau, mats in Mt.items():
    m = len(mats[0])
    M = [[sum(wCC[i]*mats[i][r][c] for i in range(n)) for c in range(m)] for r in range(m)]
    lam = 1 - hA(tau[0], tau[1]) - hA(tau[2], tau[3])
    nSS = 1 if tau in [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)] else 0
    tq = sp.Symbol('t')
    MM = sp.Matrix(m, m, lambda r, c: sp.Rational(M[r][c].numerator, M[r][c].denominator))
    cpol = (tq*sp.eye(m) - MM).det()
    tgt = tq**nSS * (tq - sp.Rational(lam.numerator, lam.denominator))**(m - nSS)
    ok &= sp.expand(cpol - tgt) == 0
    ok &= MM.is_diagonalizable()
print('o(omega_CC): eigenvalue 1 - h_AxA(tau) on non-SS sectors, 0 on (SxS)_1, diagonalizable, for every HW_tau :', ok)
# check omega_CC is idempotent-like: omega_(1)omega = 2 omega in G coordinates
w2 = [sum(wCC[i]*wCC[j]*P[i][j][k] for i in range(n) for j in range(n)) for k in range(n)]
print('omega_CC_(1)omega_CC = 2 omega_CC (G coords):', w2 == [2*x for x in wCC])
print('(omega_CC, omega_CC) =', sum(wCC[i]*wCC[j]*gram[i][j] for i in range(n) for j in range(n)))
# eigenvalue data (sector z, h_z, A-type of first factor)
hz = {'c1': Fr(1, 3), 'c2': Fr(2, 3), 'd1': Fr(6, 13), 'd5': Fr(7, 13), 'd4': Fr(5, 13), 'd6': Fr(8, 13)}
pi = {'c1': 'c2', 'c2': 'c1', 'd1': 'd5', 'd5': 'd1', 'd4': 'd6', 'd6': 'd4'}
import json, os
T = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 't12_table.json')))
for z in hz:
    assert Fr(T[z]['h']) == hz[z]
expected = {}
for z in hz:
    for k1, m1 in T[z]['R'].items():
        a, b = map(int, k1.split(','))
        for k2, m2 in T[pi[z]]['R'].items():
            c, d = map(int, k2.split(','))
            for _ in range(m1*m2):
                expected.setdefault((a, b, c, d), []).append(hz[z] - hA(a, b))
for tau in [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)]:
    expected.setdefault(tau, []).append(Fr(0))
print('types:', len(expected), ' multiplicities match HW dims:', all(len(expected[t]) == len(Mt[t][0]) for t in Mt) and set(expected) == set(Mt))
for t in sorted(expected): print('  ', t, [str(x) for x in expected[t]])
xs = sp.symbols('x0:16')
eqs_lin, eqs_nl = [], []
tt = sp.Symbol('t')
for tau in sorted(Mt):
    mats = Mt[tau]; m = len(mats[0])
    M = sp.Matrix(m, m, lambda r, c: sum(Q(mats[i][r][c])*xs[i] for i in range(n)))
    cp = (tt*sp.eye(m) - M).det(method='berkowitz')
    target = sp.prod([tt - Q(l) for l in expected[tau]])
    diff = sp.Poly(sp.expand(cp - target), tt)
    for co in diff.all_coeffs():
        co = sp.expand(co)
        if co == 0: continue
        if sp.Poly(co, *xs).total_degree() <= 1:
            eqs_lin.append(co)
        else:
            eqs_nl.append(co)
u = sp.Matrix(xs)
for k in range(n):
    e = sum(Q(P[i][j][k])*xs[i]*xs[j] for i in range(n) for j in range(n)) - 2*xs[k]
    eqs_nl.append(sp.expand(e))
eqs_nl.append(sp.expand(sum(Q(gram[i][j])*xs[i]*xs[j] for i in range(n) for j in range(n)) - sp.Rational(9, 5)))
print('linear equations:', len(eqs_lin), ' nonlinear:', len(eqs_nl))
A, bb = sp.linear_eq_to_matrix(eqs_lin, xs)
print('rank of linear system:', A.rank(), ' -> free parameters:', n - A.rank())
pickle.dump({'lin': eqs_lin, 'nl': eqs_nl, 'xs': xs, 'expected': expected}, open('system.pkl', 'wb'))
print('time %.1fs' % (time.time()-t0))
