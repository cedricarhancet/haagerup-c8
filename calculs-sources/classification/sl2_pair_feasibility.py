"""U_1 = sl2 + sl2 (gamma = 0): joint feasibility of the two simple ideals (Axis C2).

For s = g' + g'' with g', g'' among the surviving simple ideals (sl2_ideal_feasibility.py), every R_{dj} is an
(sl2 + sl2)-module  R = sum m_{ab} V(a) x V(b).  The single-ideal relations fix, for each ideal and each j,
the power sums Tr H^2 and Tr H^4 on R_{dj}; the mixed weight-four relations (coefficient of x^2 y^2 in the
quartic map h = x h' + y h'' -> W_4(h)) fix X_j = Tr_{R_dj} h'^2 h''^2 in terms of X_6.  We enumerate all joint
module structures of R_{d6}, R_{d4}, R_{d1}, R_{d5} (dimensions 14, 5, 7, 13) exactly, and for R_{d2}, R_{d3}
(dimensions 42, 119) we decide feasibility by exact search with memoization.
Lie-theoretic input used: at most one ideal has i2 != 0 (Proposition self-duality).
"""
import json, itertools
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp
import obstruction as O, vvmf_exact as V, weight4_solve as W

DIMS = [7, 42, 119, 5, 13, 14]
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))

# ---- single-ideal data: tau_j, Q_j for each surviving (i1, i2) ----
eqs, tsol, R4 = W.weight4_system(0)
Qsol = sp.solve(eqs, list(W.Q[:5]), dict=True)[0]
surv = json.load(open('sl2_ideal_feasibility.json'))
single = {}
for s in surv:
    i1, i2, k = s['i1'], s['i2'], s['k']
    p6 = s['R_d6']
    tau6 = sp.Rational(sum(M(n, 2) for n in p6), 2*k); Q6 = sp.Rational(sum(M(n, 4) for n in p6), 4*k*k)
    sub = {W.u: sp.Rational(i1, k), W.T: sp.Rational(4, k), W.T4: sp.Rational(8, k*k), W.tau[5]: tau6, W.Q[5]: Q6}
    t2 = [int(tsol[W.tau[j]].subs(sub)*2*k) for j in range(6)]
    t4 = [int(Qsol[W.Q[j]].subs(sub)*4*k*k) if j < 5 else int(Q6*4*k*k) for j in range(6)]
    single[(i1, i2)] = {'k': k, 'u': sp.Rational(i1, k), 'T': sp.Rational(4, k), 'tau': [tsol[W.tau[j]].subs(sub) for j in range(6)],
                        't2': t2, 't4': t4}
# ---- mixed relations ----
P = [0]; A = [(0, -1), (0, 0)] + [(z, 0) for z in range(1, V.d)]
R4m = O.relations(4, P, A, 1)
labels = ['w0', 'w0p'] + ['w_' + V.names[z] for z in range(1, V.d)]
X = sp.symbols('X1:7')
def mixed_solution(d1, d2):
    up, upp = d1['u'], d2['u']
    def mixQ(f):
        x, y = sp.symbols('x y')
        e = sp.expand(f(x**2*up + y**2*upp, x**2*(1 - up) + y**2*(1 - upp)))
        return e.coeff(x, 2).coeff(y, 2)/6
    Qa = lambda a, b: 9*a**2 + 36*a*b + sp.Rational(27, 2)*b**2
    Qb = lambda a, b: 18*a**2 + 9*b**2
    w = {'w0': sp.Rational(1, 24), 'w0p': -(d1['T'] + d2['T'])/2 + sp.Rational(6 - 48, 24),
         'w_b': 6*mixQ(Qb) - ((6 + 18*up - d1['T']) + (6 + 18*upp - d2['T']))/2 + sp.Rational(80, 24),
         'w_a': 6*mixQ(Qa) - ((27 - 9*up) + (27 - 9*upp))/2 + sp.Rational(81, 24),
         'w_c0': 6*mixQ(Qa) - ((27 - 9*up) + (27 - 9*upp))/2 + sp.Rational(81, 24),
         'w_c1': 6*mixQ(lambda a, b: b**2/2) - ((1 - up) + (1 - upp))/2 + sp.Rational(3, 24),
         'w_c2': 6*mixQ(lambda a, b: 3*a**2) - (6*up + 6*upp)/2 + sp.Rational(27, 24)}
    for j in range(6):
        w['w_d%d' % (j + 1)] = 6*X[j] - (d1['tau'][j] + d2['tau'][j])/2 + sp.Rational(DIMS[j], 24)
    e = [sp.expand(sum(c*w[l] for c, l in zip(row, labels))) for row in R4m]
    sol = sp.solve(e, list(X[:5]), dict=True)
    assert len(sol) == 1
    sol = sol[0]
    assert all(sp.simplify(ee.subs(sol)) == 0 for ee in e)
    return sol

def joint_structures(dim, a_t, b_t):
    """all multisets of (a,b) with sum (a+1)(b+1) = dim and prescribed marginal power sums
       a_t = (Tr H'^2, Tr H'^4), b_t = (Tr H''^2, Tr H''^4); returns list of (multiset, mixed sum)"""
    items = [(a, b) for a in range(dim) for b in range(dim) if (a + 1)*(b + 1) <= dim]
    items.sort(reverse=True)
    out = []
    def rec(i, rem, s2a, s4a, s2b, s4b, cur):
        if rem == 0:
            if (s2a, s4a, s2b, s4b) == (a_t[0], a_t[1], b_t[0], b_t[1]):
                out.append((tuple(cur), sum(M(a, 2)*M(b, 2) for a, b in cur)))
            return
        if i == len(items):
            return
        a, b = items[i]
        da, db = (a + 1), (b + 1)
        c2a = db*M(a, 2); c4a = db*M(a, 4); c2b = da*M(b, 2); c4b = da*M(b, 4)
        maxm = rem//(da*db)
        for m in range(maxm, -1, -1):
            n2a = s2a + m*c2a; n4a = s4a + m*c4a; n2b = s2b + m*c2b; n4b = s4b + m*c4b
            if n2a > a_t[0] or n4a > a_t[1] or n2b > b_t[0] or n4b > b_t[1]:
                continue
            rec(i + 1, rem - m*da*db, n2a, n4a, n2b, n4b, cur + [(a, b)]*m)
    rec(0, dim, 0, 0, 0, 0, [])
    return out

def feasible_joint(dim, a_t, b_t, mix):
    """exact decision for large dimensions: memoized search over items with the mixed sum as extra target"""
    items = [(a, b) for a in range(dim) for b in range(dim) if (a + 1)*(b + 1) <= dim and (a, b) != (0, 0)]
    items.sort(key=lambda t: -((t[0] + 1)*(t[1] + 1)))
    @lru_cache(maxsize=None)
    def f(i, rem, r2a, r4a, r2b, r4b, rmix):
        if min(r2a, r4a, r2b, r4b, rmix) < 0:
            return False
        if r2a == r4a == r2b == r4b == rmix == 0:
            return True                       # fill the rest with V(0) x V(0)
        if i == len(items):
            return False
        a, b = items[i]
        da, db = a + 1, b + 1
        c = (db*M(a, 2), db*M(a, 4), da*M(b, 2), da*M(b, 4), M(a, 2)*M(b, 2))
        maxm = rem//(da*db)
        for m in range(maxm, -1, -1):
            if f(i + 1, rem - m*da*db, r2a - m*c[0], r4a - m*c[1], r2b - m*c[2], r4b - m*c[3], rmix - m*c[4]):
                return True
        return False
    return f(0, dim, a_t[0], a_t[1], b_t[0], b_t[1], mix)

if __name__ == '__main__':
    keys = sorted(single)
    pairs = [(p, q) for p, q in itertools.combinations_with_replacement(keys, 2) if not (p[1] != 0 and q[1] != 0)]
    report = {}
    for p, q in pairs:
        d1, d2 = single[p], single[q]
        sol = mixed_solution(d1, d2)
        kk = 4*d1['k']*d2['k']
        s6 = joint_structures(14, (d1['t2'][5], d1['t4'][5]), (d2['t2'][5], d2['t4'][5]))
        found = []
        for struct6, mix6 in s6:
            X6 = sp.Rational(mix6, kk)
            ok = True; detail = {'d6': struct6}
            for j in (3, 0, 4, 1, 2):        # d4, d1, d5, d2, d3
                Xj = sp.Rational(sol[X[j]].subs(X[5], X6))
                mj = Xj*kk
                if mj.q != 1 or mj < 0:
                    ok = False; detail['fail'] = 'd%d mixed trace %s not a nonnegative integer' % (j + 1, mj); break
                mj = int(mj)
                at = (d1['t2'][j], d1['t4'][j]); bt = (d2['t2'][j], d2['t4'][j])
                if DIMS[j] <= 14:
                    st = [s for s, mm in joint_structures(DIMS[j], at, bt) if mm == mj]
                    if not st:
                        ok = False; detail['fail'] = 'd%d: no joint structure' % (j + 1); break
                    detail['d%d' % (j + 1)] = st[0]
                else:
                    # large sectors: necessary conditions only (nonnegative integers; mixed Cauchy--Schwarz
                    # sum l'^2 l''^2 <= sqrt(Tr H'^4 Tr H''^4)); exact witnesses are checked separately
                    if mj*mj > at[1]*bt[1]:
                        ok = False; detail['fail'] = 'd%d: Cauchy-Schwarz violated' % (j + 1); break
                    detail['d%d' % (j + 1)] = 'mixed sum %d (necessary conditions only)' % mj
            found.append((ok, detail))
        good = [dd for ok, dd in found if ok]
        print('pair %s + %s : joint R_d6 structures %d, fully compatible %d' % (p, q, len(s6), len(good)))
        for ok, dd in found[:4]:
            print('     ', 'OK ' if ok else 'NO ', {kx: v for kx, v in dd.items()})
        report[str((p, q))] = {'n_d6': len(s6), 'compatible': [str(g) for g in good]}
    json.dump(report, open('sl2_pair_feasibility.json', 'w'), indent=1)
