"""Weight-six relations and the surviving simple ideals (Axis C2, sextic step).

W6(h) = Tr h^6 - (5/4) E2 Tr h^4 + (5/16) E2^2 Tr h^2 - (5/576) E2^3 ch  (<h,h> = 1) lies in M_6(rho, mu)
(Heisenberg recursion, tested in weight6_formula_check.py).  Leading coefficients:
   w0 = -5/576, w0p = T6 - 5 T4/4 + 5 T/16 + (5/576)(72 - dim U_1),
   w_z = S_z - 5 Q_z/4 + 5 tau_z/16 - 5 dim R_z/576,   S_z = Tr_{R_z} h^6.
The sextic traces on R_b, R_a, R_c0, R_c1, R_c2 depend on the embedding (e6 has a sextic invariant) and are
computed from the explicit weights (labels of the neutral element, sl2_e6_superset.json).  For each surviving
ideal the unknowns are S_1..S_6; we solve the exact relations, test consistency, integrality and the joint
existence of sl2-modules R_{dj} with prescribed (dim, Tr H^2, Tr H^4, Tr H^6) and highest weights <= level.
"""
import json, sys
from functools import lru_cache
from fractions import Fraction as Fr
import sympy as sp
import obstruction as O, vvmf_exact as V, e6_tools as E, weight4_solve as W

DIMS = [7, 42, 119, 5, 13, 14]
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
names = V.names; d = V.d
labels = ['w0', 'w0p'] + ['w_' + names[z] for z in range(1, d)]
P = [0]; A = [(0, -1), (0, 0)] + [(z, 0) for z in range(1, d)]

def sector_power_sums(lab, i2, p):
    """Tr H^p on R_b, R_a, R_c0, R_c1, R_c2 for the neutral element H of the ideal (labels lab, sl3 index i2)"""
    Mi = sp.Matrix(E.C).inv()
    c = [Fr(int(sp.fraction(t)[0]), int(sp.fraction(t)[1])) for t in sp.Matrix([lab])*Mi]
    w27 = [E.ev(c, w) for w in E.W27]
    rts = [E.ev(c, r) for r in E.ROOTS]
    w3 = {0: [0, 0, 0], 1: [1, -1, 0], 4: [2, 0, -2]}[i2]
    rts3 = [x - y for x in w3 for y in w3 if True]       # includes 3 zeros (Cartan of sl3 has dim 2: fix below)
    rts3 = [x - y for i, x in enumerate(w3) for j, y in enumerate(w3) if i != j] + [0, 0]
    s_ideal = [2, -2, 0]
    Rb = sum(x**p for x in rts) + sum(x**p for x in rts3) - sum(x**p for x in s_ideal)
    Ra = sum((x + y)**p for x in w27 for y in w3)
    Rc1 = sum(y**p for y in w3); Rc2 = sum(x**p for x in w27)
    return {'b': Rb, 'a': Ra, 'c0': Ra, 'c1': Rc1, 'c2': Rc2}

def single_data():
    eqs4, tsol, R4 = W.weight4_system(0)
    Qsol = sp.solve(eqs4, list(W.Q[:5]), dict=True)[0]
    surv = json.load(open('sl2_ideal_feasibility_level.json'))
    sup = {o['index']: o for o in json.load(open('sl2_e6_superset.json'))}
    out = {}
    for s in surv:
        i1, i2, k = s['i1'], s['i2'], s['k']
        p6 = s['R_d6']
        tau6 = sp.Rational(sum(M(n, 2) for n in p6), 2*k); Q6 = sp.Rational(sum(M(n, 4) for n in p6), 4*k*k)
        sub = {W.u: sp.Rational(i1, k), W.T: sp.Rational(4, k), W.T4: sp.Rational(8, k*k), W.tau[5]: tau6, W.Q[5]: Q6}
        tau = [sp.Rational(tsol[W.tau[j]].subs(sub)) for j in range(6)]
        Q = [sp.Rational(Qsol[W.Q[j]].subs(sub)) for j in range(5)] + [Q6]
        out[(i1, i2)] = {'k': k, 'lab': sup[i1]['labels'], 'tau': tau, 'Q': Q, 'u': sp.Rational(i1, k), 'R_d6': p6}
    return out

@lru_cache(maxsize=None)
def feas3(dim, t2, t4, t6, nmax):
    if t2 == 0:
        return t4 == 0 and t6 == 0
    if dim <= 0 or nmax == 0 or t4 < t2 or t6 < t4:
        return False
    for n in range(min(nmax, dim - 1), 0, -1):
        if M(n, 2) <= t2 and M(n, 4) <= t4 and M(n, 6) <= t6:
            if feas3(dim - n - 1, t2 - M(n, 2), t4 - M(n, 4), t6 - M(n, 6), n):
                return True
    return False
def wit3(dim, t2, t4, t6, nmax):
    if t2 == 0:
        return [0]*dim
    for n in range(min(nmax, dim - 1), 0, -1):
        if M(n, 2) <= t2 and M(n, 4) <= t4 and M(n, 6) <= t6 and feas3(dim - n - 1, t2 - M(n, 2), t4 - M(n, 4), t6 - M(n, 6), n):
            return [n] + wit3(dim - n - 1, t2 - M(n, 2), t4 - M(n, 4), t6 - M(n, 6), n)

if __name__ == '__main__':
    R6 = O.relations(6, P, A, 1)
    F6 = O.exact_forms(6, P, A, 1)
    print('weight six: %d exact relations; exact forms give rank %d (sum %d)' % (len(R6), V.rank(F6, len(A)), len(R6) + V.rank(F6, len(A))))
    json.dump({'labels': labels, 'relations': [[str(c) for c in r] for r in R6]}, open('weight6_relations_exact_gamma0.json', 'w'), indent=1)
    Ssym = sp.symbols('S1:7')
    data = single_data()
    report = {}
    for key, D in sorted(data.items()):
        k = D['k']; n2 = 2*k
        ps = {p: sector_power_sums(D['lab'], key[1], p) for p in (2, 4, 6)}
        # consistency checks with weight two and four data (normalization h = H/sqrt(2k))
        u = D['u']
        assert sp.Rational(ps[2]['c2'], n2) == 6*u and sp.Rational(ps[2]['c1'], n2) == 1 - u
        assert sp.Rational(ps[4]['c2'], n2**2) == 3*u**2
        assert sp.Rational(ps[4]['a'], n2**2) == 9*u**2 + 36*u*(1 - u) + sp.Rational(27, 2)*(1 - u)**2
        T, T4, T6 = sp.Rational(4, k), sp.Rational(8, k*k), sp.Rational(128, n2**3)
        w = {'w0': sp.Rational(-5, 576), 'w0p': T6 - 5*T4/4 + 5*T/16 + sp.Rational(5*(72 - 6), 576)}
        for z, dimz in (('b', 80), ('a', 81), ('c0', 81), ('c1', 3), ('c2', 27)):
            w['w_' + z] = sp.Rational(ps[6][z], n2**3) - sp.Rational(5, 4)*sp.Rational(ps[4][z], n2**2) + sp.Rational(5, 16)*sp.Rational(ps[2][z], n2) - sp.Rational(5*dimz, 576)
        # weight-four consistency of b-sector power sums: Q_b = 18u^2 + 9(1-u)^2 - T4
        assert sp.Rational(ps[4]['b'], n2**2) == 18*u**2 + 9*(1 - u)**2 - T4
        for j in range(6):
            w['w_d%d' % (j + 1)] = Ssym[j] - sp.Rational(5, 4)*D['Q'][j] + sp.Rational(5, 16)*D['tau'][j] - sp.Rational(5*DIMS[j], 576)
        eqs = [sp.expand(sum(c*w[l] for c, l in zip(row, labels))) for row in R6]
        Mx, b = sp.linear_eq_to_matrix(eqs, list(Ssym))
        sol = sp.solve(eqs, list(Ssym), dict=True)
        line = 'ideal (i1,i2)=%s level %d: rank in S = %d; ' % (key, k, Mx.rank())
        if not sol:
            print(line + 'INCONSISTENT (no solution)'); report[str(key)] = 'inconsistent'; continue
        sol = sol[0]
        free = [s for s in Ssym if s not in sol]
        print(line + 'free S: %s' % free)
        if free:
            report[str(key)] = 'free parameters %s' % free
            continue
        t6 = [sol[Ssym[j]]*n2**3 for j in range(6)]
        t2 = [D['tau'][j]*n2 for j in range(6)]; t4 = [D['Q'][j]*n2**2 for j in range(6)]
        print('    Tr H^6 on R_d1..R_d6 =', [str(x) for x in t6])
        ok = True; wits = {}
        for j in range(6):
            if any(sp.Rational(x).q != 1 or x < 0 for x in (t2[j], t4[j], t6[j])):
                ok = False; print('    d%d: non-integral or negative power sum' % (j + 1)); break
            if not feas3(DIMS[j], int(t2[j]), int(t4[j]), int(t6[j]), min(k, DIMS[j] - 1)):
                ok = False; print('    d%d: no sl2-module with (dim, Tr H^2, Tr H^4, Tr H^6) = (%d, %s, %s, %s)' % (j + 1, DIMS[j], t2[j], t4[j], t6[j])); break
            wits['d%d' % (j + 1)] = wit3(DIMS[j], int(t2[j]), int(t4[j]), int(t6[j]), min(k, DIMS[j] - 1))
        print('    =>', 'SURVIVES' if ok else 'EXCLUDED', wits if ok else '')
        report[str(key)] = 'survives' if ok else 'excluded'
    json.dump(report, open('weight6_single.json', 'w'), indent=1)

# ---- second pass: enumerate the module structures of R_d4, R_d5, R_d6 (free sextic traces) ----
def structures(dim, t2, t4, nmax):
    out = []
    def rec(rem, r2, r4, mx, cur):
        if r2 == 0 and r4 == 0:
            out.append(tuple(cur + [0]*rem)); return
        for n in range(min(mx, rem - 1), 0, -1):
            if M(n, 2) <= r2 and M(n, 4) <= r4:
                rec(rem - n - 1, r2 - M(n, 2), r4 - M(n, 4), n, cur + [n])
    rec(dim, t2, t4, nmax, [])
    return out

def second_pass():
    R6 = [[sp.Rational(c) for c in r] for r in json.load(open('weight6_relations_exact_gamma0.json'))['relations']]
    Ssym = sp.symbols('S1:7')
    data = single_data()
    res = {}
    for key, D in sorted(data.items()):
        k = D['k']; n2 = 2*k
        ps = {p: sector_power_sums(D['lab'], key[1], p) for p in (2, 4, 6)}
        T, T4, T6 = sp.Rational(4, k), sp.Rational(8, k*k), sp.Rational(128, n2**3)
        w = {'w0': sp.Rational(-5, 576), 'w0p': T6 - 5*T4/4 + 5*T/16 + sp.Rational(5*66, 576)}
        for z, dimz in (('b', 80), ('a', 81), ('c0', 81), ('c1', 3), ('c2', 27)):
            w['w_' + z] = sp.Rational(ps[6][z], n2**3) - sp.Rational(5, 4)*sp.Rational(ps[4][z], n2**2) + sp.Rational(5, 16)*sp.Rational(ps[2][z], n2) - sp.Rational(5*dimz, 576)
        for j in range(6):
            w['w_d%d' % (j + 1)] = Ssym[j] - sp.Rational(5, 4)*D['Q'][j] + sp.Rational(5, 16)*D['tau'][j] - sp.Rational(5*DIMS[j], 576)
        eqs = [sp.expand(sum(c*w[l] for c, l in zip(row, labels))) for row in R6]
        sol = sp.solve(eqs, list(Ssym[:3]), dict=True)[0]
        t2 = [int(D['tau'][j]*n2) for j in range(6)]; t4 = [int(D['Q'][j]*n2**2) for j in range(6)]
        s4 = structures(5, t2[3], t4[3], min(k, 4)); s5 = structures(13, t2[4], t4[4], min(k, 12))
        s6 = [tuple(D['R_d6'])]
        good = []
        for a4 in s4:
            for a5 in s5:
                for a6 in s6:
                    sub = {Ssym[3]: sp.Rational(sum(M(n, 6) for n in a4), n2**3), Ssym[4]: sp.Rational(sum(M(n, 6) for n in a5), n2**3),
                           Ssym[5]: sp.Rational(sum(M(n, 6) for n in a6), n2**3)}
                    t6 = [sp.Rational(sol[Ssym[j]].subs(sub))*n2**3 for j in range(3)]
                    ok = True; wits = {}
                    for j in range(3):
                        if t6[j].q != 1 or t6[j] < 0 or not feas3(DIMS[j], t2[j], t4[j], int(t6[j]), min(k, DIMS[j] - 1)):
                            ok = False; break
                        wits['d%d' % (j + 1)] = wit3(DIMS[j], t2[j], t4[j], int(t6[j]), min(k, DIMS[j] - 1))
                    if ok:
                        good.append({'d4': a4, 'd5': a5, 'd6': a6, **wits})
        print('ideal %s level %d: structures R_d4 %d, R_d5 %d ; compatible with the sextic relations: %d' % (key, k, len(s4), len(s5), len(good)))
        for g in good[:3]:
            print('    ', {kk: [x for x in v if x] for kk, v in g.items()})
        res[str(key)] = len(good)
    json.dump(res, open('weight6_single_pass2.json', 'w'), indent=1)

if __name__ == '__main__' and '--pass2' in sys.argv:
    second_pass()
