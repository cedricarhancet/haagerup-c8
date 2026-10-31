"""Exclusion of the Evans--Gannon vector gamma = 1 (exact arithmetic).

Input: the weight-two relations (weight2_relations_gamma*.json), the explicit traces on the lowest
spaces of the modules b, a, c0, c1, c2 coming from U < V(E6+A2) < V(E8), and tau_j >= 0.
For h real in a Cartan subalgebra of s = U_1 with <h,h> = 1:
   u = (h1,h1) in [0,1]   (e6 component), T = Tr_s (ad h)^2,
   tau_j = Tr h(0)^2 on the lowest space of M_{dj} (>= 0); for gamma = 1, tau_3 = chi(h)^2 and
   tau_3 = 0 when h lies in [s,s] (the lowest space of M_{d3} is one-dimensional).
Simple ideal g of s with Dynkin indices (i1, i2) in e6 + sl3, level k = i1 + i2:
   T = 2 h^vee/k, u = i1/k.  i2 in {0} unless g = sl2, then i2 in {0, 1, 4} (self-duality of 3).
Steps: (1) feasibility for simple ideals; (2) u-range for central elements (T = 0);
(3) enumeration of all reductive s of dimension 6 + 13 gamma with rank p1(s) <= 4.
Usage: python3 gamma1_exclusion.py
"""
import json
from fractions import Fraction as Fr
from itertools import product
import sympy as sp

u, T = sp.symbols('u T', real=True)
tau = sp.symbols('tau1:7', real=True)

def system(gamma):
    """Return expressions of all tau_j in terms of (u, T, tau5, tau6) (gamma=1) or (u, T, tau6) (gamma=0)."""
    R = json.load(open('weight2_relations_gamma%d.json' % gamma))
    xd = [7 - gamma, 42 + 16*gamma, 1 if gamma == 1 else 119, 5 - 3*gamma, 13 - 5*gamma, 14 + 7*gamma]
    w = {'w0': sp.Rational(-1, 12),
         'w0p': T - sp.Rational(1, 12)*((6 + 13*gamma) - 24),
         'w_b': 24*u + 6*(1 - u) - T - sp.Rational(80 - 13*gamma, 12),
         'w_a': 18*u + 27*(1 - u) - sp.Rational(81, 12),
         'w_c0': 18*u + 27*(1 - u) - sp.Rational(81, 12),
         'w_c1': (1 - u) - sp.Rational(3, 12),
         'w_c2': 6*u - sp.Rational(27, 12)}
    for j in range(6):
        w['w_d%d' % (j + 1)] = tau[j] - sp.Rational(xd[j], 12)
    eqs = [sum(sp.Rational(c)*w[l] for c, l in zip(row, R['labels'])) for row in R['relations']]
    free = [tau[4], tau[5]] if gamma == 1 else [tau[5]]
    solve_for = [t for t in tau if t not in free]
    sol = sp.solve(eqs, solve_for, dict=True)
    assert len(sol) == 1
    sol = sol[0]
    # all equations satisfied identically after substitution
    assert all(sp.simplify(e.subs(sol)) == 0 for e in eqs)
    full = [sp.expand(sol.get(t, t)) for t in tau]
    return full, free

def feasible_interval(exprs, var, subs):
    """exprs: affine expressions in var (after subs); return (lo, hi) of var with all exprs >= 0 and var >= 0."""
    lo, hi = Fr(0), None
    for e in exprs:
        e = sp.expand(e.subs(subs))
        a = e.coeff(var); b = e.subs(var, 0)
        a = Fr(str(a)); b = Fr(str(b))
        if a == 0:
            if b < 0:
                return None
        elif a > 0:
            lo = max(lo, -b/a)
        else:
            hi = -b/a if hi is None else min(hi, -b/a)
    if hi is not None and lo > hi:
        return None
    return (lo, hi)

SIMPLE = {'A1': (2, 3, 1), 'A2': (3, 8, 2), 'B2': (3, 10, 2), 'G2': (4, 14, 2), 'A3': (4, 15, 3)}

def allowed_simple(gamma, full, free):
    out = []
    for name, (hv, dim, rk) in SIMPLE.items():
        i2s = (0, 1, 4) if name == 'A1' else (0,)
        for i2 in i2s:
            for i1 in range(0, 40):
                k = i1 + i2
                if k == 0:
                    continue
                subs = {T: Fr(2*hv, k), u: Fr(i1, k)}
                exprs = [full[j] for j in range(6)]
                if gamma == 1:
                    # tau3 = 0 forces tau5 = T - 9u + 3 (solve tau3 = 0 for tau5); then interval in tau6
                    t5 = sp.solve(sp.Eq(full[2], 0), tau[4])[0]
                    sub2 = dict(subs); sub2[tau[4]] = t5.subs(subs)
                    exprs2 = [e.subs(tau[4], t5) for e in exprs] + [t5]
                    iv = feasible_interval(exprs2, tau[5], subs)
                else:
                    iv = feasible_interval(exprs, tau[5], subs)
                if iv is not None:
                    out.append((name, i1, i2, k, iv))
    return out

def center_u_range(gamma, full, free):
    """Fourier--Motzkin: for T = 0 and u in [0,1], find the set of u admitting tau >= 0 (tau3 free >= 0)."""
    # sample u on a fine rational grid and also report the exact threshold by bisection on rationals
    ok = []
    for num in range(0, 1801):
        uu = Fr(num, 1800)
        if gamma == 1:
            # variables tau5, tau6 >= 0 ; constraints tau1..tau4 >= 0 ; 2D feasibility via vertex enumeration
            cons = [sp.expand(full[j].subs({T: 0, u: uu})) for j in range(4)] + [tau[4], tau[5]]
            feas = False
            lines = []
            for c in cons:
                a = Fr(str(c.coeff(tau[4]))); b = Fr(str(c.coeff(tau[5]))); c0 = Fr(str(c.subs({tau[4]: 0, tau[5]: 0})))
                lines.append((a, b, c0))
            for (a1, b1, c1), (a2, b2, c2) in product(lines, lines):
                det = a1*b2 - a2*b1
                if det == 0:
                    continue
                x = (-c1*b2 + c2*b1)/det; y = (-a1*c2 + a2*c1)/det
                if all(a*x + b*y + c >= 0 for a, b, c in lines):
                    feas = True; break
            ok.append((uu, feas))
        else:
            iv = feasible_interval([full[j] for j in range(6)], tau[5], {T: 0, u: uu})
            ok.append((uu, iv is not None))
    feas_u = [x for x, f in ok if f]
    return (min(feas_u), max(feas_u)) if feas_u else None

if __name__ == '__main__':
    for gamma in (1, 0):
        full, free = system(gamma)
        print('=== gamma =', gamma)
        for j in range(6):
            print('   tau%d =' % (j + 1), full[j])
        al = allowed_simple(gamma, full, free)
        print('allowed simple ideals (type, i1, i2, level, tau6-interval):')
        for a in al:
            print('   ', a[:4], 'tau6 in [%s, %s]' % (a[4][0], a[4][1]))
        ur = center_u_range(gamma, full, free)
        print('central elements: feasible u in [%s, %s] (grid 1/1800)' % ur if ur else 'central elements: none')
        if gamma == 1:
            # exact threshold for central elements: u >= 1/9 (derived by hand, checked here at 1/9 and 1/9 - 1/1800)
            # enumeration of s with dim 19, rank p1(s) <= 4
            e6_ideals = [(a[0], SIMPLE[a[0]][1], SIMPLE[a[0]][2]) for a in al if a[2] == 0]
            diag = [a for a in al if a[2] != 0]
            e6_types = sorted(set(e6_ideals))
            print('e6-ideals allowed:', e6_types, ' diagonal sl2 allowed:', [(a[1], a[2]) for a in diag])
            sols = []
            maxn = 3
            for counts in product(range(maxn + 1), repeat=len(e6_types)):
                dim_e6 = sum(c*t[1] for c, t in zip(counts, e6_types))
                rk_e6 = sum(c*t[2] for c, t in zip(counts, e6_types))
                for delta in (0, 1) if diag else (0,):
                    for zdim in range(0, 20):
                        dim = dim_e6 + 3*delta + zdim
                        # centres: u >= 1/9 > 0 so p1 is injective on z; diagonal sl2 has i1 >= 1
                        rank_p1 = rk_e6 + delta + zdim
                        if dim == 19 and rank_p1 <= 4:
                            sols.append((counts, delta, zdim))
            print('candidate s of dimension 19 with rank p1(s) <= 4:', sols if sols else 'NONE')
