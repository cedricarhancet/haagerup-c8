"""Simple ideals of s = U_1 for gamma = 0 (Axis C1-C2): exact finite search.

A simple ideal g of s is sl2 (dim s = 6).  Its neutral element H acts on every lowest space R_z with integer
eigenvalues; R_z is a finite-dimensional sl2-module.  With h = H/sqrt(2k), k = i1 + i2:
   u = i1/k, T = 4/k, T4 = Tr_s (ad h)^4 = 8/k^2,
   tau_j = Tr_{R_dj} H^2 / (2k),  Q_j = Tr_{R_dj} H^4 / (2k)^2 .
Weight two (exact relations): tau_1..tau_5 are affine in (u, T, tau_6).  Weight four (exact relations):
Q_1..Q_5 are affine in (u, T, T4, tau_6, Q_6).  For every sl2-module structure of R_{d6} (dimension 14)
we compute (tau_6, Q_6), deduce (tau_j, Q_j), and test whether an sl2-module of dimension dim R_{dj}
with these power sums exists (exact dynamic programming).  i1 runs over a computed SUPERSET of the Dynkin
indices of sl2 in e6 (sl2_e6_superset.py), i2 in {0, 1, 4}.
Level bound (option --level): every sl2 highest weight b occurring in a lowest space R_z satisfies b <= k.
Indeed a highest weight vector v in R_z is annihilated by the positive modes x(n), n >= 1 (R_z is the lowest
L_0-eigenspace) and by e(0); it is a singular vector of the integrable level-k module M_z, hence b <= k.
This is applied to R_{dj} and, using the explicit embedding (sl2_e6_superset.json), to R_b, R_a, R_c0, R_c1, R_c2.
Output: the list of surviving (i1, i2) with the compatible module structures (as multisets of highest weights).
"""
import json, sys
from fractions import Fraction as Fr
from functools import lru_cache
import sympy as sp
import weight4_solve as W

GAMMA = 0
DIMS = [7, 42, 119, 5, 13, 14]
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))

def partitions_sl2(dim):
    """all sl2-modules of the given dimension, as tuples of highest weights (non-increasing)"""
    out = []
    def rec(rem, maxn, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        for n in range(min(maxn, rem - 1), -1, -1):
            cur.append(n); rec(rem - n - 1, n, cur); cur.pop()
    rec(dim, dim - 1, [])
    return out

@lru_cache(maxsize=None)
def feasible(dim, t2, t4, nmax):
    """is there an sl2-module of dimension dim with sum m^2 = t2, sum m^4 = t4, highest weights <= nmax ?"""
    if dim == 0:
        return t2 == 0 and t4 == 0
    if t2 < 0 or t4 < 0:
        return False
    if t2 == 0:
        return t4 == 0           # all trivial
    if nmax == 0:
        return False
    # bounds: each irrep V(n) with n >= 1 has M2 >= 2; t4 >= t2 (integer eigenvalues)
    if t4 < t2:
        return False
    for n in range(nmax, 0, -1):
        if n + 1 <= dim and M(n, 2) <= t2 and M(n, 4) <= t4:
            if feasible(dim - n - 1, t2 - M(n, 2), t4 - M(n, 4), n):
                return True
    return False

def witness(dim, t2, t4, nmax):
    if dim == 0:
        return []
    if t2 == 0 and t4 == 0:
        return [0]*dim
    for n in range(nmax, 0, -1):
        if n + 1 <= dim and M(n, 2) <= t2 and M(n, 4) <= t4 and feasible(dim - n - 1, t2 - M(n, 2), t4 - M(n, 4), n):
            return [n] + witness(dim - n - 1, t2 - M(n, 2), t4 - M(n, 4), n)
    return None

LEVEL = '--level' in sys.argv

def embedding_ok(entry, i2, k):
    """level bound on the explicitly known lowest spaces R_b, R_a, R_c0, R_c1, R_c2"""
    m27 = max(int(x) for x in entry['dec27'])
    m3 = {0: 0, 1: 1, 4: 2}[i2]
    ad = {int(x): v for x, v in entry['dec_ad'].items()}
    if ad.get(2, 0) == 0:           # the image of the ideal is a copy of V(2) in e6: spurious labelling
        return False
    adsl3 = {0: {0: 8}, 1: {2: 1, 1: 2, 0: 1}, 4: {4: 1, 2: 1}}[i2]
    Rb = dict(ad)
    for x, v in adsl3.items():
        Rb[x] = Rb.get(x, 0) + v
    Rb[2] -= 1                      # remove the ideal itself (one copy of the adjoint)
    mb = max(x for x, v in Rb.items() if v > 0)
    return max(m27, m3, m27 + m3, mb) <= k

def main():
    eqs, tsol, R4 = W.weight4_system(GAMMA)
    unknowns = list(W.Q[:5])
    Qsol = sp.solve(eqs, unknowns, dict=True)[0]
    tau_f = [sp.lambdify((W.u, W.T, W.tau[5]), tsol[W.tau[j]]) for j in range(6)]
    sup = json.load(open('sl2_e6_superset.json'))
    i1s = sorted(set(o['index'] for o in sup))
    entry_of = {o['index']: o for o in sup}
    assert len(entry_of) == len(sup)
    P14 = partitions_sl2(14)
    survivors = []
    for i1 in i1s:
        for i2 in (0, 1, 4):
            k = i1 + i2
            u = Fr(i1, k); T = Fr(4, k); T4 = Fr(8, k*k)
            if T > 9*u - 5:
                continue
            if LEVEL and not embedding_ok(entry_of[i1], i2, k):
                continue
            for p in P14:
                if LEVEL and max(p) > k:
                    continue
                t2_6 = sum(M(n, 2) for n in p); t4_6 = sum(M(n, 4) for n in p)
                tau6 = Fr(t2_6, 2*k); Q6 = Fr(t4_6, 4*k*k)
                sub = {W.u: sp.Rational(u.numerator, u.denominator), W.T: sp.Rational(T.numerator, T.denominator),
                       W.T4: sp.Rational(T4.numerator, T4.denominator), W.tau[5]: sp.Rational(tau6.numerator, tau6.denominator),
                       W.Q[5]: sp.Rational(Q6.numerator, Q6.denominator)}
                ok = True; mods = {6: p}
                for j in range(5):
                    tj = sp.Rational(tsol[W.tau[j]].subs(sub)); Qj = sp.Rational(Qsol[W.Q[j]].subs(sub))
                    t2 = tj*2*k; t4 = Qj*4*k*k
                    if t2.q != 1 or t4.q != 1 or t2 < 0 or t4 < 0:
                        ok = False; break
                    t2 = int(t2); t4 = int(t4)
                    nmax = min(DIMS[j] - 1, k) if LEVEL else DIMS[j] - 1
                    if not feasible(DIMS[j], t2, t4, nmax):
                        ok = False; break
                    mods[j + 1] = tuple(witness(DIMS[j], t2, t4, nmax))
                if ok:
                    survivors.append({'i1': i1, 'i2': i2, 'k': k, 'R_d6': list(p),
                                      'witnesses': {('d%d' % j): list(v) for j, v in mods.items()}})
    pairs = sorted(set((s['i1'], s['i2']) for s in survivors))
    print('gamma = 0, level bound %s: surviving simple ideals (i1, i2):' % ('ON' if LEVEL else 'OFF'), pairs)
    for pr in pairs:
        ss = [s for s in survivors if (s['i1'], s['i2']) == pr]
        print('  (i1, i2) = %s, level %d: %d compatible R_d6 structures' % (pr, ss[0]['k'], len(ss)))
        for s in ss[:6]:
            print('      R_d6 =', s['R_d6'], ' e.g.', {kk: v for kk, v in s['witnesses'].items() if kk != 'd6'})
    json.dump(survivors, open('sl2_ideal_feasibility%s.json' % ('_level' if LEVEL else ''), 'w'), indent=0)

if __name__ == '__main__':
    main()
