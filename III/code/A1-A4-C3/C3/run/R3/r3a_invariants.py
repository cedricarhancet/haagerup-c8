"""R3a step 3b: identification of S x 1 (inside W) with the algebra S of Article I in weight two, through the
scale-free coupling invariants of Article I (Calculs/generation-forte/invariants_K.py):
   N_tau = (u_tau, f'^a f''^b u_tau),  T = ((u_tau)_(1) w', f'^.. u_rho),  I(tau, sigma, rho) = T^2 / (N_tau N_sigma N_rho),
w' = f'^i f''^j u_sigma of weight rho - tau.  These are invariant under rescaling of the primary lines and under
isomorphisms (inner on S_1).  Exact over Q(sqrt(-6)); compared with invariants_K_0.json (exact, Article I).
Usage: python3 -I r3a_invariants.py CI.  Output: invariants_CI.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time
from fractions import Fraction as Fr
from gmpy2 import mpq
import sympy as sp
import e8e8_voa as V
import common as CM
import qs6 as Q
t0 = time.time()
CI = int(sys.argv[1]) if len(sys.argv) > 1 else 0
S2 = pickle.load(open(os.path.join(CM.HERE, 'S2_%d.pkl' % CI), 'rb'))
Fp, Fpp = CM.SL['g1p'][2], CM.SL['g1pp'][2]
types = sorted(t for t, v in S2['prims'].items() if v and t != (0, 0))
MT = {}
for tau in types:
    a, b = tau; p = S2['prims'][tau][0]
    mods = {}
    x = p
    for i in range(a + 1):
        y = x
        for j in range(b + 1):
            mods[(i, j)] = y
            y = Q.mode_rq(Fpp, 0, y)
        x = Q.mode_rq(Fp, 0, x)
    MT[tau] = mods
N = {tau: Q.vac(Q.mode(MT[tau][(tau[0], tau[1])], 3, MT[tau][(0, 0)])) for tau in types}
print('N_tau:', {str(t): Q.s(x) for t, x in N.items()}, flush=True)
res = {}
for tau in types:
    for sigma in types:
        for rho in types:
            m1, m2 = rho[0] - tau[0], rho[1] - tau[1]
            a, b = sigma
            if (a - m1) % 2 or (b - m2) % 2 or abs(m1) > a or abs(m2) > b: continue
            w = MT[sigma][((a - m1)//2, (b - m2)//2)]
            x = Q.mode(MT[tau][(0, 0)], 1, w)
            T = Q.vac(Q.mode(MT[rho][(rho[0], rho[1])], 3, x))
            if Q.iszero(T): continue
            I = Q.div(Q.mul(T, T), Q.mul(Q.mul(N[tau], N[sigma]), N[rho]))
            res[(tau, sigma, rho)] = I
    print('  tau', tau, '(%.0fs)' % (time.time() - t0), flush=True)
# compare with Article I
A = json.load(open(os.path.join(CM.HERE, '..', 'data', 'invariants_K_0.json')))
art = {}
for r in A['invariants']:
    key = (eval(r['tau']), eval(r['sigma']), eval(r['rho']))
    art[key] = [sp.Rational(x) for x in r['I']]
match = mism = 0; notreal = 0
missing_here = [k for k in art if k not in res]
missing_art = [k for k in res if k not in art]
for k, I in res.items():
    if I[1] != 0: notreal += 1
    if k in art:
        a = art[k]
        ok = (a[1] == 0 and a[2] == 0 and a[3] == 0 and sp.Rational(str(I[0])) == a[0] and I[1] == 0)
        if ok: match += 1
        else:
            mism += 1
            if mism <= 10: print('   MISMATCH', k, Q.s(I), [str(x) for x in a])
print('nonzero couplings here: %d ; in Article I: %d ; equal (exactly): %d ; different: %d ; non-real here: %d' %
      (len(res), len(art), match, mism, notreal))
print('present in Article I but zero here:', missing_here[:10], '; present here but absent in Article I:', missing_art[:10])
json.dump({'CI': CI, 'N': {str(t): Q.s(x) for t, x in N.items()}, 'n_here': len(res), 'n_art': len(art), 'match': match,
           'mismatch': mism, 'nonreal': notreal, 'I': {str(k): Q.s(v) for k, v in res.items()}},
          open(os.path.join(CM.HERE, 'invariants_%d.json' % CI), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
