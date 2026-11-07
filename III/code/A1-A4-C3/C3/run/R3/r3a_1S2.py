"""R3a step 3d (symmetric control): 1 x S_2 primaries inside W (exact over Q(sqrt(-6))), for candidate CI:
P'_tau = A x A-primary vectors of weight (0, 0, a, b), hw for (g'_2, g''_2); the 1 x S_2 primaries are the eigenvectors of
o(omega_{1 x C}) = o(omega_CC - u_CI) with eigenvalue 2 - h_A(a, b).  Checks: multiplicities = U2 primaries, commutant
with S x 1 (omega_{C x 1}_(n) p = 0, n = 1, 2, 3), phi_chi(p) = sign(N)/(a! b!) p_lw, and the 165 coupling invariants
(with the lowering word f'_2^a f''_2^b) compared with Article I.  Output: S2bis_CI.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time
from math import factorial
from gmpy2 import mpq
import sympy as sp
import e8e8_voa as V
import common as CM
import qs6 as Q
import r2_signature as R2
from r3b_verify_core import phi
t0 = time.time()
CI = R2.CI
DG = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
om1C = Q.vclean(Q.vadd(Q.vreal(DG['omCC']), R2.wC1, (mpq(-1), mpq(0))))
E2p, H2p, F2p = CM.SL['g2p']; E2pp, H2pp, F2pp = CM.SL['g2pp']
ops = [(CM.E_OPS[i], 0) for i in range(4)] + [(x, n) for x in CM.CUR for n in (1, 2)]
prims = {}
for tau in [(0, 4), (1, 1), (1, 3), (1, 5), (1, 7), (2, 2), (2, 6), (3, 3), (0, 2), (2, 0), (0, 6), (2, 4), (3, 1)]:
    a, b = tau
    P = CM.cond_kernel(2, (0, 0, a, b), ops)
    if not P: prims[tau] = []; continue
    Pq = [Q.vreal(p) for p in P]
    M = R2.op_matrix(Pq, lambda x: Q.mode(om1C, 1, x))
    lam = 2 - CM.hA(a, b)
    L = (mpq(lam.numerator, lam.denominator), mpq(0))
    ker = R2.qkernel([[Q.sub(M[i][j], L) if i == j else M[i][j] for j in range(len(P))] for i in range(len(P))])
    vecs = []
    for c in ker:
        v = Q.vzero()
        for p, x in zip(Pq, c):
            if not Q.iszero(x): Q.vadd(v, p, x)
        vecs.append(Q.vclean(v))
    prims[tau] = vecs
    print('type', tau, ': dim P\' =', len(P), '; 1 x S_2 primaries:', len(vecs), flush=True)
bad = sum(1 for vv in prims.values() for p in vv for n in (1, 2, 3) if not Q.iszero_v(Q.mode(R2.wC1, n, p)))
print('omega_{C x 1}_(n) p = 0 (n = 1,2,3): failures', bad)
types = sorted(t for t, v in prims.items() if v)
MT = {}; N = {}; okphi = True
for tau in types:
    a, b = tau; p = prims[tau][0]; mods = {}; x = p
    for i in range(a + 1):
        y = x
        for j in range(b + 1):
            mods[(i, j)] = y; y = Q.mode_rq(F2pp, 0, y)
        x = Q.mode_rq(F2p, 0, x)
    MT[tau] = mods
    N[tau] = Q.vac(Q.mode(mods[(a, b)], 3, p))
    c = mpq(1 if N[tau][0] > 0 else -1, factorial(a)*factorial(b))
    e = Q.eq(Q.vclean(phi(p)), Q.vscale(mods[(a, b)], (c, mpq(0))))
    okphi &= e and N[tau][1] == 0
print('phi_chi(p) = sign(N)/(a!b!) p_lw for all 1 x S_2 primaries (N real):', okphi)
res = {}
for tau in types:
    for sigma in types:
        for rho in types:
            m1, m2 = rho[0] - tau[0], rho[1] - tau[1]; a, b = sigma
            if (a - m1) % 2 or (b - m2) % 2 or abs(m1) > a or abs(m2) > b: continue
            x = Q.mode(MT[tau][(0, 0)], 1, MT[sigma][((a - m1)//2, (b - m2)//2)])
            T = Q.vac(Q.mode(MT[rho][(rho[0], rho[1])], 3, x))
            if Q.iszero(T): continue
            res[(tau, sigma, rho)] = Q.div(Q.mul(T, T), Q.mul(Q.mul(N[tau], N[sigma]), N[rho]))
A = json.load(open(os.path.join(CM.HERE, '..', 'data', 'invariants_K_0.json')))
art = {(eval(r['tau']), eval(r['sigma']), eval(r['rho'])): [sp.Rational(x) for x in r['I']] for r in A['invariants']}
match = sum(1 for k, I in res.items() if k in art and I[1] == 0 and art[k][1:] == [0, 0, 0] and sp.Rational(str(I[0])) == art[k][0])
print('coupling invariants of 1 x S: %d nonzero, %d in Article I, %d equal exactly' % (len(res), len(art), match))
json.dump({'CI': CI, 'prim_counts': {str(t): len(v) for t, v in prims.items()}, 'commutant_failures': bad, 'phi_ok': bool(okphi),
           'invariants_here': len(res), 'invariants_art': len(art), 'match': match}, open(os.path.join(CM.HERE, 'S2bis_%d.json' % CI), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
