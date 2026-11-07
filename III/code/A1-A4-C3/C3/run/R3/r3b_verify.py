"""R3b verification (exact, directly in W): phi_chi(c mono e^beta) = (-1)^{#mono} chi(beta) conj(c) mono e^{-beta},
chi(beta) = 2^{-(c_8(beta^(1)) + c_8(beta^(2)))} (c_8 = coefficient of the simple root alpha_8 of each E8 copy; r3b_chi.json).
Checks for candidate CI:
 (1) phi_chi is an involution and maps the 12 currents of sl2^4 by E -> -F, H -> -H;
 (2) phi_chi(p_tau) = c_tau p_tau^lw for the 8 primary types (c_tau = sign(N_tau)/(a! b!): the PCT forced by positivity,
     used in r2_*.py), and phi_chi(C_2 x 1) = C_2 x 1;
 (3) phi_chi(S_2 x 1) = S_2 x 1 (all 120 basis vectors) and phi_chi(omega_{C x 1}) = omega_{C x 1};
 (4) the Hermitian form (phi_chi x, y) restricted to the R2 slices R_z x w0 (z = d4, d1, d5, d6, c2) is a positive
     multiple of the form found in r2_signature.py, and is positive definite;
 (5) sample automorphism / positivity checks of phi_chi on W (random basis states, weights <= 2).
Usage: python3 -I r3b_verify.py CI.  Output: r3b_verify_CI.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time, random
from math import factorial
from gmpy2 import mpq
import e8e8_voa as V
import common as CM
import qs6 as Q
import r2_signature as R2
t0 = time.time()
CI = R2.CI
DG = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
S2F = pickle.load(open(os.path.join(CM.HERE, 'S2full_%d.pkl' % CI), 'rb')) if os.path.exists(os.path.join(CM.HERE, 'S2full_%d.pkl' % CI)) else None
_C = {}
def chi(ib):
    if ib not in _C:
        y = V.VECS[ib]
        e = V.coords8(y[:8])[7] + V.coords8(y[8:])[7]
        _C[ib] = mpq(1, 2**e) if e >= 0 else mpq(2**(-e))
    return _C[ib]
def phi_rat(v):
    return {(V.NEGI[b], m): c*chi(b)*(-1)**len(m) for (b, m), c in v.items()}
def phi(v):
    """antilinear: (re, im) -> (phi_rat(re), -phi_rat(im))"""
    return (phi_rat(v[0]), {k: -c for k, c in phi_rat(v[1]).items()})
res = {}
# (1)
ok1 = True
for nm in CM.NAMES:
    E, H, F = CM.SL[nm]
    ok1 &= V.clean(phi_rat(E)) == V.scale(F, -1) and V.clean(phi_rat(H)) == V.scale(H, -1) and V.clean(phi_rat(F)) == V.scale(E, -1)
    ok1 &= V.clean(phi_rat(phi_rat(E))) == E
print('(1) phi_chi(E) = -F, phi_chi(H) = -H, phi_chi(F) = -E, phi_chi^2 = 1 on the 12 currents:', ok1, flush=True)
res['currents'] = bool(ok1)
# (2)
ok2 = True
for tau, p, plw, c, N in R2.GENS2:
    lhs = Q.vclean(phi(p)); rhs = Q.vscale(plw, (c, mpq(0)))
    e = Q.eq(lhs, rhs); ok2 &= e
    print('   phi_chi(p_%s) = c_tau p_lw (c_tau = %s):' % (tau, c), e, flush=True)
res['primaries'] = bool(ok2)
# (3)
if S2F:
    from span import Span
    SP = Span(S2F['B'])
    bad = sum(1 for b in S2F['B'] if SP.coords(Q.vclean(phi(b))) is None)
    print('(3) phi_chi(S_2 x 1) in S_2 x 1: failures %d / %d' % (bad, len(S2F['B'])), flush=True)
    res['S2_stable_failures'] = bad
ok3 = Q.eq(Q.vclean(phi(R2.wC1)), R2.wC1)
print('    phi_chi(omega_{C x 1}) = omega_{C x 1}:', ok3, flush=True)
res['omegaC1_fixed'] = bool(ok3)
# (4)
def herm_form(x, y):      # <x, y> = (phi x, y) for weight-1 Q(s)-vectors (bilinear form via V.form on parts)
    px = phi(x)
    rr = V.form(px[0], y[0]) if px[0] and y[0] else mpq(0)
    ii = V.form(px[1], y[1]) if px[1] and y[1] else mpq(0)
    ri = V.form(px[0], y[1]) if px[0] and y[1] else mpq(0)
    ir = V.form(px[1], y[0]) if px[1] and y[0] else mpq(0)
    return (rr - 6*ii, ri + ir)
res['slices'] = {}
for lab, (c0, d0), hz, dimR in R2.CASES:
    if lab == 'c1': continue
    Yq, Dm, X = R2.build_slice(c0, d0, hz)
    sols, _ = R2.hermitian_solutions(X)
    Hs = sols[0]
    Hc = [[herm_form(X[i], X[j]) for j in range(len(X))] for i in range(len(X))]
    # proportionality
    i0 = next(i for i in range(len(X)) if not Q.iszero(Hs[i][i]))
    lam = Q.div(Hc[i0][i0], Hs[i0][i0])
    prop = all(Hc[i][j] == Q.mul(lam, Hs[i][j]) for i in range(len(X)) for j in range(len(X)))
    ine = R2.inertia(Hc)
    print('(4) %s: (phi_chi x, y) on R_%s x w0 = %s x (R2 form) : %s ; inertia %s' % (lab, lab, Q.s(lam), prop, ine), flush=True)
    res['slices'][lab] = {'proportional': bool(prop), 'factor': Q.s(lam), 'inertia': ine}
# (5) samples
random.seed(5)
B1, B2 = V.basis(1), V.basis(2)
fails = 0; tests = 0
for _ in range(40):
    u = {random.choice(B1 + random.sample(B2, 20)): mpq(1)}; v = {random.choice(B1 + random.sample(B2, 20)): mpq(1)}
    wu, wv = V.weight(next(iter(u))), V.weight(next(iter(v)))
    for n in range(max(0, wu + wv - 3), wu + wv):
        if wu + wv - n - 1 > 2: continue
        a = V.clean(phi_rat(V.mode_vec(u, n, v))); b = V.clean(V.mode_vec(phi_rat(u), n, phi_rat(v))); tests += 1
        if a != b: fails += 1
pos_bad = 0
for k in random.sample(B1, 30) + random.sample(B2, 30):
    if V.form(phi_rat({k: mpq(1)}), {k: mpq(1)}) <= 0: pos_bad += 1
print('(5) phi_chi automorphism samples: %d tests, %d failures ; positivity on 60 basis states: %d failures' % (tests, fails, pos_bad))
res['auto_samples'] = (tests, fails); res['positivity_samples_fail'] = pos_bad
json.dump(res, open(os.path.join(CM.HERE, 'r3b_verify_%d.json' % CI), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
