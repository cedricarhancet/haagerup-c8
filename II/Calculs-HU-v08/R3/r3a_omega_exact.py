"""R3a step 2d: exact candidates for omega_{C x 1} over Q(sqrt(-6)), recognized from the numerical solutions
(r3a_idem_complex.py) and VERIFIED EXACTLY: in G-coordinates (u.u = 2u, (u,u) = 9/5, trace and determinant data on
W_1) and directly in W (u_(1)u = 2u, u_(3)u = 9/5, (omega_CC - u)_(1) u = 0).  Output: omegaC1.pkl"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time
from fractions import Fraction as Fr
from gmpy2 import mpq
import numpy as np
import e8e8_voa as V
import common as CM
import qs6 as Q
t0 = time.time()
DG = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
D = pickle.load(open(os.path.join(CM.HERE, 'G_griess.pkl'), 'rb'))
G = DG['G']; prod, form, Mk, comp = D['prod'], D['form'], D['Mk'], D['comp']
n = len(G)
num = json.load(open(os.path.join(CM.HERE, 'idem_complex.json')))
def recog(x):
    f = Fr(x).limit_denominator(10**6)
    assert abs(float(f) - x) < 1e-12, x
    return mpq(f.numerator, f.denominator)
cands = []
for re, im in zip(num['re'], num['im']):
    u = [(recog(a), recog(b/6**0.5)) for a, b in zip(re, im)]
    cands.append(u)
def gmul(x, y):
    out = [Q.Z]*n
    for i in range(n):
        if Q.iszero(x[i]): continue
        for j in range(n):
            if Q.iszero(y[j]): continue
            c = Q.mul(x[i], y[j])
            for k in range(n):
                if prod[(i, j)][k]: out[k] = Q.add(out[k], Q.scal(prod[(i, j)][k], c))
    return out
def gform(x, y):
    t = Q.Z
    for i in range(n):
        for j in range(n):
            if form[i][j]: t = Q.add(t, Q.scal(form[i][j], Q.mul(x[i], y[j])))
    return t
def mat(u, kappa):
    k = len(comp[kappa])
    M = [[Q.Z]*k for _ in range(k)]
    for g in range(n):
        if Q.iszero(u[g]): continue
        for i in range(k):
            for j in range(k):
                if Mk[kappa][g][i][j]: M[i][j] = Q.add(M[i][j], Q.scal(Mk[kappa][g][i][j], u[g]))
    return M
def charpoly_ok(M, lams):
    """check prod (M - lam) = 0 exactly"""
    k = len(M)
    P = [[Q.ONE if i == j else Q.Z for j in range(k)] for i in range(k)]
    for lam in lams:
        L = (mpq(lam.numerator, lam.denominator), mpq(0))
        Ml = [[Q.sub(M[i][j], L) if i == j else M[i][j] for j in range(k)] for i in range(k)]
        P = [[sum_q([Q.mul(P[i][l], Ml[l][j]) for l in range(k)]) for j in range(k)] for i in range(k)]
    return all(Q.iszero(P[i][j]) for i in range(k) for j in range(k))
def sum_q(xs):
    t = Q.Z
    for x in xs: t = Q.add(t, x)
    return t
results = []
omCC = DG['omCC']
for ci, u in enumerate(cands):
    ok_idem = gmul(u, u) == [Q.scal(mpq(2), x) for x in u]
    nrm = gform(u, u)
    ok_spec = all(charpoly_ok(mat(u, kappa), [CM.omega_C1_eig(kappa, s) for s in comp[kappa]]) for kappa in comp)
    # in W
    w = Q.vzero()
    for g, c in zip(G, u):
        if not Q.iszero(c): Q.vadd(w, Q.vreal(g), c)
    w = Q.vclean(w)
    ww = Q.mode(w, 1, w)
    okW = Q.eq(ww, Q.vscale(w, (mpq(2), mpq(0))))
    c2 = Q.vac(Q.mode(w, 3, w))
    om1C = Q.vclean(Q.vadd(Q.vreal(omCC), w, (mpq(-1), mpq(0))))
    orth = Q.iszero_v(Q.mode(w, 1, om1C))
    print('candidate', ci, ': u.u = 2u', ok_idem, '; (u,u) =', Q.s(nrm), '; spectra of o(u) on all W_1 components = h_z - h_A :',
          ok_spec, '; in W: u_(1)u = 2u', okW, ', u_(3)u =', Q.s(c2), ', omega_C1 _(1) omega_1C = 0 :', orth, flush=True)
    results.append({'coords': [(str(a), str(b)) for a, b in u], 'ok': bool(ok_idem and ok_spec and okW and orth)})
    cands[ci] = (u, w)
# relations between candidates
conj = lambda u: [Q.conj(x) for x in u]
for a in range(len(cands)):
    for b in range(a + 1, len(cands)):
        print('  candidates', a, b, ': complex conjugate', conj(cands[a][0]) == cands[b][0])
pickle.dump({'cands': cands, 'results': results}, open(os.path.join(CM.HERE, 'omegaC1.pkl'), 'wb'))
json.dump(results, open(os.path.join(CM.HERE, 'omegaC1.json'), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
