"""Checks for axis E (exact, over Q):
 (i)  phi_0 : e^lam (x) h-monomial  ->  (-1)^{#factors} e^{-lam} (x) same monomial  is a VOA automorphism of V' on the
      rational form (hence phi_0 o complex conjugation is an antilinear automorphism): phi(u_n v) = phi(u)_n phi(v)
      for all basis states u of weight <= 2, v of weight <= 2 (sampled) and all modes n landing in weight >= 0;
 (ii) H_0(u, v) = (u, phi_0 v) is positive definite on V'_n, n <= 3 (block by lattice vector beta: the Gram matrix of
      the monomials of e^beta is checked by exact Cholesky / leading minors)."""
import random, time, json
from gmpy2 import mpq
import lattice_voa as V
from flint import fmpq_mat
def phi0(v):
    return {(V.neg(l), m): c*(-1)**len(m) for (l, m), c in v.items()}
random.seed(1)
t0 = time.time()
res = {}
B1, B2 = V.basis(1), V.basis(2)
fails = 0; tests = 0
for u in random.sample(B1, 40) + random.sample(B2, 25):
    for v in random.sample(B1, 10) + random.sample(B2, 10):
        wu, wv = V.weight(u), V.weight(v)
        for n in range(-1, wu + wv):
            a = V.clean(V.mode_vec({u: mpq(1)}, n, {v: mpq(1)}))
            b = V.clean(V.mode_vec(phi0({u: mpq(1)}), n, phi0({v: mpq(1)})))
            tests += 1
            if phi0(a) != b: fails += 1
print('(i) automorphism tests:', tests, 'failures:', fails, '(%.0fs)' % (time.time() - t0), flush=True)
res['automorphism_tests'] = tests; res['automorphism_failures'] = fails
for wt in (1, 2, 3):
    keys = V.basis(wt); bybeta = {}
    for k in keys: bybeta.setdefault(k[0], []).append(k)
    worst = None; ok = True; offblock = 0
    for beta, ks in bybeta.items():
        n = len(ks)
        M = fmpq_mat(n, n)
        for i in range(n):
            for j in range(n):
                f = V.form({ks[i]: mpq(1)}, phi0({ks[j]: mpq(1)}))
                M[i, j] = int(f.numerator) if f.denominator == 1 else fmpq_mat([[f]])[0, 0]
        # leading principal minors positive
        for m in range(1, n + 1):
            d = fmpq_mat([[M[i, j] for j in range(m)] for i in range(m)]).det()
            if d <= 0: ok = False
    # orthogonality between different beta sectors is automatic (charge conservation)
    print('(ii) weight', wt, ': H_0 positive definite on all', len(bybeta), 'lattice sectors:', ok, '(%.0fs)' % (time.time() - t0), flush=True)
    res['positive_weight_%d' % wt] = ok
json.dump(res, open('phi_automorphism_check.json', 'w'), indent=1)
