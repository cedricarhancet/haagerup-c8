"""Exhaustive check (exact) that the standard antilinear involution phi_0 of W = V_{E8+E8} (our cocycle) is positive in
weights <= 2: for every lattice vector beta of norm <= 4, the Gram matrix of (phi_0 x, y) on the monomial basis of the
sector e^beta in weight <= 2 is positive definite (leading principal minors), and phi_0 maps sectors beta -> -beta.
Since phi_chi = g phi_0 g^{-1} (g = chi^{-1/2} positive torus element), phi_chi is positive as well."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time
from gmpy2 import mpq
from flint import fmpq_mat
import e8e8_voa as V
t0 = time.time()
def phi0(v): return {(V.NEGI[b], m): c*(-1)**len(m) for (b, m), c in v.items()}
bad = 0; nsec = 0
for wt in (0, 1, 2):
    keys = V.basis(wt); by = {}
    for k in keys: by.setdefault(k[0], []).append(k)
    for b, ks in by.items():
        if len(ks) == 1:
            f = V.form(phi0({ks[0]: mpq(1)}), {ks[0]: mpq(1)})
            nsec += 1
            if f <= 0: bad += 1
            continue
        n = len(ks)
        M = fmpq_mat(n, n)
        for i in range(n):
            for j in range(n):
                f = V.form(phi0({ks[i]: mpq(1)}), {ks[j]: mpq(1)})
                M[i, j] = fmpq_mat([[f]])[0, 0] if False else int(f.numerator) if f.denominator == 1 else 0
                if f.denominator != 1: M[i, j] = fmpq_mat([[0]])[0, 0] + __import__('flint').fmpq(int(f.numerator), int(f.denominator))
        for m in range(1, n + 1):
            if fmpq_mat([[M[i, j] for j in range(m)] for i in range(m)]).det() <= 0: bad += 1; break
        nsec += 1
    print('weight', wt, ': sectors checked', len(by), '; failures so far', bad, '(%.0fs)' % (time.time() - t0), flush=True)
print('phi_0 positive definite on W_0 + W_1 + W_2 (all %d lattice sectors):' % nsec, bad == 0)
