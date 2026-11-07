"""phi_chi on W (exact): phi_chi(c mono e^beta) = (-1)^{#mono} chi(beta) conj(c) mono e^{-beta},
chi(beta) = 2^{-(c_8(beta^(1)) + c_8(beta^(2)))} (see r3b_chi.py / r3b_chi.json)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gmpy2 import mpq
import e8e8_voa as V
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
    return (phi_rat(v[0]), {k: -c for k, c in phi_rat(v[1]).items()})
