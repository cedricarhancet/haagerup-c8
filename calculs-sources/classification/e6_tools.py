"""Exact E6 data: Cartan matrix, roots, weights of the 27 (fundamental-weight coordinates), quartic traces.
An element h of the Cartan subalgebra is written h = sum_i c_i alpha_i^vee (simple coroots); the eigenvalue
of h on a weight w = sum w_i omega_i is sum_i c_i w_i, and (h, h) = c^T C c (normalized form, roots of norm 2).
"""
from fractions import Fraction as Fr
import itertools
C = [[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]
def orbit(w0):
    seen = {tuple(w0)}; st = [tuple(w0)]
    while st:
        w = st.pop()
        for i in range(6):
            if w[i] != 0:
                v = tuple(w[k] - w[i]*C[i][k] for k in range(6))
                if v not in seen:
                    seen.add(v); st.append(v)
    return sorted(seen)
W27 = orbit((1,0,0,0,0,0)); assert len(W27) == 27
# roots in fundamental-weight coordinates: orbit of the highest root (= omega_2 in Bourbaki labelling)
ROOTS = orbit((0,1,0,0,0,0)); assert len(ROOTS) == 72
def ev(c, w):
    return sum(Fr(ci)*wi for ci, wi in zip(c, w))
def form(c, c2=None):
    c2 = c if c2 is None else c2
    return sum(Fr(c[i])*C[i][j]*Fr(c2[j]) for i in range(6) for j in range(6))
def tr27(c, k):
    return sum(ev(c, w)**k for w in W27)
def trad(c, k):
    return sum(ev(c, r)**k for r in ROOTS)
def kappas():
    out = set()
    for c in [(1,0,0,0,0,0), (0,0,0,1,0,0), (1,2,-1,0,3,1), (2,-1,1,1,0,-3)]:
        hh = form(c)
        assert tr27(c, 2) == 6*hh and trad(c, 2) == 24*hh
        assert tr27(c, 3) == 0 and trad(c, 3) == 0
        out.add((tr27(c, 4)/hh**2, trad(c, 4)/hh**2))
    assert len(out) == 1
    return out.pop()
if __name__ == '__main__':
    k27, kad = kappas()
    print('Tr_27 h^4 = %s (h,h)^2 ;  Tr_e6 (ad h)^4 = %s (h,h)^2' % (k27, kad))
