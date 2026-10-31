"""Axis E: search for a positive lattice character chi such that the antilinear involution
   phi_chi( h_1(-n_1)...h_l(-n_l) e^lam ) = (-1)^l chi(lam) conj(coef) h_1(-n_1)...h_l(-n_l) e^{-lam}
(Riesen, Proposition 3.1 / FLM: the standard unitary structure transported by the torus automorphism
g(e^lam) = chi(lam)^{-1/2} e^lam) maps e' -> -kappa' f', e'' -> -kappa'' f'' (kappa > 0).  Then phi_chi(U_1) = U_1.
Linear system in log chi(alpha_i), log kappa."""
import json, math
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, cases as C
T = C.build('II')
eqs = []
for e, f, name in ((T['e1'], T['f1'], 'kap1'), (T['e2'], T['f2'], 'kap2')):
    supp_e = {k[0]: c for k, c in e.items()}
    supp_f = {k[0]: c for k, c in f.items()}
    for a, c in supp_e.items():
        na = V.neg(a)
        if na not in supp_f:
            print('support mismatch', a); continue
        ratio = -supp_f[na]/c
        eqs.append((a, ratio, name))
print('equations:', len(eqs), '; all ratios positive:', all(r > 0 for _, r, _ in eqs))
# unknowns: x_i = log chi(alpha_i) (8), y1 = log kappa1, y2 = log kappa2 ; chi(a) = kappa * ratio
xs = sp.symbols('x0:8'); y1, y2 = sp.symbols('y1 y2')
rows = []
for a, r, name in eqs:
    if r <= 0: continue
    lhs = sum(a[i]*xs[i] for i in range(8)) - (y1 if name == 'kap1' else y2)
    rows.append(sp.Eq(lhs, sp.log(sp.Rational(int(r.numerator), int(r.denominator)))))
sol = sp.solve(rows, list(xs) + [y1, y2], dict=True)
print('solution exists:', bool(sol))
if sol:
    s = sol[0]
    print({str(k): sp.nsimplify(sp.exp(v)) for k, v in s.items()})
    json.dump({str(k): str(sp.simplify(sp.exp(v))) for k, v in s.items()}, open('unitary_character.json', 'w'), indent=1)
