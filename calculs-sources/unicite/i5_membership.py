"""Axis I, final step: are the images of Jiao's seeds u26, u17 (under the explicit automorphism psi, torus) in one of
our four subspaces U_2?  Exact membership over K (exact_engine.solve_block).  Scalars: u17 maps to sqrt(-5) times a
K-vector (irrelevant for membership).  Subalgebra 10 = sigma(subalgebra 9), sigma: sqrt13 -> -sqrt13."""
import sys, os, json, pickle
HERE = os.path.abspath('.'); sys.path.insert(0, HERE); sys.path.insert(0, os.path.abspath('../generation-forte')); os.chdir('../generation-forte')
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, exact_engine as E, Kfield as K
R = json.load(open(HERE + '/i4_torus.json'))
s13, s3 = sp.sqrt(13), sp.sqrt(3)
def to_K(expr):
    x = sp.expand(sp.sympify(expr).subs(sp.sqrt(39), s13 * s3).subs(sp.sqrt(195), s13 * s3 * sp.sqrt(5))
                  .subs(sp.sqrt(65), s13 * sp.sqrt(5)).subs(sp.sqrt(15), s3 * sp.sqrt(5)))
    return x
A13, A3, A5 = sp.symbols('A13 A3 A5', positive=True)
def radsub(e):
    e = sp.sympify(e)
    for n, r in ((195, A13*A3*A5), (65, A13*A5), (39, A13*A3), (15, A3*A5), (13, A13), (5, A5), (3, A3)):
        e = e.subs(sp.sqrt(n), r)
    return sp.expand(e)
import mpmath
mpmath.mp.dps = 60
JJ = sp.Symbol('JJ')
def kvec(d, scale):
    out = {}
    for ks, c in d.items():
        x = radsub(c).subs(sp.I, JJ) / radsub(scale).subs(sp.I, JJ)
        fn = sp.lambdify((A13, A3, A5, JJ), x, 'mpmath')
        r13, r3, r5, i = mpmath.sqrt(13), mpmath.sqrt(3), mpmath.sqrt(5), mpmath.mpc(0, 1)
        xs = [fn(s1 * r13, r3, r5, s2 * i) for s1, s2 in ((1, 1), (-1, 1), (1, -1), (-1, -1))]
        comps = [(xs[0] + xs[1] + xs[2] + xs[3]) / 4, (xs[0] - xs[1] + xs[2] - xs[3]) / (4 * r13),
                 (xs[0] + xs[1] - xs[2] - xs[3]) / (4 * i * r3), (xs[0] - xs[1] - xs[2] + xs[3]) / (4 * r13 * i * r3)]
        vals = []
        for z in comps:
            assert abs(mpmath.im(z)) < mpmath.mpf(10) ** -40, z
            q = sp.nsimplify(sp.Float(mpmath.nstr(mpmath.re(z), 50), 50), rational=True, tolerance=1e-40)
            q = sp.Rational(q).limit_denominator(10 ** 6)
            assert abs(sp.Float(q, 50) - sp.Float(mpmath.nstr(mpmath.re(z), 50), 50)) < 1e-40
            vals.append(q)
        key = (tuple(int(t) for t in ks.strip('()').split(',')), ())
        out[key] = tuple(mpq(int(v.p), int(v.q)) for v in vals)
    return out
u26 = kvec(R['U26 coefficients'], 1)
u17 = kvec(R['U17 coefficients'], sp.sqrt(5) * sp.I)
def sigma(v): return {k: (a[0], -a[1], a[2], -a[3]) for k, a in v.items()}
res = {}
for chi in (0, 1):
    S = E.setup(chi); T = S['T']; B2 = E.Blocks(T, 2)
    piv = {}
    for v in S['U2span']:
        mu = B2.wt(v); piv.setdefault(mu, []).append(v)
    # reduce to independent pivots per block
    basis = {}
    for mu, vs in piv.items():
        idx = E.select_pivots(vs, B2.pos[mu], len(B2.blk[mu])); basis[mu] = [vs[i] for i in idx]
    for name, vec in (('u26', u26), ('u17', u17)):
        for sg in (0, 1):
            w = sigma(vec) if sg else vec
            mu = B2.wt(w)
            ok, bad = E.solve_block(basis.get(mu, []), [w], B2.pos[mu], len(B2.blk[mu]))
            res['%s in U_2(subalgebra %d, choice %d)' % (name, 10 if sg else 9, chi)] = bool(ok and not bad)
print(json.dumps(res, indent=1))
json.dump(res, open(HERE + '/i5_membership.json', 'w'), indent=1)
