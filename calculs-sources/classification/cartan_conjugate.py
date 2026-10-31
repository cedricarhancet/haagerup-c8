"""Put the two commuting neutral elements H' and H''_e6 simultaneously into the standard Cartan t of e6.
We search H'_t in t with integral values alpha_i(H'_t) in [-3, 3] such that the joint multiset of
(H'_t, H''_e6)-weights on the 27 and on the roots equals the exactly computed joint weights (e8_joint).
Any such H'_t commutes with H''_e6 and generates, together with H''_e6, a torus conjugate to the one
spanned by (H', H''_e6) (equality of the joint characters of a faithful representation)."""
import itertools, json
from collections import Counter
from fractions import Fraction as Fr
import sympy as sp
import e6_tools as E, e8_joint as J
def find(lab, i1p):
    ad6, w27 = J.joint_weights(lab, i1p)
    target27 = Counter(w27)
    targetad = Counter([w for w in ad6 if w != (0, 0)] )
    Mi = sp.Matrix(E.C).inv()
    cH = [Fr(int(sp.fraction(t)[0]), int(sp.fraction(t)[1])) for t in sp.Matrix([lab])*Mi]
    h2_27 = [E.ev(cH, w) for w in E.W27]
    h2_r = [E.ev(cH, r) for r in E.ROOTS]
    sols = []
    for l in itertools.product(range(-3, 4), repeat=6):
        c = [Fr(x) for x in (sp.Matrix([l])*Mi)]
        w = Counter((E.ev(c, x), h) for x, h in zip(E.W27, h2_27))
        if w != target27:
            continue
        r = Counter((E.ev(c, x), h) for x, h in zip(E.ROOTS, h2_r))
        r[(0, 0)] += 0
        rr = Counter({k: v for k, v in r.items()})
        tad = Counter(ad6); tad[(0, 0)] -= 6
        if rr != Counter({k: v for k, v in tad.items() if v}):
            continue
        sols.append(l)
        if len(sols) >= 3:
            break
    return sols
if __name__ == '__main__':
    out = {}
    for name, lab, i1p in (('A11+A139', (2, 1, 1, 0, 1, 2), 1), ('A13+A113', (1, 0, 0, 1, 0, 1), 3)):
        s = find(lab, i1p)
        print(name, 'H\'_t labels (alpha_i(H\'_t)):', s[:3])
        out[name] = [list(x) for x in s]
    json.dump(out, open('cartan_conjugate.json', 'w'), indent=1)
