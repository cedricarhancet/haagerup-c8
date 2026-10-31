"""Explicit e6 (Frenkel--Kac construction on the E6 root lattice) and centralizers of sl2-triples (Axis C3).

g = h + sum_alpha C e_alpha, [h, e_a] = a(h) e_a, [e_a, e_b] = eps(a,b) e_{a+b} if a+b is a root,
[e_a, e_{-a}] = eps(a,-a) a (identifying h with h* via the normalized form), eps bimultiplicative with
eps(a,b) eps(b,a) = (-1)^{(a,b)}.  For a weighted Dynkin diagram, we take H with alpha_i(H) = label_i, a
random integral e in g_2, solve [e,f] = H (f in g_{-2}) exactly, check the sl2 relations, and compute the
centralizer c of the triple and the Dynkin index of c in e6 (when c is simple of rank one):
   index(c) = (Tr_{e6} (ad x)^2 / 24) / (Tr_c (ad x|_c)^2 / 4)  for x in c with Tr_c (ad x)^2 != 0.
"""
import random, sys, itertools
from fractions import Fraction as Fr
import sympy as sp
C = [[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]
n = 6
def ip(a, b):
    return sum(a[i]*C[i][j]*b[j] for i in range(n) for j in range(n))
# roots in simple-root coordinates
roots = set()
frontier = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
roots.update(frontier)
while frontier:
    new = []
    for r in frontier:
        for i in range(n):
            s = tuple(r[k] - (ip(r, [1 if t == i else 0 for t in range(n)]))*(1 if k == i else 0) for k in range(n))
            if s not in roots:
                roots.add(s); new.append(s)
    frontier = new
roots = sorted(roots)
assert len(roots) == 72 and all(ip(r, r) == 2 for r in roots)
ridx = {r: i for i, r in enumerate(roots)}
DIM = n + 72
def eps(a, b):
    e = 0
    for i in range(n):
        for j in range(n):
            if i < j:
                e += a[i]*b[j]*C[i][j]
            elif i == j:
                e += a[i]*b[j]*1
    return -1 if e % 2 else 1
# structure constants: basis 0..5 = coroots alpha_i^vee (h_i), 6.. = e_alpha
def bracket_basis(x, y):
    """returns dict index->coeff for [b_x, b_y]"""
    if x < n and y < n:
        return {}
    if x < n:                       # [h_i, e_a] = a(h_i) e_a = (alpha_i, a) e_a
        a = roots[y - n]
        c = ip([1 if k == x else 0 for k in range(n)], a)
        return {y: Fr(c)} if c else {}
    if y < n:
        r = bracket_basis(y, x)
        return {k: -v for k, v in r.items()}
    a, b = roots[x - n], roots[y - n]
    s = tuple(p + q for p, q in zip(a, b))
    if all(v == 0 for v in s):      # [e_a, e_{-a}] = eps(a,-a) a, a = sum a_i alpha_i -> coroot coords
        return {i: Fr(eps(a, b)*a[i]) for i in range(n) if a[i]}
    if s in ridx:
        return {n + ridx[s]: Fr(eps(a, b))}
    return {}
STR = [[bracket_basis(x, y) for y in range(DIM)] for x in range(DIM)]
def br(u, v):
    out = [Fr(0)]*DIM
    for x, ux in enumerate(u):
        if ux == 0: continue
        for y, vy in enumerate(v):
            if vy == 0: continue
            for k, c in STR[x][y].items():
                out[k] += ux*vy*c
    return out
def ad(u):
    return [[sum(u[x]*STR[x][y].get(k, 0) for x in range(DIM) if u[x]) for y in range(DIM)] for k in range(DIM)]
def check_jacobi(trials=30):
    random.seed(1)
    for _ in range(trials):
        u, v, w = [[Fr(random.randint(-2, 2)) for _ in range(DIM)] for _ in range(3)]
        j = [p + q + r for p, q, r in zip(br(u, br(v, w)), br(v, br(w, u)), br(w, br(u, v)))]
        assert all(t == 0 for t in j)
    return True
def trace_ad2(x):
    A = sp.Matrix(ad(x))
    return (A*A).trace()

def triple(labels, seed=0):
    Mi = sp.Matrix(C).inv()
    c = list(sp.Matrix([labels])*Mi)
    H = [Fr(int(sp.fraction(t)[0]), int(sp.fraction(t)[1])) for t in c] + [Fr(0)]*72
    g2 = [k for k, r in enumerate(roots) if sum(r[i]*labels[i] for i in range(n)) == 2]
    gm2 = [k for k, r in enumerate(roots) if sum(r[i]*labels[i] for i in range(n)) == -2]
    random.seed(seed)
    e = [Fr(0)]*DIM
    for k in g2:
        e[n + k] = Fr(random.randint(1, 7))
    # solve [e, f] = H with f = sum y_k e_{-..}
    ys = sp.symbols('y0:%d' % len(gm2))
    eqs = [sp.Integer(0)]*DIM
    for t, k in enumerate(gm2):
        b = [Fr(0)]*DIM; b[n + k] = Fr(1)
        r = br(e, b)
        for i in range(DIM):
            if r[i]:
                eqs[i] += sp.Rational(r[i].numerator, r[i].denominator)*ys[t]
    eqs = [eqs[i] - sp.Rational(H[i].numerator, H[i].denominator) for i in range(DIM)]
    sol = sp.solve([q for q in eqs if q != 0], ys, dict=True)
    if not sol:
        return None
    sol = sol[0]
    f = [Fr(0)]*DIM
    for t, k in enumerate(gm2):
        val = sol.get(ys[t], 0)
        val = sp.Rational(val)
        f[n + k] = Fr(int(val.p), int(val.q))
    assert br(H, e) == [2*x for x in e] and br(H, f) == [-2*x for x in f] and br(e, f) == H
    return H, e, f

def centralizer(H, e, f):
    rows = []
    for u in (H, e, f):
        A = ad(u)
        rows += [[A[k][y] for y in range(DIM)] for k in range(DIM)]
    Msp = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in r] for r in rows])
    ns = Msp.nullspace()
    return [[Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in v] for v in ns]

if __name__ == '__main__':
    assert check_jacobi()
    print('e6 built: dim', DIM, '; Jacobi identity checked on random triples')
    cases = {'index 9 (2A2+A1)': (1, 0, 0, 1, 0, 1), 'index 35 (A5)': (2, 1, 1, 0, 1, 2), 'index 3 (3A1)': (0, 0, 0, 1, 0, 0),
             'index 1 (A1)': (0, 1, 0, 0, 0, 0)}
    for name, lab in cases.items():
        t = None
        for seed in range(50):        # a random integral e in g_2 is generic for most seeds
            t = triple(lab, seed)
            if t is not None:
                break
        if t is None:
            print(name, ': no triple found'); continue
        H, e, f = t
        idx = Fr(ip([0]*6, [0]*6))
        hc = [x for x in H[:n]]
        index = sp.Rational(trace_ad2(H), 24)/2
        cz = centralizer(H, e, f)
        print('%s: Dynkin index of the triple = %s, centralizer dim = %d' % (name, index, len(cz)))
        if len(cz) == 3:
            # structure of c: compute brackets in the basis cz
            B = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in v] for v in cz]).T
            def coords(v):
                vv = sp.Matrix([sp.Rational(x.numerator, x.denominator) for x in v])
                sol = B.solve_least_squares(vv) if False else B.pinv()*vv
                return sol
            # Killing form of c and index
            def adc(x):
                return sp.Matrix.hstack(*[coords(br(x, y)) for y in cz])
            for x in cz + [[a + b for a, b in zip(cz[0], cz[1])], [a + 2*b + 3*c for a, b, c in zip(cz[0], cz[1], cz[2])]]:
                tc = (adc(x)**2).trace()
                if tc != 0:
                    print('    centralizer: [c,c] = c ?', sp.Matrix.hstack(*[coords(br(cz[i], cz[j])) for i in range(3) for j in range(i + 1, 3)]).rank() == 3,
                          '; index of c in e6 =', sp.nsimplify(sp.Rational(trace_ad2(x), 24)/(tc/4)))
                    break
