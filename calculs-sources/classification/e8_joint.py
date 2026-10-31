"""Explicit e8 (Frenkel--Kac) containing e6 + sl3, and joint traces of commuting elements (Axis C3/C2).

For the two candidate weight-one algebras s = g' + g'' (g'' diagonal in e6 + sl3 with p2(g'') the principal so3,
g' the centralizer in e6 of p1(g'')), we compute exactly the joint power sums
   Tr_R (H')^p (H'')^q   on R = e6 (adjoint), 27 (+27bar), (27,3), (1,3), sl3 (adjoint),
from traces of ad on e8 = (78,1) + (1,8) + (27,3) + (27bar,3bar):  for even total degree,
   Tr_{(27,3)} = (Tr_{e8} - Tr_{e6} - Tr_{sl3})/2 ,  and for y in e6:  Tr_27 = (Tr_{e8} - Tr_{e6})/6 .
H' is a hyperbolic element x of the centralizer c, rescaled: x ~ t H' with (H',H') = 2 i1'; the rescaling of
power sums is rational for even powers of x.
"""
import random, itertools, sys
from fractions import Fraction as Fr
import sympy as sp

CE8 = [[2,0,-1,0,0,0,0,0],[0,2,0,-1,0,0,0,0],[-1,0,2,-1,0,0,0,0],[0,-1,-1,2,-1,0,0,0],
       [0,0,0,-1,2,-1,0,0],[0,0,0,0,-1,2,-1,0],[0,0,0,0,0,-1,2,-1],[0,0,0,0,0,0,-1,2]]
n = 8
def ip(a, b):
    return sum(a[i]*CE8[i][j]*b[j] for i in range(n) for j in range(n))
def root_system():
    roots = set(); fr = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
    roots.update(fr)
    while fr:
        new = []
        for r in fr:
            for i in range(n):
                c = ip(r, [1 if t == i else 0 for t in range(n)])
                s = tuple(r[k] - (c if k == i else 0) for k in range(n))
                if s not in roots:
                    roots.add(s); new.append(s)
        fr = new
    return sorted(roots)
roots = root_system(); assert len(roots) == 240
ridx = {r: i for i, r in enumerate(roots)}
DIM = n + 240
def eps(a, b):
    e = sum(a[i]*b[j]*(CE8[i][j] if i < j else (1 if i == j else 0)) for i in range(n) for j in range(n))
    return -1 if e % 2 else 1
# sparse structure constants
def bracket_pairs():
    S = {}
    for x in range(DIM):
        for y in range(DIM):
            if x < n and y < n:
                continue
            if x < n:
                a = roots[y - n]; c = ip([1 if k == x else 0 for k in range(n)], a)
                if c: S[(x, y)] = {y: c}
            elif y < n:
                a = roots[x - n]; c = ip([1 if k == y else 0 for k in range(n)], a)
                if c: S[(x, y)] = {x: -c}
            else:
                a, b = roots[x - n], roots[y - n]; s = tuple(p + q for p, q in zip(a, b))
                if all(v == 0 for v in s):
                    S[(x, y)] = {i: eps(a, b)*a[i] for i in range(n) if a[i]}
                elif s in ridx:
                    S[(x, y)] = {n + ridx[s]: eps(a, b)}
    return S
STR = bracket_pairs()
def br(u, v):
    out = {}
    for x, ux in u.items():
        for y, vy in v.items():
            for k, c in STR.get((x, y), {}).items():
                out[k] = out.get(k, 0) + ux*vy*c
    return {k: v for k, v in out.items() if v != 0}
def admat(u):
    """sparse ad(u) as dict of columns"""
    return {y: br(u, {y: 1}) for y in range(DIM)}
def matpow_trace(Ms, word):
    """trace of the product of the sparse matrices Ms[w] for w in word (applied right to left)"""
    tot = 0
    for y in range(DIM):
        v = {y: 1}
        for w in reversed(word):
            nv = {}
            for k, c in v.items():
                for kk, cc in Ms[w][k].items():
                    nv[kk] = nv.get(kk, 0) + c*cc
            v = {k: c for k, c in nv.items() if c != 0}
            if not v: break
        tot += v.get(y, 0)
    return tot
# subalgebras: e6 = roots supported on nodes 1..6 ; sl3 = roots orthogonal to e6
E6_roots = [r for r in roots if all(r[i] == 0 for i in (6, 7))]
assert len(E6_roots) == 72
A2_roots = [r for r in roots if all(ip(r, [1 if t == i else 0 for t in range(n)]) == 0 for i in range(6))]
assert len(A2_roots) == 6
E6_basis = list(range(6)) + [n + ridx[r] for r in E6_roots]
# Cartan of sl3: the coroots of A2_roots (span)
def cartan_vec(c):          # element sum c_i alpha_i^vee of the e8 Cartan
    return {i: Fr(c[i]) for i in range(n) if c[i] != 0}

CE6 = [row[:6] for row in CE8[:6]]
def lin_solve(rows, rhs, nvar):
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    sol, params = M.gauss_jordan_solve(b)
    sol = sol.subs({p: 0 for p in params})
    return [Fr(int(sp.fraction(t)[0]), int(sp.fraction(t)[1])) for t in sol]
def to_sp(x):
    return sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else sp.Integer(x)

def e6_triple(labels, seed):
    c = list(sp.Matrix([labels])*sp.Matrix(CE6).inv())
    H = {i: Fr(int(sp.fraction(c[i])[0]), int(sp.fraction(c[i])[1])) for i in range(6) if c[i] != 0}
    lev = lambda r: sum(r[i]*labels[i] for i in range(6))
    g2 = [n + ridx[r] for r in E6_roots if lev(r) == 2]
    gm2 = [n + ridx[r] for r in E6_roots if lev(r) == -2]
    random.seed(seed)
    e = {k: Fr(random.randint(1, 9)) for k in g2}
    cols = [br(e, {k: 1}) for k in gm2]
    keys = sorted(set().union(*[set(c_) for c_ in cols]) | set(H))
    rows = [[to_sp(cols[t].get(key, 0)) for t in range(len(gm2))] for key in keys]
    rhs = [to_sp(H.get(key, 0)) for key in keys]
    try:
        y = lin_solve(rows, rhs, len(gm2))
    except ValueError:
        return None
    f = {gm2[t]: y[t] for t in range(len(gm2)) if y[t] != 0}
    if br(e, f) != {k: v for k, v in H.items()}:
        return None
    return H, e, f

def centralizer_e6(H, e, f):
    B = E6_basis
    cols = {}
    rows = []
    for u in (H, e, f):
        imgs = [br(u, {b: 1}) for b in B]
        keys = sorted(set().union(*[set(i) for i in imgs]))
        rows += [[to_sp(imgs[t].get(key, 0)) for t in range(len(B))] for key in keys]
    ns = sp.Matrix(rows).nullspace()
    out = []
    for v in ns:
        out.append({B[t]: Fr(int(sp.fraction(v[t])[0]), int(sp.fraction(v[t])[1])) for t in range(len(B)) if v[t] != 0})
    return out

def tr(sub, word, basis=None):
    Ms = {w: admat(w_) for w, w_ in sub.items()}
    if basis is None:
        return matpow_trace(Ms, word)
    tot = 0
    for y in basis:
        v = {y: 1}
        for w in reversed(word):
            nv = {}
            for k, cc in v.items():
                for kk, c2 in Ms[w][k].items():
                    nv[kk] = nv.get(kk, 0) + cc*c2
            v = {k: c_ for k, c_ in nv.items() if c_ != 0}
            if not v: break
        tot += v.get(y, 0)
    return tot

def joint_sums(labels, i1p, degs=((2, 0), (0, 2), (4, 0), (0, 4), (2, 2), (6, 0), (0, 6), (4, 2), (2, 4)), seed0=0):
    """joint power sums Tr (H')^p (H'')^q on the sectors e6, 27, (27,3), (1,3), sl3 and on s, for
    H'' = (H''_e6 with the given labels) + (principal so3 neutral element of sl3), H' = neutral element of
    the centralizer (in e6) of the triple of H''_e6, normalized by (H', H') = 2 i1p (exact rescaling)."""
    for seed in range(seed0, seed0 + 60):
        T = e6_triple(labels, seed)
        if T: break
    Hd, e, f = T
    C = centralizer_e6(Hd, e, f)
    assert len(C) == 3
    random.seed(7)
    while True:
        coef = [random.randint(-3, 3) for _ in C]
        x = {}
        for c_, v in zip(coef, C):
            for k, val in v.items():
                x[k] = x.get(k, 0) + c_*val
        x = {k: v for k, v in x.items() if v != 0}
        kappa = Fr(tr({'x': x}, ['x', 'x'], E6_basis), 24)
        if kappa > 0:
            break
    # sl3 principal: 2 beta1 + 2 beta2 with beta1 = alpha8, beta2 = theta - alpha8 (coroot coordinates)
    beta1 = [0]*7 + [1]; beta2 = [2, 3, 4, 6, 5, 4, 3, 1]
    Hs = {i: Fr(2*beta1[i] + 2*beta2[i]) for i in range(n) if 2*beta1[i] + 2*beta2[i] != 0}
    Hsum = dict(Hd)
    for k, v in Hs.items():
        Hsum[k] = Hsum.get(k, 0) + v
    out = {}
    for p, q in degs:
        sc = (Fr(2*i1p)/kappa)**(p//2)
        word = ['x']*p + ['H']*q
        e8tot = tr({'x': x, 'H': Hsum}, word) if q or p else 248
        e6tot = tr({'x': x, 'H': Hd}, word, E6_basis)
        # sl3 adjoint under the principal so3: eigenvalues 2,2,-2,-2,4,-4,0,0 ; x acts trivially
        sl3 = 0 if p else sum(ev**q for ev in (2, 2, -2, -2, 4, -4, 0, 0))
        R27_3 = Fr(e8tot - e6tot - sl3, 2)
        # 27 alone: use H''_e6 only
        e8_e6only = tr({'x': x, 'H': Hd}, word)
        R27 = Fr(e8_e6only - e6tot, 6)
        R13 = 0 if p else sum(ev**q for ev in (2, 0, -2))
        sdat = (2*2**p if q == 0 else 0) + (2*2**q if p == 0 else 0)
        out[(p, q)] = {'e6': e6tot*sc, 'sl3': sl3, '27': R27*sc, '27x3': R27_3*sc, '1x3': R13, 's': sdat}
    return out

def joint_weights(labels, i1p, seed0=0):
    """joint weights (H'-eigenvalue, H''_e6-eigenvalue) on e6 (adjoint) and on 27, exact.
    H' = t^{-1} x with x hyperbolic in the centralizer, t^2 = kappa/(2 i1p)."""
    for seed in range(seed0, seed0 + 60):
        T = e6_triple(labels, seed)
        if T: break
    Hd, e, f = T
    C = centralizer_e6(Hd, e, f)
    random.seed(7)
    while True:
        coef = [random.randint(-3, 3) for _ in C]
        x = {}
        for c_, v in zip(coef, C):
            for k, val in v.items():
                x[k] = x.get(k, 0) + c_*val
        x = {k: v for k, v in x.items() if v != 0}
        kappa = Fr(tr({'x': x}, ['x', 'x'], E6_basis), 24)
        if kappa > 0:
            break
    r = sp.Rational(kappa.numerator, kappa.denominator)/(2*i1p)        # t^2
    # H''_e6 eigenvalue of each basis vector of e8
    def hev(b):
        if b < n:
            return 0
        a = roots[b - n]
        return sum(Hd.get(i, 0)*ip([1 if t_ == i else 0 for t_ in range(n)], a) for i in range(n))
    Mx = admat(x)
    X = sp.Symbol('X')
    def weights_on(basis):
        out = []
        groups = {}
        for b in basis:
            groups.setdefault(hev(b), []).append(b)
        for lam, bs in groups.items():
            idx = {b: i for i, b in enumerate(bs)}
            A_ = sp.zeros(len(bs), len(bs))
            for b in bs:
                for k, v in Mx[b].items():
                    assert k in idx, 'x does not preserve the eigenspace'
                    A_[idx[k], idx[b]] = sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else v
            P = sp.Poly(A_.charpoly(X).as_expr(), X)
            rem = P
            mults = {}
            for m in range(30, 0, -1):
                fac = sp.Poly(X**2 - r*m*m, X)
                while True:
                    qq, rr = sp.div(rem, fac)
                    if rr.is_zero:
                        rem = qq; mults[m] = mults.get(m, 0) + 1
                    else:
                        break
            z0 = rem.degree()
            assert sp.expand(rem.monic().as_expr() - X**z0) == 0, ('non-integral eigenvalues', rem, r)
            out += [(0, int(lam))]*z0
            for m, c_ in mults.items():
                out += [(m, int(lam))]*c_ + [(-m, int(lam))]*c_
        return out
    ad6 = weights_on(E6_basis)
    rest = [b for b in range(DIM) if b not in set(E6_basis)]
    w_rest = weights_on(rest)                 # = 3*27 + 3*27bar + 8*(0,0)  (sl3 part trivial for x and H''_e6)
    from collections import Counter
    cnt = Counter(w_rest); cnt[(0, 0)] -= 8
    assert all(v % 6 == 0 for v in cnt.values()), 'unexpected multiplicities'
    w27 = []
    for k_, v in cnt.items():
        w27 += [k_]*(v//6)
    assert len(w27) == 27
    return ad6, w27
