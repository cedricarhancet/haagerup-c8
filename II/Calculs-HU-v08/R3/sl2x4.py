"""Explicit sl2^4 = (S x S)_1 inside e8 + e8 = W_1 (R3a step 1), exact over Q.

Model (T12, t12_e8split.py): in each copy e8^(k) = so(16) + 128, the level-3 factor g'_k is the 2A1 subalgebra so(3)
acting on C^3 = <f_0^+, f_0^-, u>, u = f_1^+ + f_1^- ; its centraliser so(13) acts on C^13 = <u'> + <f_2..f_7>,
u' = f_1^+ - f_1^-, and the three other factors act on C^13 through U' = C^4 + C^9:
   C^4 = <f_2^pm, f_3^pm> = V1 x V1 (roots e2+e3 and e2-e3 of so(4)),  C^9 = <u', f_4..f_7> = V2 x V2.
copy 1 (indices (2,6,1,7)): g'_2 = root e2+e3 ; g''_2 = root e2-e3 + second factor of V2xV2 ; g''_1 = first factor.
copy 2 (indices (1,7,2,6)): g'_1 = root e2+e3 ; g''_1 = root e2-e3 + second factor ; g''_2 = first factor.
C^16 has the hyperbolic basis f_k^+ (index 2k), f_k^- (index 2k+1), B(f_k^+, f_k^-) = 1.  so(16) is identified with
the D8 part of e8 by T_k = f_k^+ ^ f_k^- -> eps_k(-1)1 and X(s e_k + t e_l) = f_k^s ^ f_l^t -> sigma e^{s e_k + t e_l},
with signs sigma determined by BFS and then checked on ALL brackets.  Output: sl2x4.pkl (Chevalley triples
(E, H, F) of g'_1, g''_1, g'_2, g''_2 as W_1-vectors, exact)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, itertools
from fractions import Fraction as Fr
from gmpy2 import mpq
import e8e8_voa as V

N16 = 16
def zero(): return [[Fr(0)]*N16 for _ in range(N16)]
Bm = zero()
for k in range(8): Bm[2*k][2*k+1] = Bm[2*k+1][2*k] = Fr(1)
def vec(d):
    v = [Fr(0)]*N16
    for i, c in d.items(): v[i] = Fr(c)
    return v
def Bf(a, b): return sum(a[i]*Bm[i][j]*b[j] for i in range(N16) for j in range(N16) if a[i] and b[j])
def wedge(a, b):
    """(a ^ b) v = B(b, v) a - B(a, v) b"""
    Ba = [sum(a[i]*Bm[i][j] for i in range(N16)) for j in range(N16)]
    Bb = [sum(b[i]*Bm[i][j] for i in range(N16)) for j in range(N16)]
    return [[a[r]*Bb[c] - b[r]*Ba[c] for c in range(N16)] for r in range(N16)]
def mmul(A, B_): return [[sum(A[i][k]*B_[k][j] for k in range(N16) if A[i][k]) for j in range(N16)] for i in range(N16)]
def madd(A, B_, c=1): return [[A[i][j] + c*B_[i][j] for j in range(N16)] for i in range(N16)]
def br(A, B_): return madd(mmul(A, B_), mmul(B_, A), -1)
def fp(k): return vec({2*k: 1})
def fm(k): return vec({2*k+1: 1})
def f(k, s): return fp(k) if s > 0 else fm(k)
# so(16) basis
ROOTS16 = []
for k in range(8):
    for l in range(k+1, 8):
        for s in (1, -1):
            for t in (1, -1):
                ROOTS16.append((k, s, l, t))
def Xm(r):
    k, s, l, t = r
    return wedge(f(k, s), f(l, t))
def Tm(k): return wedge(fp(k), fm(k))
def decompose(A):
    """A in so(16) -> ({k: d_k}, {root: c})"""
    d = {k: A[2*k][2*k] for k in range(8) if A[2*k][2*k]}
    c = {}
    for r in ROOTS16:
        k, s, l, t = r
        row = 2*k if s > 0 else 2*k+1
        col = 2*l + (1 if t > 0 else 0)          # f_l^{-t}
        if A[row][col]: c[r] = A[row][col]
    # check
    Bk = zero()
    for kk, x in d.items(): Bk = madd(Bk, Tm(kk), x)
    for r, x in c.items(): Bk = madd(Bk, Xm(r), x)
    assert Bk == A, 'decomposition failed'
    return d, c
def rootvec(r):             # doubled coordinates in one E8 copy
    k, s, l, t = r
    y = [0]*8; y[k] = 2*s; y[l] = 2*t
    return tuple(y)

# ------------------------------------------------- signs sigma (copy-independent: cocycle is a product)
def lat_root_idx(r, copy):
    y = rootvec(r)
    full = y + (0,)*8 if copy == 0 else (0,)*8 + y
    return V.IDX[full]
def lattice_bracket(x, y): return V.clean(V.mode_vec(x, 0, y))
def find_sigma(copy):
    # structure: matrix brackets
    rd = {rootvec(r): r for r in ROOTS16}
    sigma = {}
    # start: simple roots e_k - e_{k+1} (k<7) and e_6 + e_7 get sigma = +1
    simple = [(k, 1, k+1, -1) for k in range(7)] + [(6, 1, 7, 1)]
    for r in simple: sigma[r] = 1
    def lat_e(r): return V.st_e(lat_root_idx(r, copy))
    def ok_pair(r1, r2):
        pass
    changed = True
    while changed:
        changed = False
        known = list(sigma.items())
        for (r1, s1), (r2, s2) in itertools.product(known, known):
            y = tuple(a + b for a, b in zip(rootvec(r1), rootvec(r2)))
            if y in rd and rd[y] not in sigma:
                r3 = rd[y]
                M = br(Xm(r1), Xm(r2)); d, c = decompose(M)
                N = c.get(r3, 0); assert N != 0 and len(c) == 1 and not d
                L = lattice_bracket(lat_e(r1), lat_e(r2))
                (key, val), = L.items(); assert key == (lat_root_idx(r3, copy), ())
                # psi([X1,X2]) = N sigma3 e^{3} ; [psi X1, psi X2] = s1 s2 val e^{3}
                sigma[r3] = int(Fr(int(s1*s2*val.numerator), int(val.denominator))/N)
                assert sigma[r3] in (1, -1)
                changed = True
            if all(a + b == 0 for a, b in zip(rootvec(r1), rootvec(r2))) and False:
                pass
        # negatives of known roots through [X, X_-] = Cartan
        for r1, s1 in list(sigma.items()):
            y = tuple(-a for a in rootvec(r1)); r2 = rd[y]
            if r2 in sigma: continue
            M = br(Xm(r1), Xm(r2)); d, c = decompose(M); assert not c
            L = lattice_bracket(lat_e(r1), lat_e(r2))
            # psi(sum d_k T_k) = sum d_k eps_k(-1) ; lattice: s1 s2 L
            off = 8*copy
            lat = {k - off: v for (b, m), v in L.items() for (n, k) in m}
            ratio = None
            for kk, dv in d.items():
                rr = Fr(int(lat[kk].numerator), int(lat[kk].denominator))/dv
                ratio = rr if ratio is None else ratio
                assert ratio == rr
            sigma[r2] = int(s1*ratio); assert sigma[r2] in (1, -1)
            changed = True
    assert len(sigma) == len(ROOTS16)
    return sigma
def psi(A, copy, sigma):
    d, c = decompose(A)
    out = {}
    for k, x in d.items():
        V.vadd(out, V.st_h(k + 8*copy), mpq(x.numerator, x.denominator))
    for r, x in c.items():
        V.vadd(out, V.st_e(lat_root_idx(r, copy)), mpq(x.numerator*sigma[r], x.denominator))
    return out
def check_psi(copy, sigma):
    gens = [Tm(k) for k in range(8)] + [Xm(r) for r in ROOTS16]
    bad = 0
    for A, B_ in itertools.combinations(gens, 2):
        if lattice_bracket(psi(A, copy, sigma), psi(B_, copy, sigma)) != psi(br(A, B_), copy, sigma): bad += 1
    return bad

# ------------------------------------------------- small representations (sl2 triples on C^3, C^9)
def sl2_V2():
    """adjoint rep of sl2 on (E, H, F) basis, form B(E,F)=1, B(H,H)=2"""
    E = [[0, -2, 0], [0, 0, 1], [0, 0, 0]]          # ad E: H -> -2E, F -> H
    H = [[2, 0, 0], [0, 0, 0], [0, 0, -2]]
    F = [[0, 0, 0], [-1, 0, 0], [0, 2, 0]]          # ad F: E -> -H, H -> 2F
    Bq = [[0, 0, 1], [0, 2, 0], [1, 0, 0]]
    return [[[Fr(x) for x in row] for row in M] for M in (E, H, F)], [[Fr(x) for x in row] for row in Bq]
def kron(A, B_):
    n, m = len(A), len(B_)
    return [[A[i//m][j//m]*B_[i % m][j % m] for j in range(n*m)] for i in range(n*m)]
def ident(n): return [[Fr(int(i == j)) for j in range(n)] for i in range(n)]
def push(Msmall, images):
    """Msmall acting on small space with basis b_j; images[j] = vector in C^16 of b_j (an isometric embedding).
    Returns the 16x16 matrix acting as Msmall on the image and 0 on its B-orthogonal complement."""
    n = len(images)
    # dual vectors: g_j with B(g_j, image_i) = delta_ij within the image span
    Gram = [[Bf(images[i], images[j]) for j in range(n)] for i in range(n)]
    import sympy as sp
    Ginv = sp.Matrix(Gram).inv()
    # P v = sum_ij image_i Ginv_ij B(image_j, v)
    out = zero()
    for i in range(n):
        for j in range(n):
            # contribution: Msmall[i][j] * image_i (x) dual_j, dual_j(v) = sum_l Ginv[j,l] B(image_l, v)
            if not Msmall[i][j]: continue
            for l in range(n):
                g = Fr(int(sp.fraction(Ginv[j, l])[0]), int(sp.fraction(Ginv[j, l])[1]))
                if not g: continue
                coef = Msmall[i][j]*g
                Bl = [sum(images[l][a]*Bm[a][b] for a in range(N16)) for b in range(N16)]
                for r in range(N16):
                    if images[i][r]:
                        for cc in range(N16):
                            if Bl[cc]: out[r][cc] += coef*images[i][r]*Bl[cc]
    return out
def check_isometry(Bsmall, images, lam):
    n = len(images)
    for i in range(n):
        for j in range(n):
            assert Bf(images[i], images[j]) == lam*Bsmall[i][j], (i, j)

def build_copy():
    """returns dict name -> (E, H, F) 16x16 matrices for (A, B, Cf, Df) =
       (so(3) on C^3 [2A1], root e2+e3 sl2, root e2-e3 + second V2 factor, first V2 factor)"""
    (E3, H3, F3), B3 = sl2_V2()
    u = vec({2: 1, 3: 1}); up = vec({2: 1, 3: -1})
    # C^3: E -> f_0^+, H -> u, F -> lam f_0^- with B(u,u) = 2 = 2 lam  => lam = 1
    imgs3 = [fp(0), u, fm(0)]; check_isometry(B3, imgs3, 1)
    sA = tuple(push(M, imgs3) for M in (E3, H3, F3))
    # C^9 = V2 x V2: basis (a, b) a, b in (E, H, F) ; form B3 x B3 * lam9 ; zero weight (H, H) -> u', B(u',u') = -2 = 4 lam9
    lam9 = Fr(-1, 2)
    B9 = kron(B3, B3)
    # weights: (E,H)->(2,0) f_4 ; (H,E)->(0,2) f_5 ; (E,E)->(2,2) f_6 ; (E,F)->(2,-2) f_7
    names = ['E', 'H', 'F']
    target = {('E', 'H'): fp(4), ('H', 'E'): fp(5), ('E', 'E'): fp(6), ('E', 'F'): fp(7), ('H', 'H'): up}
    imgs9 = [None]*9
    for a in range(3):
        for b in range(3):
            key = (names[a], names[b])
            if key in target: imgs9[3*a + b] = target[key]
    # negative partners: (a', b') with B9 nonzero against the positive one
    for a in range(3):
        for b in range(3):
            if imgs9[3*a + b] is not None: continue
            # partner p with B9[p][3a+b] != 0
            p = [q for q in range(9) if B9[q][3*a + b] != 0][0]
            k = imgs9[p].index(Fr(1))//2          # f_k^+
            imgs9[3*a + b] = [x*lam9*B9[p][3*a + b] for x in fm(k)]
    check_isometry(B9, imgs9, lam9)
    I3 = ident(3)
    first = tuple(push(kron(M, I3), imgs9) for M in (E3, H3, F3))
    second = tuple(push(kron(I3, M), imgs9) for M in (E3, H3, F3))
    # C^4 roots
    def root_sl2(r):
        E = Xm(r); rn = (r[0], -r[1], r[2], -r[3]); Fx = Xm(rn)
        H = br(E, Fx)
        # normalise: [H, E] = 2E
        d, c = decompose(br(H, E)); ce = c[r]
        Fx = [[x/ce*2 for x in row] for row in Fx]
        H = br(E, Fx)
        assert br(H, E) == [[2*x for x in row] for row in E]
        return (E, H, Fx)
    sB = root_sl2((2, 1, 3, 1))
    sC4 = root_sl2((2, 1, 3, -1))
    sC = tuple(madd(x, y) for x, y in zip(sC4, second))
    sD = first
    for (E, H, Fx) in (sA, sB, sC, sD):
        assert br(E, Fx) == H and br(H, E) == madd(E, E) and br(H, Fx) == madd(zero(), Fx, -2)
    for P, Q in itertools.combinations((sA, sB, sC, sD), 2):
        for x in P:
            for y in Q: assert br(x, y) == zero()
    return sA, sB, sC, sD

if __name__ == '__main__':
    sig = [find_sigma(0), find_sigma(1)]
    print('sigma found; bracket check failures copy1, copy2:', check_psi(0, sig[0]), check_psi(1, sig[1]))
    sA, sB, sC, sD = build_copy()
    # copy 1 (index 0): g'_1 = A, g'_2 = B, g''_2 = C, g''_1 = D ; copy 2 (index 1): g'_2 = A, g'_1 = B, g''_1 = C, g''_2 = D
    assign = {0: {'g1p': sA, 'g2p': sB, 'g2pp': sC, 'g1pp': sD}, 1: {'g2p': sA, 'g1p': sB, 'g1pp': sC, 'g2pp': sD}}
    out = {}
    for name in ('g1p', 'g1pp', 'g2p', 'g2pp'):
        trip = []
        for i in range(3):
            v = {}
            for copy in (0, 1):
                V.vadd(v, psi(assign[copy][name][i], copy, sig[copy]))
            trip.append(V.clean(v))
        out[name] = trip
    # checks in the VOA
    names = ['g1p', 'g1pp', 'g2p', 'g2pp']; levels = {'g1p': 3, 'g1pp': 13, 'g2p': 3, 'g2pp': 13}
    for nm in names:
        E, H, F = out[nm]
        assert lattice_bracket(E, F) == H and lattice_bracket(H, E) == V.scale(E, 2) and lattice_bracket(H, F) == V.scale(F, -2)
        lev = V.mode_vec(E, 1, F)
        print(nm, ': sl2 triple OK ; E_(1) F =', dict(lev), ' (expected level', levels[nm], ') ; (H,H)/2 =',
              sum(c*c for (b, m), c in H.items())/2)
        assert lev == {V.VAC: mpq(levels[nm])}
    for a, b in itertools.combinations(names, 2):
        for x in out[a]:
            for y in out[b]:
                assert lattice_bracket(x, y) == {} and V.mode_vec(x, 1, y) == {}
    print('the four sl2 commute and are orthogonal (x_(0) y = x_(1) y = 0)')
    # Dynkin indices per copy: (H restricted to copy k, itself)/2
    for nm in names:
        H = out[nm][1]
        ind = [sum(c*c for (b, m), c in H.items() if (m[0][1] < 8) == (cp == 0))/2 for cp in (0, 1)]
        print('  Dynkin indices of', nm, 'in e8^(1), e8^(2):', [str(x) for x in ind])
    # phi_0 compatibility: phi_0(E) = -F, phi_0(H) = -H ?
    def phi0(v): return {(V.NEGI[b], m): c*(-1)**len(m) for (b, m), c in v.items()}
    for nm in names:
        E, H, F = out[nm]
        print('  phi_0(E) = -F for', nm, ':', phi0(E) == V.scale(F, -1), '; phi_0(H) = -H :', phi0(H) == V.scale(H, -1))
    pickle.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sl2x4.pkl'), 'wb'))
    print('written sl2x4.pkl')
