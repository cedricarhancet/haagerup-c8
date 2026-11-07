"""Exact lattice vertex operator algebra V_L for L = E6 + A2 (rank 8), up to small weights.

Conventions (FLM / Kac).  L has basis the simple roots alpha_0..alpha_5 (E6, Bourbaki order, node 2 attached
to node 4; here indices 0..5) and alpha_6, alpha_7 (A2), Gram matrix G = diag(C_E6, C_A2).
States: e^beta (x) monomial in alpha_i(-n), n >= 1.  A basis state is a key (beta, mono) with beta a tuple of
8 integers (simple-root coordinates) and mono a sorted tuple of pairs (n, i) meaning alpha_i(-n).
Vectors are dicts key -> gmpy2.mpq.
Vertex operators:
   [h(m), h'(n)] = m (h,h') delta_{m+n,0},  h(0) e^beta = (h,beta) e^beta,
   Y(e^a, z) = exp(sum_{n>0} a(-n) z^n / n) exp(-sum_{n>0} a(n) z^{-n} / n) e_a z^{a(0)},
   e_a e^b = eps(a,b) e^{a+b},  eps bimultiplicative, eps(a_i,a_j) = (-1)^{G_ij} (i<j), (-1)^{G_ii/2} (i=j), 1 (i>j).
Modes of general states by the iterate formula with a = alpha_i(-1)1:
   (a_{-m} w)_k = sum_{j>=0} C(m+j-1, j) [ a_{-m-j} w_{k+j} - (-1)^m w_{k-m-j} a_j ].
"""
from gmpy2 import mpq
from math import comb
import itertools
CE6 = [[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]
G = [[0]*8 for _ in range(8)]
for i in range(6):
    for j in range(6):
        G[i][j] = CE6[i][j]
G[6][6] = G[7][7] = 2; G[6][7] = G[7][6] = -1
R = 8
ZERO = (0,)*R
def ip(a, b):
    return sum(a[i]*G[i][j]*b[j] for i in range(R) for j in range(R) if a[i] and b[j])
def ipb(a, i):      # (a, alpha_i)
    return sum(a[j]*G[j][i] for j in range(R) if a[j])
def eps(a, b):
    e = 0
    for i in range(R):
        if not a[i]: continue
        e += a[i]*b[i]
        for j in range(i+1, R):
            e += a[i]*b[j]*G[i][j]
    return -1 if e % 2 else 1
def add(a, b): return tuple(x+y for x, y in zip(a, b))
def neg(a): return tuple(-x for x in a)

# lattice vectors of norm <= 2N by BFS over roots
def lattice_vectors(maxnorm):
    unit = [tuple(1 if k == i else 0 for k in range(R)) for i in range(R)]
    roots = set()
    fr = list(unit); roots.update(fr)
    while fr:
        new = []
        for r in fr:
            for i in range(R):
                c = ipb(r, i)
                s = tuple(r[k] - (c if k == i else 0) for k in range(R))
                if s not in roots: roots.add(s); new.append(s)
        fr = new
    roots = sorted(roots)
    vecs = {ZERO}; fr = [ZERO]
    while fr:
        new = []
        for v in fr:
            for r in roots:
                w = add(v, r)
                if w not in vecs and ip(w, w) <= maxnorm:
                    vecs.add(w); new.append(w)
        fr = new
    return roots, sorted(vecs, key=lambda v: (ip(v, v), v))
ROOTS, VECS = lattice_vectors(8)
assert len(ROOTS) == 78
NORM = {v: ip(v, v) for v in VECS}

def partitions_modes(deg):
    """multisets of (n, i) with sum n = deg, as sorted tuples"""
    out = []
    def rec(rem, maxpair, cur):
        if rem == 0:
            out.append(tuple(sorted(cur))); return
        for n in range(1, rem+1):
            for i in range(R):
                p = (n, i)
                if maxpair is not None and p > maxpair: continue
                rec(rem-n, p, cur+[p])
    rec(deg, None, [])
    return sorted(set(out))
_PM = {}
def monos(deg):
    if deg not in _PM: _PM[deg] = partitions_modes(deg)
    return _PM[deg]
def weight(key):
    b, m = key
    return NORM[b]//2 + sum(n for n, i in m)
def basis(wt):
    out = []
    for b in VECS:
        h = NORM[b]//2
        if h > wt: continue
        for m in monos(wt - h):
            out.append((b, m))
    return out

def vadd(acc, v, c=1):
    for k, x in v.items():
        y = acc.get(k)
        y = x*c if y is None else y + x*c
        if y: acc[k] = y
        else: acc.pop(k, None)
    return acc
def clean(v): return {k: x for k, x in v.items() if x}

# Heisenberg action of alpha_i(j) on a basis key
def heis(i, j, key):
    b, m = key
    if j < 0:
        return {(b, tuple(sorted(m + ((-j, i),)))): mpq(1)}
    if j == 0:
        c = ipb(b, i)
        return {key: mpq(c)} if c else {}
    out = {}
    for t, (n, l) in enumerate(m):
        if n == j and G[i][l]:
            nm = m[:t] + m[t+1:]
            k2 = (b, nm)
            out[k2] = out.get(k2, 0) + mpq(j*G[i][l])
    return clean(out)
def heis_vec(i, j, v):
    out = {}
    for k, c in v.items():
        vadd(out, heis(i, j, k), c)
    return out

# Schur polynomials p_m(a) in the a(-n), expanded in the basis alpha_i(-n)
_SCH = {}
def schur(a, m):
    """dict mono -> coeff, for the coefficient of z^m in exp(sum a(-n) z^n / n)"""
    key = (a, m)
    if key in _SCH: return _SCH[key]
    if m == 0:
        res = {(): mpq(1)}
    else:
        # m p_m = sum_{n=1}^m a(-n) p_{m-n}
        res = {}
        for n in range(1, m+1):
            prev = schur(a, m-n)
            for mono, c in prev.items():
                for i in range(R):
                    if a[i]:
                        nm = tuple(sorted(mono + ((n, i),)))
                        res[nm] = res.get(nm, 0) + c*a[i]
        res = {k: v/m for k, v in res.items() if v}
    _SCH[key] = res
    return res
def mono_mul(m1, m2): return tuple(sorted(m1 + m2))

_EMODE = {}
def emode(a, k, key):
    """(e^a)_k applied to basis key"""
    ck = (a, k, key)
    if ck in _EMODE: return _EMODE[ck]
    b, m = key
    out = {}
    if NORM[add(a, b)] if add(a, b) in NORM else True:
        pass
    ab = add(a, b)
    if ab not in NORM:
        _EMODE[ck] = {}; return {}
    e = eps(a, b); s = ip(a, b)
    nf = len(m)
    cs = [ipb(a, i) for (n, i) in m]
    for T in itertools.product((0, 1), repeat=nf):
        coef = mpq(e); rest = []; shift = 0; ok = True
        for t in range(nf):
            if T[t]:
                if cs[t] == 0: ok = False; break
                coef *= -cs[t]; shift += m[t][0]
            else:
                rest.append(m[t])
        if not ok: continue
        mm = -k - 1 - s + shift
        if mm < 0: continue
        rest = tuple(rest)
        for mono, c in schur(a, mm).items():
            kk = (ab, mono_mul(rest, mono))
            out[kk] = out.get(kk, 0) + coef*c
    out = clean(out)
    _EMODE[ck] = out
    return out

_MODE = {}
def mode(u, k, v):
    """u_k v for basis keys u, v"""
    wres = weight(u) + weight(v) - k - 1
    if wres < 0: return {}
    ck = (u, k, v)
    r = _MODE.get(ck)
    if r is not None: return r
    a, m = u
    if not m:
        r = emode(a, k, v)
    else:
        (mm, i) = m[0]
        w = (a, m[1:])
        r = {}
        wtw, wtv = weight(w), weight(v)
        # first sum: a_{-mm-j} (w_{k+j} v)
        j = 0
        while wtw + wtv - k - j - 1 >= 0:
            x = mode(w, k + j, v)
            if x:
                c = comb(mm + j - 1, j)
                for kk, cc in x.items():
                    vadd(r, heis(i, -mm - j, kk), cc*c)
            j += 1
        # second sum: -(-1)^mm w_{k-mm-j} (a_j v), j >= 0
        sgn = -1 if mm % 2 == 0 else 1
        maxj = max([n for n, l in v[1]] + [0])
        for j in range(0, maxj + 1):
            y = heis(i, j, v)
            if not y: continue
            c = sgn*comb(mm + j - 1, j)
            for kk, cc in y.items():
                vadd(r, mode(w, k - mm - j, kk), cc*c)
    _MODE[ck] = r
    return r

def mode_vec(u, k, v):
    out = {}
    for ku, cu in u.items():
        for kv, cv in v.items():
            x = mode(ku, k, kv)
            if x: vadd(out, x, cu*cv)
    return out

VAC = (ZERO, ())
def state_h(i, n=1):        # alpha_i(-n) 1
    return {(ZERO, ((n, i),)): mpq(1)}
def state_e(a):
    return {(a, ()): mpq(1)}

# conformal vector omega = 1/2 sum G^{-1}_{ij} alpha_i(-1) alpha_j(-1) 1
import fractions
def _ginv():
    import sympy as sp
    M = sp.Matrix(G).inv()
    return [[mpq(int(sp.fraction(M[i, j])[0]), int(sp.fraction(M[i, j])[1])) for j in range(R)] for i in range(R)]
GINV = _ginv()
def omega():
    out = {}
    for i in range(R):
        for j in range(R):
            if GINV[i][j]:
                k = (ZERO, tuple(sorted(((1, i), (1, j)))))
                out[k] = out.get(k, 0) + GINV[i][j]/2
    return clean(out)
OMEGA = omega()
def Lmode(n, v):            # L_n = omega_{n+1}
    return mode_vec(OMEGA, n + 1, v)

def vac_coeff(v):
    return v.get(VAC, mpq(0))
def form(u, v):
    """invariant bilinear form <u, v> with <1,1> = 1 (u, v homogeneous of the same weight)"""
    if not u or not v: return mpq(0)
    n = weight(next(iter(u)))
    if weight(next(iter(v))) != n: return mpq(0)
    tot = mpq(0); x = dict(u); fact = 1
    for j in range(0, n + 1):
        if not x: break
        tot += vac_coeff(mode_vec(x, 2*n - j - 1, v)) / mpq(fact)
        x = Lmode(1, x); fact *= (j + 1)
    return tot*(-1)**n
