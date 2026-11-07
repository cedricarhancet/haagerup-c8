"""Exact lattice vertex operator algebra W = V_{E8+E8} in weights <= 2 (R3a, v04-HU).

Written from scratch for this route (the engine of Calculs/cas-II is hard-wired to E6+A2 in a root basis); the
conventions are the same (FLM / Kac), but the Heisenberg part uses the ORTHONORMAL basis eps_0..eps_15 of
h = C^16 (standard coordinates of E8 + E8), which makes Gram = identity.
  Lattice vectors are stored in DOUBLED coordinates (16 integers y, actual vector y/2).  E8 = {y in Z^8 all even
  or all odd, sum y = 0 mod 4}.  Inner product (a, b) = sum a_i b_i / 4.
  States: key = (b, mono), b = index of a lattice vector of norm <= 4 (weight <= 2), mono = sorted tuple of (n, i)
  meaning eps_i(-n).  Vectors: dict key -> coefficient (gmpy2.mpq, or int mod p via the MOD switch).
  Cocycle: eps(a, b) = (-1)^{c(a)^T U c(b)}, c = coordinates in a Z-basis (E8 simple roots, Bourbaki), U upper
  triangular with U_ii = 1, U_ij = G_ij (i < j): bimultiplicative, eps(a,a) = (-1)^{|a|^2/2}, eps(a,-a) = -1 on roots.
  Y(e^a, z) = E^-(-a,z) E^+(-a,z) e_a z^{a(0)}, e_a e^b = eps(a,b) e^{a+b}; general states by the iterate formula
  (a_{-m} w)_k = sum_j C(m+j-1, j) [ a_{-m-j} w_{k+j} - (-1)^m w_{k-m-j} a_j ]   (u_k = u_(k)).
Only products landing in weight <= 2 are allowed (asserted)."""
import itertools
from math import comb
from gmpy2 import mpq

R = 16
# ------------------------------------------------------------------ lattice
def e8_vectors(maxsq):
    """doubled coordinates of E8 vectors with sum y^2 <= maxsq (maxsq = 4 * norm)"""
    out = []
    # all even: y = 2x
    rng = range(-2, 3)
    def rec_even(i, cur, sq):
        if i == 8:
            if sum(cur) % 4 == 0: out.append(tuple(cur))
            return
        for x in rng:
            s = sq + 4*x*x
            if s <= maxsq: rec_even(i + 1, cur + [2*x], s)
    rec_even(0, [], 0)
    def rec_odd(i, cur, sq):
        if i == 8:
            if sum(cur) % 4 == 0: out.append(tuple(cur))
            return
        for x in (-3, -1, 1, 3):
            s = sq + x*x
            if s + (7 - i) <= maxsq: rec_odd(i + 1, cur + [x], s)
    rec_odd(0, [], 0)
    return out
E8_SIMPLE = [(1, -1, -1, -1, -1, -1, -1, 1), (2, 2, 0, 0, 0, 0, 0, 0), (-2, 2, 0, 0, 0, 0, 0, 0),
             (0, -2, 2, 0, 0, 0, 0, 0), (0, 0, -2, 2, 0, 0, 0, 0), (0, 0, 0, -2, 2, 0, 0, 0),
             (0, 0, 0, 0, -2, 2, 0, 0), (0, 0, 0, 0, 0, -2, 2, 0)]
def _dot8(a, b): return sum(x*y for x, y in zip(a, b))
G8 = [[_dot8(a, b)//4 for b in E8_SIMPLE] for a in E8_SIMPLE]
def _inv8():
    import sympy as sp
    M = sp.Matrix([[E8_SIMPLE[j][i] for j in range(8)] for i in range(8)])   # columns = simple roots
    return M.inv()
_MINV = _inv8()
def coords8(y):
    import sympy as sp
    c = _MINV*sp.Matrix(y)
    assert all(x.is_integer for x in c), y
    return tuple(int(x) for x in c)
# U mod 2 as bit masks
U8 = [[(1 if i == j else (G8[i][j] % 2 if i < j else 0)) for j in range(8)] for i in range(8)]

E8V = e8_vectors(16)
E8V.sort(key=lambda y: (_dot8(y, y), y))
assert sum(1 for y in E8V if _dot8(y, y) == 8) == 240 and sum(1 for y in E8V if _dot8(y, y) == 16) == 2160
_C8 = {}
def _bits(y8):
    if y8 not in _C8:
        c = coords8(y8)
        cm = sum((c[i] % 2) << i for i in range(8))
        ucm = 0
        for i in range(8):
            if sum(U8[i][j]*(c[j] % 2) for j in range(8)) % 2: ucm |= 1 << i
        _C8[y8] = (cm, ucm)
    return _C8[y8]
# E8 + E8 vectors of norm <= 4
VECS = []
for a in E8V:
    for b in E8V:
        if _dot8(a, a) + _dot8(b, b) <= 16: VECS.append(a + b)
VECS.sort(key=lambda y: (_dot8(y, y), y))
IDX = {v: i for i, v in enumerate(VECS)}
NSQ = [_dot8(v, v) for v in VECS]           # 4 * norm
ZERO = (0,)*R
Z_IDX = IDX[ZERO]
NEGI = [IDX[tuple(-x for x in v)] for v in VECS]
BITS = []
for v in VECS:
    c1, u1 = _bits(v[:8]); c2, u2 = _bits(v[8:])
    BITS.append((c1 | (c2 << 8), u1 | (u2 << 8)))
def eps_i(ia, ib):
    """cocycle eps(a, b) for lattice indices"""
    ca = BITS[ia][0]; ub = BITS[ib][1]
    return -1 if bin(ca & ub).count('1') % 2 else 1
def ipd(ia, ib):            # 4 * (a, b)
    a, b = VECS[ia], VECS[ib]
    return sum(x*y for x, y in zip(a, b))
def add_i(ia, ib):
    v = tuple(x + y for x, y in zip(VECS[ia], VECS[ib]))
    return IDX.get(v)
ROOT_IDX = [i for i in range(len(VECS)) if NSQ[i] == 8]
N4_IDX = [i for i in range(len(VECS)) if NSQ[i] == 16]
assert len(ROOT_IDX) == 480 and len(N4_IDX) == 61920

# ------------------------------------------------------------------ coefficient ring
MOD = None                  # None: exact mpq ; int p: integers mod p
def C(x):
    if MOD is None: return mpq(x)
    if isinstance(x, int): return x % MOD
    return int(x.numerator) * pow(int(x.denominator), -1, MOD) % MOD
def set_mod(p):
    global MOD
    MOD = p
    clear_caches()
HALF = None

# ------------------------------------------------------------------ states
def weight(key):
    b, m = key
    return NSQ[b]//8 + sum(n for n, i in m)
def partitions_modes(deg):
    out = set()
    def rec(rem, maxpair, cur):
        if rem == 0:
            out.add(tuple(sorted(cur))); return
        for n in range(1, rem + 1):
            for i in range(R):
                p = (n, i)
                if maxpair is not None and p > maxpair: continue
                rec(rem - n, p, cur + [p])
    rec(deg, None, [])
    return sorted(out)
_PM = {}
def monos(deg):
    if deg not in _PM: _PM[deg] = partitions_modes(deg)
    return _PM[deg]
def basis(wt):
    out = []
    for b in range(len(VECS)):
        h = NSQ[b]//8
        if NSQ[b] % 8 or h > wt: continue
        for m in monos(wt - h): out.append((b, m))
    return out

def vadd(acc, v, c=1):
    if MOD is None:
        for k, x in v.items():
            y = acc.get(k)
            y = x*c if y is None else y + x*c
            if y: acc[k] = y
            else: acc.pop(k, None)
    else:
        for k, x in v.items():
            y = (acc.get(k, 0) + x*c) % MOD
            if y: acc[k] = y
            else: acc.pop(k, None)
    return acc
def clean(v):
    if MOD is None: return {k: x for k, x in v.items() if x}
    return {k: x % MOD for k, x in v.items() if x % MOD}
def scale(v, c):
    return clean({k: x*c for k, x in v.items()})
def lin(terms):
    out = {}
    for c, v in terms: vadd(out, v, c)
    return out

def coord(ib, i):           # actual coordinate i of lattice vector ib, as a ring element
    y = VECS[ib][i]
    if MOD is None: return mpq(y, 2)
    return y*pow(2, -1, MOD) % MOD

def heis(i, j, key):
    """eps_i(j) applied to a basis key"""
    b, m = key
    if j < 0:
        return {(b, tuple(sorted(m + ((-j, i),)))): C(1)}
    if j == 0:
        c = VECS[b][i]
        return {key: coord(b, i)} if c else {}
    out = {}
    for t, (n, l) in enumerate(m):
        if n == j and l == i:
            k2 = (b, m[:t] + m[t+1:])
            out[k2] = out.get(k2, 0) + C(j)
    return clean(out)

_SCH = {}
def schur(ia, mdeg):
    """coefficient of z^m in exp(sum_n a(-n) z^n / n): dict mono -> coeff"""
    key = (ia, mdeg)
    r = _SCH.get(key)
    if r is not None: return r
    if mdeg == 0:
        res = {(): C(1)}
    else:
        res = {}
        nz = [(i, coord(ia, i)) for i in range(R) if VECS[ia][i]]
        for n in range(1, mdeg + 1):
            prev = schur(ia, mdeg - n)
            for mono, c in prev.items():
                for i, ai in nz:
                    nm = tuple(sorted(mono + ((n, i),)))
                    res[nm] = res.get(nm, 0) + c*ai
        inv = C(mpq(1, mdeg)) if MOD is None else pow(mdeg, -1, MOD)
        res = clean({k: v*inv for k, v in res.items()})
    _SCH[key] = res
    return res

_EMODE = {}
def emode(ia, k, key):
    """(e^a)_k applied to basis key"""
    ck = (ia, k, key)
    r = _EMODE.get(ck)
    if r is not None: return r
    b, m = key
    out = {}
    ab = add_i(ia, b)
    if ab is None:
        # e^{a+b} of norm > 4: the result has weight >= 3 or vanishes; we never ask for weight > 2
        wres = NSQ[ia]//8 + weight(key) - k - 1
        assert wres <= 2 and True, "weight > 2 requested"
        _EMODE[ck] = {}; return {}
    e = eps_i(ia, b); s4 = ipd(ia, b)
    assert s4 % 4 == 0
    s = s4//4
    nf = len(m)
    cs = [coord(ia, i) for (n, i) in m]
    for T in itertools.product((0, 1), repeat=nf):
        coef = C(e); rest = []; shift = 0; ok = True
        for t in range(nf):
            if T[t]:
                if not VECS[ia][m[t][1]]: ok = False; break
                coef = coef*(-cs[t]); shift += m[t][0]
            else:
                rest.append(m[t])
        if not ok: continue
        mm = -k - 1 - s + shift
        if mm < 0: continue
        rest = tuple(rest)
        for mono, c in schur(ia, mm).items():
            kk = (ab, tuple(sorted(rest + mono)))
            out[kk] = out.get(kk, 0) + coef*c
    out = clean(out)
    _EMODE[ck] = out
    return out

_MODE = {}
def mode(u, k, v):
    """u_(k) v for basis keys u, v"""
    wres = weight(u) + weight(v) - k - 1
    if wres < 0: return {}
    assert wres <= 2, "weight > 2 requested"
    ab = add_i(u[0], v[0])
    if ab is None or NSQ[ab] > 8*wres: return {}      # lattice part e^{a+b} has weight >= |a+b|^2/2
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
        j = 0
        while wtw + wtv - k - j - 1 >= 0:
            x = mode(w, k + j, v)
            if x:
                c = comb(mm + j - 1, j)
                for kk, cc in x.items():
                    vadd(r, heis(i, -mm - j, kk), cc*c)
            j += 1
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

def clear_caches():
    _MODE.clear(); _EMODE.clear(); _SCH.clear()

VAC = (Z_IDX, ())
def one(): return {VAC: C(1)}
def st_h(i, n=1): return {(Z_IDX, ((n, i),)): C(1)}
def st_e(ib): return {(ib, ()): C(1)}
def omega():
    out = {}
    for i in range(R):
        out[(Z_IDX, ((1, i), (1, i)))] = C(mpq(1, 2)) if MOD is None else pow(2, -1, MOD)
    return out
def Lmode(n, v): return mode_vec(omega(), n + 1, v)
def form(u, v):
    """invariant bilinear form (u, v), (1, 1) = 1, u, v homogeneous of the same weight (FHL / engine formula)"""
    if not u or not v: return C(0)
    n = weight(next(iter(u)))
    if weight(next(iter(v))) != n: return C(0)
    tot = C(0); x = dict(u); fact = 1
    for j in range(0, n + 1):
        if not x: break
        r = mode_vec(x, 2*n - j - 1, v).get(VAC, C(0))
        if MOD is None: tot += r / fact
        else: tot = (tot + r*pow(fact, -1, MOD)) % MOD
        x = Lmode(1, x); fact *= (j + 1)
    return tot*(-1)**n if MOD is None else tot*(-1)**n % MOD
