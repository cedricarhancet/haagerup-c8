# Lattice E8+E8, roots, norm-4 vectors, cocycle, sl2^4 embedding data.
# All lattice vectors stored DOUBLED (integer tuples of length 16): v2 = 2*v.
import itertools
from fractions import Fraction as Fr

def ip(a, b):
    """inner product of doubled vectors -> exact Fraction (always integer for lattice vectors)"""
    s = sum(x*y for x, y in zip(a, b))
    assert s % 4 == 0, (a, b, s)
    return s // 4

def e8_vectors(norm):
    """all vectors of E8 (doubled coords, length 8) with given norm (2 or 4)"""
    out = []
    # integer type: coords in Z, sum even ; doubled coords even
    # half-integer type: coords in Z+1/2, sum even ; doubled coords odd
    target = 4*norm  # sum of squares of doubled coords
    rng = range(-4, 5)
    for v in itertools.product(rng, repeat=8):
        if sum(x*x for x in v) != target:
            continue
        par = [x % 2 for x in v]
        if all(p == 0 for p in par):
            if (sum(v)//2) % 2 == 0:
                out.append(v)
        elif all(p == 1 for p in par):
            if (sum(v)//2) % 2 == 0:  # sum of coords (=sum(v)/2) even
                out.append(v)
    return out

# Simple roots of E8 (Bourbaki, even coordinate system), doubled
def simple_roots_e8():
    def e(i):
        v = [0]*8; v[i] = 2; return v
    a1 = [1, -1, -1, -1, -1, -1, -1, 1]
    sr = [a1]
    sr.append([x+y for x, y in zip(e(0), e(1))])
    for i in range(1, 7):
        sr.append([x-y for x, y in zip(e(i), e(i-1))])
    return [tuple(s) for s in sr]

SR8 = simple_roots_e8()

def mat_inv_fr(M):
    n = len(M)
    A = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x/pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f*y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]

# coordinates of a vector in simple root basis: v = sum a_i SR[i]  => v^T = a^T S (rows of S = SR)
_S = [list(r) for r in SR8]
_Sinv = mat_inv_fr(_S)  # a^T = v^T S^{-1}

def coords8(v8):
    a = [sum(Fr(v8[k])*_Sinv[k][j] for k in range(8)) for j in range(8)]
    for x in a:
        assert x.denominator == 1
    return [int(x) for x in a]

# Gram matrix of simple roots & cocycle matrix B (mod 2), per copy
_G8 = [[ip(SR8[i], SR8[j]) for j in range(8)] for i in range(8)]
_B8 = [[(_G8[i][j] % 2 if i < j else ((_G8[i][i]//2) % 2 if i == j else 0)) for j in range(8)] for i in range(8)]

def coords16(v):
    return coords8(v[:8]) + coords8(v[8:])

_cache_mask = {}
def masks(v):
    """returns (cmask, bmask): cmask = coords mod 2 bitmask, bmask = (a^T B) mod 2 bitmask"""
    r = _cache_mask.get(v)
    if r is not None:
        return r
    a = coords16(v)
    cm = 0
    bm = 0
    for blk in range(2):
        aa = a[8*blk:8*blk+8]
        for j in range(8):
            if aa[j] % 2:
                cm |= 1 << (8*blk+j)
            s = sum(aa[i]*_B8[i][j] for i in range(8)) % 2
            if s:
                bm |= 1 << (8*blk+j)
    r = (cm, bm)
    _cache_mask[v] = r
    return r

def eps(a, b):
    """bimultiplicative cocycle eps(a,b) = (-1)^{a^T B b}"""
    return -1 if bin(masks(a)[1] & masks(b)[0]).count('1') % 2 else 1

def vadd(a, b):
    return tuple(x+y for x, y in zip(a, b))
def vsub(a, b):
    return tuple(x-y for x, y in zip(a, b))
def vneg(a):
    return tuple(-x for x in a)
def vscale(k, a):
    return tuple(k*x for x in a)

R8 = e8_vectors(2)
N8 = e8_vectors(4)
assert len(R8) == 240 and len(N8) == 2160, (len(R8), len(N8))
Z8 = (0,)*8
ROOTS = [r + Z8 for r in R8] + [Z8 + r for r in R8]
NORM4 = [v + Z8 for v in N8] + [Z8 + v for v in N8] + [a + b for a in R8 for b in R8]
ZERO = (0,)*16

# ---------------- sl2^4 embedding -------------------------------------------
# patterns in ONE E8 copy (undoubled integer coords, positions 1..8 -> index 0..7)
def E(i, s=1):
    v = [0]*8; v[i-1] = s; return v
def add(*vs):
    return [sum(x) for x in zip(*vs)]
# pattern: list of (root, k) with H = sum k*root
PAT = {
 'a': [(add(E(1), E(2, -1)), 1), (add(E(1), E(2)), 1)],                      # index 2, H=2e1
 'b': [(add(E(3), E(6, -1)), 2), (add(E(5), E(6)), 2),
       (add(E(4), E(2, -1)), 1), (add(E(4), E(2)), 1)],                      # index 6, H=2(e3+e4+e5)
 'c': [(add(E(7), E(8)), 1)],                                                 # index 1, H=e7+e8
 'd': [(add(E(7), E(8, -1)), 1), (add(E(6), E(2)), 1), (add(E(6), E(2, -1)), 1),
       (add(E(3), E(4, -1)), 2), (add(E(4), E(5, -1)), 2)],                   # index 7
}
# assignment: sl2 X in order (g'_1, g''_1, g'_2, g''_2)
ASSIGN = [('a', 'c'), ('b', 'd'), ('c', 'a'), ('d', 'b')]   # (pattern in copy1, pattern in copy2)
LEVELS = [3, 13, 3, 13]

def support(X):
    """list of (doubled 16-vector root alpha, k_alpha) for sl2 number X"""
    p1, p2 = ASSIGN[X]
    out = []
    for r, k in PAT[p1]:
        out.append((tuple(2*x for x in r) + Z8, k))
    for r, k in PAT[p2]:
        out.append((Z8 + tuple(2*x for x in r), k))
    return out

def Hvec(X):
    """H_X as doubled 16-vector"""
    h = [0]*16
    for r, k in support(X):
        for i in range(16):
            h[i] += k*r[i]
    return tuple(h)
