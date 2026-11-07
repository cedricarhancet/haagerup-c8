"""a4_voa.py -- exact lattice vertex operator algebra engine, conventions of Dong-Jiao-Xu (DJX, arXiv:2611.04778v1).

Independent second implementation (written for the A4 verification, 5 Nov 2026).  Nothing is imported from the
Article-I engine (the original verification files); the conventions are those of
DJX Section 4.1 and Appendix A.1:

  * ordered basis alpha_1..alpha_n of a positive definite even lattice M, Gram matrix G, coordinates x = (x_1..x_n),
    (x, y) = x^T G y ; internally index i-1 represents alpha_i ;
  * bimultiplicative cocycle (4.2):  eps(x, y) = (-1)^{ sum_i x_i y_i + sum_{i>j} x_i G_ij y_j } ;
  * Fock monomial  P e^x, P a product of oscillators alpha_i(-k), k >= 1 ; weight = x^2/2 + sum k ;
  * Y(e^x, z) e^y = eps(x, y) z^{(x,y)} exp( sum_{n>=1} x(-n) z^n / n ) e^{x+y}  and the general formulas
    (A.1) (exponential acting on a monomial) and (A.2) (oscillators on the first vector, iterate formula) ;
  * Hermitian form (A.3) with a positive lattice character kappa : (P e^x, Q e^y)_kappa = kappa(x) delta_{x,y} <P,Q>,
    linear in the first variable, coefficients of the second vector conjugated ;
  * antilinear involution theta_kappa (h_1(-n_1)...h_m(-n_m) e^x) = (-1)^m kappa(x) h_1(-n_1)...h_m(-n_m) e^{-x}.

A vector is a dict {(x, mono): K} where x is a tuple of n ints and mono a sorted tuple of pairs (i, k) meaning
alpha_{i+1}(-k) (0-based i), with repetitions.  Coefficients are in K = Q(sqrt13, sqrt-3) (a4_kfield.K).
"""
from gmpy2 import mpq
from math import comb
from functools import lru_cache
from a4_kfield import K, ONE, ZERO

__all__ = ["Lattice", "Vec", "LatticeVOA"]


# ----------------------------------------------------------------------------------------------------------------
# exact linear algebra helpers (rational)
def mat_inverse(M):
    """exact inverse of a square rational matrix (list of lists), Gauss-Jordan with mpq"""
    n = len(M)
    A = [[mpq(M[i][j]) for j in range(n)] + [mpq(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            raise ValueError("singular matrix")
        A[col], A[piv] = A[piv], A[col]
        p = A[col][col]
        A[col] = [v / p for v in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return [row[n:] for row in A]


def mat_det(M):
    n = len(M)
    A = [[mpq(M[i][j]) for j in range(n)] for i in range(n)]
    det = mpq(1)
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return mpq(0)
        if piv != col:
            A[col], A[piv] = A[piv], A[col]
            det = -det
        p = A[col][col]
        det *= p
        for r in range(col + 1, n):
            if A[r][col] != 0:
                f = A[r][col] / p
                A[r] = [a - f * b for a, b in zip(A[r], A[col])]
    return det


def mat_mul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(mpq(A[i][k]) * mpq(B[k][j]) for k in range(m)) for j in range(p)] for i in range(n)]


def mat_vec(A, v):
    return tuple(sum(mpq(A[i][k]) * mpq(v[k]) for k in range(len(v))) for i in range(len(A)))


def transpose(A):
    return [list(r) for r in zip(*A)]


# ----------------------------------------------------------------------------------------------------------------
class Lattice(object):
    """positive definite even lattice with ordered basis and DJX cocycle (4.2)"""

    def __init__(self, G):
        self.n = len(G)
        self.G = [[int(v) for v in row] for row in G]
        for i in range(self.n):
            assert self.G[i][i] % 2 == 0, "lattice must be even"
            for j in range(self.n):
                assert self.G[i][j] == self.G[j][i]
        self.Ginv = mat_inverse(self.G)
        self.det = mat_det(self.G)
        self.zero = (0,) * self.n

    def ip(self, x, y):
        G = self.G
        s = 0
        for i, xi in enumerate(x):
            if xi:
                Gi = G[i]
                for j, yj in enumerate(y):
                    if yj:
                        s += xi * Gi[j] * yj
        return s

    def ipb(self, x, i):
        """(x, alpha_{i+1})"""
        Gi = self.G[i]
        return sum(xj * Gi[j] for j, xj in enumerate(x) if xj)

    def norm(self, x):
        return self.ip(x, x)

    def eps(self, x, y):
        """DJX cocycle (4.2): (-1)^{ sum_i x_i y_i + sum_{i>j} x_i G_ij y_j }"""
        G = self.G
        e = 0
        for i, xi in enumerate(x):
            if xi:
                e += xi * y[i]
                Gi = G[i]
                for j in range(i):
                    if y[j]:
                        e += xi * Gi[j] * y[j]
        return -1 if e % 2 else 1

    @staticmethod
    def add(x, y):
        return tuple(a + b for a, b in zip(x, y))

    @staticmethod
    def sub(x, y):
        return tuple(a - b for a, b in zip(x, y))

    @staticmethod
    def neg(x):
        return tuple(-a for a in x)

    @staticmethod
    def smul(c, x):
        return tuple(c * a for a in x)


# ----------------------------------------------------------------------------------------------------------------
def _insert_sorted(mono, pair):
    lst = list(mono)
    lst.append(pair)
    lst.sort()
    return tuple(lst)


class Vec(object):
    """sparse vector of V_M with coefficients in K"""
    __slots__ = ("t",)

    def __init__(self, t=None):
        self.t = {} if t is None else t

    def copy(self):
        return Vec(dict(self.t))

    def is_zero(self):
        return not self.t

    def __bool__(self):
        return bool(self.t)

    def __len__(self):
        return len(self.t)

    def addterm(self, key, c):
        if not c:
            return
        t = self.t
        if key in t:
            s = t[key] + c
            if s:
                t[key] = s
            else:
                del t[key]
        else:
            t[key] = c

    def __iadd__(self, other):
        for key, c in other.t.items():
            self.addterm(key, c)
        return self

    def __add__(self, other):
        out = self.copy()
        out += other
        return out

    def __neg__(self):
        return Vec({k: -c for k, c in self.t.items()})

    def __sub__(self, other):
        out = self.copy()
        for key, c in other.t.items():
            out.addterm(key, -c)
        return out

    def scale(self, c):
        if not isinstance(c, K):
            c = K(c)
        if not c:
            return Vec()
        return Vec({k: v * c for k, v in self.t.items()})

    def __mul__(self, c):
        return self.scale(c)

    __rmul__ = __mul__

    def __eq__(self, other):
        return (self - other).is_zero()

    def __ne__(self, other):
        return not self.__eq__(other)

    def conj_coeffs(self):
        return Vec({k: v.conj() for k, v in self.t.items()})

    def momenta(self):
        return sorted(set(k[0] for k in self.t))

    def keys_sorted(self):
        return sorted(self.t.keys())

    def __repr__(self):
        if not self.t:
            return "0"
        parts = []
        for (x, mono) in self.keys_sorted():
            c = self.t[(x, mono)]
            osc = "".join("a%d(-%d)" % (i + 1, k) for (i, k) in mono)
            parts.append("[%s] %se^%s" % (c, osc, x))
        return "\n".join(parts)


# ----------------------------------------------------------------------------------------------------------------
class LatticeVOA(object):
    def __init__(self, lattice):
        self.L = lattice
        self.n = lattice.n
        self._schur_cache = {}

    # ---- constructors
    def vacuum(self):
        return Vec({(self.L.zero, ()): ONE})

    def exp(self, x, c=ONE):
        if not isinstance(c, K):
            c = K(c)
        return Vec({(tuple(int(v) for v in x), ()): c})

    def heis(self, h, k=1, x=None, c=ONE):
        """h(-k) e^x for h = sum h_j alpha_j given by rational coordinates h (tuple), default x = 0"""
        if x is None:
            x = self.L.zero
        x = tuple(int(v) for v in x)
        if not isinstance(c, K):
            c = K(c)
        out = Vec()
        for j, hj in enumerate(h):
            if hj:
                out.addterm((x, ((j, k),)), c * K(mpq(hj)))
        return out

    def heis_create(self, h, k, v):
        """h(-k) v for h with rational coordinates, k >= 1"""
        out = Vec()
        for (x, mono), c in v.t.items():
            for j, hj in enumerate(h):
                if hj:
                    out.addterm((x, _insert_sorted(mono, (j, k))), c * K(mpq(hj)))
        return out

    def heis_annih(self, h, j, v):
        """h(j) v for j >= 0 (h rational coordinates)"""
        L = self.L
        out = Vec()
        if j == 0:
            for (x, mono), c in v.t.items():
                val = sum(mpq(hi) * L.ipb(x, i) for i, hi in enumerate(h) if hi)
                if val:
                    out.addterm((x, mono), c * K(val))
            return out
        # j > 0: h(j) a_{l,j} = j (h, alpha_l)
        hG = [sum(mpq(hi) * L.G[i][l] for i, hi in enumerate(h) if hi) for l in range(L.n)]
        for (x, mono), c in v.t.items():
            seen = set()
            for pos, (l, k) in enumerate(mono):
                if k != j or (l, k) in seen:
                    continue
                seen.add((l, k))
                mult = mono.count((l, k))
                coef = j * hG[l] * mult
                if coef:
                    lst = list(mono)
                    lst.remove((l, k))
                    out.addterm((x, tuple(lst)), c * K(coef))
        return out

    def heis_mode(self, h, j, v):
        if j < 0:
            return self.heis_create(h, -j, v)
        return self.heis_annih(h, j, v)

    # ---- weights
    def weight_key(self, key):
        x, mono = key
        return mpq(self.L.norm(x), 2) + sum(k for (_, k) in mono)

    def weights(self, v):
        return sorted(set(self.weight_key(k) for k in v.t))

    def max_weight(self, v):
        if not v.t:
            return -1
        return max(self.weight_key(k) for k in v.t)

    def is_homogeneous(self, v, w):
        return all(self.weight_key(k) == w for k in v.t)

    # ---- Schur polynomials E_x[m]
    def schur(self, x, m):
        """E_x[m] = coefficient of t^m in exp(sum_{k>=1} x(-k) t^k / k), as dict mono -> mpq"""
        key = (x, m)
        if key in self._schur_cache:
            return self._schur_cache[key]
        if m < 0:
            return {}
        # polynomial in t: list of dicts (degree 0..m)
        poly = [dict() for _ in range(m + 1)]
        poly[0][()] = mpq(1)
        for k in range(1, m + 1):
            # factor exp(x(-k) t^k / k) = sum_j (x(-k)/k)^j / j! t^{kj}
            lin = {((i, k),): mpq(xi, k) for i, xi in enumerate(x) if xi}
            if not lin:
                continue
            powers = [{(): mpq(1)}]
            jmax = m // k
            for j in range(1, jmax + 1):
                prev = powers[-1]
                nxt = {}
                for mono1, c1 in prev.items():
                    for mono2, c2 in lin.items():
                        mm = _insert_sorted(mono1, mono2[0])
                        nxt[mm] = nxt.get(mm, mpq(0)) + c1 * c2 / j   # divide by j to build 1/j!
                powers.append(nxt)
            newpoly = [dict() for _ in range(m + 1)]
            for deg, pd in enumerate(poly):
                if not pd:
                    continue
                for j in range(0, (m - deg) // k + 1):
                    pj = powers[j]
                    tgt = newpoly[deg + k * j]
                    for mono1, c1 in pd.items():
                        for mono2, c2 in pj.items():
                            mm = tuple(sorted(mono1 + mono2))
                            tgt[mm] = tgt.get(mm, mpq(0)) + c1 * c2
            poly = newpoly
        res = {mm: c for mm, c in poly[m].items() if c}
        self._schur_cache[key] = res
        return res

    # ---- the mode formula
    def exp_mode(self, x, n, b):
        """(e^x)_n b by (A.1)"""
        L = self.L
        out = Vec()
        for (y, mono), cb in b.t.items():
            d0 = -n - 1 - L.ip(x, y)
            e = L.eps(x, y)
            xy = L.add(x, y)
            r = len(mono)
            cvals = [L.ipb(x, i) for (i, k) in mono]
            for mask in range(1 << r):
                coefS = 1
                kS = 0
                remaining = []
                for pos in range(r):
                    if mask >> pos & 1:
                        if cvals[pos] == 0:
                            coefS = 0
                            break
                        coefS *= -cvals[pos]
                        kS += mono[pos][1]
                    else:
                        remaining.append(mono[pos])
                if coefS == 0:
                    continue
                m = d0 + kS
                if m < 0:
                    continue
                Em = self.schur(x, m)
                base = cb * K(e * coefS)
                for monoE, ce in Em.items():
                    newmono = tuple(sorted(remaining + list(monoE)))
                    out.addterm((xy, newmono), base * K(ce))
        return out

    def mode_mono(self, x, mono, n, b):
        """(mono e^x)_n b by (A.2) recursively"""
        if not mono:
            return self.exp_mode(x, n, b)
        (i, k) = mono[0]
        rest = mono[1:]
        wa = self.weight_key((x, rest))
        wb = self.max_weight(b)
        out = Vec()
        ei = tuple(1 if l == i else 0 for l in range(self.n))
        # first sum: m >= 0, a'_{n+m} b nonzero requires n+m <= wa + wb - 1
        mmax = int(wa + wb - 1 - n)
        for m in range(0, mmax + 1):
            coef = comb(m + k - 1, k - 1)
            inner = self.mode_mono(x, rest, n + m, b)
            if inner:
                out += self.heis_create(ei, k + m, inner).scale(coef)
        # second sum: j >= 0 with alpha_i(j) b nonzero
        js = {0}
        for (y, mb) in b.t:
            for (_, kk) in mb:
                js.add(kk)
        sign = (-1) ** (k - 1)
        for j in sorted(js):
            hb = self.heis_annih(ei, j, b)
            if hb:
                coef = sign * comb(j + k - 1, k - 1)
                inner = self.mode_mono(x, rest, n - j - k, hb)
                if inner:
                    out += inner.scale(coef)
        return out

    def mode(self, a, n, b):
        """a_n b"""
        out = Vec()
        for (x, mono), ca in a.t.items():
            res = self.mode_mono(x, mono, n, b)
            if res:
                out += res.scale(ca)
        return out

    def normal(self, a, b):
        """:ab: = a_{-1} b"""
        return self.mode(a, -1, b)

    def deriv(self, a):
        """L(-1) a = a_{-2} 1"""
        return self.mode(a, -2, self.vacuum())

    # ---- Hermitian form (A.3)
    def _pair_mono(self, P, Q):
        """sum over bijections of prod n_a delta G_{i_a j_pi(a)} (rational)"""
        if len(P) != len(Q):
            return mpq(0)
        if not P:
            return mpq(1)
        (i, n_) = P[0]
        rest = P[1:]
        total = mpq(0)
        Ql = list(Q)
        used = set()
        for pos, (j, m_) in enumerate(Ql):
            if m_ != n_ or (j, m_) in used:
                continue
            used.add((j, m_))
            g = self.L.G[i][j]
            if g == 0:
                continue
            mult = Ql.count((j, m_))
            Qrest = list(Ql)
            Qrest.remove((j, m_))
            total += n_ * g * mult * self._pair_mono(rest, tuple(Qrest))
        return total

    def herm(self, v, w, kappa):
        """(v, w)_kappa, linear in v, antilinear in w; kappa: x -> positive rational"""
        total = ZERO
        # group by momentum
        byx = {}
        for (x, mono), c in w.t.items():
            byx.setdefault(x, []).append((mono, c))
        for (x, mono), cv in v.t.items():
            if x not in byx:
                continue
            kx = K(mpq(kappa(x)))
            for (monow, cw) in byx[x]:
                p = self._pair_mono(mono, monow)
                if p:
                    total = total + cv * cw.conj() * kx * K(p)
        return total

    def theta(self, v, kappa):
        """antilinear involution theta_kappa"""
        out = Vec()
        for (x, mono), c in v.t.items():
            sgn = (-1) ** len(mono)
            out.addterm((self.L.neg(x), mono), c.conj() * K(sgn * mpq(kappa(x))))
        return out

    # ---- conformal vector (4.3)
    def omega(self):
        out = Vec()
        Ginv = self.L.Ginv
        for i in range(self.n):
            for j in range(self.n):
                c = Ginv[i][j] / 2
                if c:
                    out.addterm((self.L.zero, tuple(sorted(((i, 1), (j, 1))))), K(c))
        return out

    # ---- lattice character involutions exp(2 pi i (v, .)) on momenta
    def momentum_phase_check(self, v, h, denom):
        """return the set of values (h, x) * denom (rational) over momenta x of v -- for involution checks"""
        return sorted(set(mpq(sum(mpq(hi) * self.L.ipb(x, i) for i, hi in enumerate(h))) * denom for x in v.momenta()))

    # ---- transport by a lattice map (rational matrix A acting on coordinates) with scalar factor on exponentials
    def transport(self, v, A, factor):
        """v -> sum c * (A h_1)(-n_1)...(A h_m)(-n_m) factor(A y) e^{A y}.  A: n x n rational matrix (list of lists)
        acting on coordinate columns; factor: function (target momentum tuple) -> K."""
        out = Vec()
        cols = [tuple(mpq(A[r][c]) for r in range(self.n)) for c in range(self.n)]
        for (y, mono), c in v.t.items():
            Ay = mat_vec(A, y)
            assert all(val.denominator == 1 for val in Ay), "transport: image momentum not integral %r" % (Ay,)
            Ay = tuple(int(val) for val in Ay)
            f = factor(Ay)
            base = Vec({(Ay, ()): c * f})
            for (i, k) in reversed(mono):
                base = self.heis_create(cols[i], k, base)
            out += base
        return out


# ----------------------------------------------------------------------------------------------------------------
def sl2_checks(voa, E, F, H, k, label=""):
    """return list of (name, bool) for E_0F=H, H_0E=2E, H_0F=-2F, E_1F=k1, H_1H=2k1, E_0E=0, F_0F=0, E_1E=0, F_1F=0"""
    vac = voa.vacuum()
    res = []
    res.append(("%sE_0F = H" % label, voa.mode(E, 0, F) == H))
    res.append(("%sH_0E = 2E" % label, voa.mode(H, 0, E) == E.scale(2)))
    res.append(("%sH_0F = -2F" % label, voa.mode(H, 0, F) == F.scale(-2)))
    res.append(("%sE_1F = k1 (k=%s)" % (label, k), voa.mode(E, 1, F) == vac.scale(k)))
    res.append(("%sH_1H = 2k1" % label, voa.mode(H, 1, H) == vac.scale(2 * k)))
    res.append(("%sE_0E = 0" % label, voa.mode(E, 0, E).is_zero()))
    res.append(("%sF_0F = 0" % label, voa.mode(F, 0, F).is_zero()))
    res.append(("%sE_1E = 0" % label, voa.mode(E, 1, E).is_zero()))
    res.append(("%sF_1F = 0" % label, voa.mode(F, 1, F).is_zero()))
    res.append(("%sH_1E = 0, H_1F = 0, E_1H = 0" % label,
                voa.mode(H, 1, E).is_zero() and voa.mode(H, 1, F).is_zero() and voa.mode(E, 1, H).is_zero()))
    res.append(("%sH_0H = 0" % label, voa.mode(H, 0, H).is_zero()))
    return res


def commute_checks(voa, triple1, triple2, label=""):
    """all a_n b = 0 for a in triple1, b in triple2, n = 0, 1 (weight-one vectors: higher modes vanish)"""
    ok = True
    for a in triple1:
        for b in triple2:
            for n in (0, 1):
                if not voa.mode(a, n, b).is_zero():
                    ok = False
    return [("%scommutation of the two triples (all a_n b = 0, n=0,1)" % label, ok)]


def sugawara(voa, E, F, H, k):
    """Omega_k(E,F,H) = [ (1/2) H_{-1}H + E_{-1}F + F_{-1}E ] / (2(k+2))   (4.10)"""
    v = voa.normal(H, H).scale(mpq(1, 2)) + voa.normal(E, F) + voa.normal(F, E)
    return v.scale(mpq(1, 2 * (k + 2)))


def eta_k(voa, E, F, H, k):
    """(7.1): eta_k = [ -k H_{-2}1 - H_{-1}H + 2k E_{-1}F ] / (2k(k+2))"""
    v = voa.deriv(H).scale(-k) - voa.normal(H, H) + voa.normal(E, F).scale(2 * k)
    return v.scale(mpq(1, 2 * k * (k + 2)))


def virasoro_checks(voa, w, label=""):
    """w_0 w = dw, w_1 w = 2w, w_2 w = 0, w_3 w = (c/2) 1 ; returns (list of checks, central charge c or None)"""
    res = []
    res.append(("%sw_0 w = L(-1)w" % label, voa.mode(w, 0, w) == voa.deriv(w)))
    res.append(("%sw_1 w = 2w" % label, voa.mode(w, 1, w) == w.scale(2)))
    res.append(("%sw_2 w = 0" % label, voa.mode(w, 2, w).is_zero()))
    w3 = voa.mode(w, 3, w)
    c = None
    if len(w3.t) == 1 and (voa.L.zero, ()) in w3.t:
        c = w3.t[(voa.L.zero, ())] * 2
    res.append(("%sw_3 w = (c/2) 1" % label, c is not None))
    return res, c


def casimir(voa, E, F, H, v):
    """C = H_0^2 + 2 E_0 F_0 + 2 F_0 E_0 applied to v"""
    out = voa.mode(H, 0, voa.mode(H, 0, v))
    out += voa.mode(E, 0, voa.mode(F, 0, v)).scale(2)
    out += voa.mode(F, 0, voa.mode(E, 0, v)).scale(2)
    return out
