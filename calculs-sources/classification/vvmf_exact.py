"""Exact q-expansion toolkit for vector-valued modular forms for rho = rho_8 (Z(H3), c = 8).

A rho-form is stored as a list V of 12 dicts: V[z][n] = coefficient of q^{Lambda_z + n} (Fractions).
Scalar series are dicts {n: coefficient of q^n} (integer exponents).  Everything is exact.
Coefficients are kept for n <= NT only; all operations below are exact up to that order
(the inputs are known far beyond it: Xi to q^60).
"""
from fractions import Fraction as Fr
import xi_exact as X

d, Lam, names, Xi = X.d, X.Lam, X.names, X.Xi
NT = 6                      # keep coefficients with n <= NT
NS = 40                     # length of scalar series

def sigma(k, n):
    return sum(m**k for m in range(1, n + 1) if n % m == 0)

def smul(a, b, nt=NS):
    c = {}
    for i, x in a.items():
        if x == 0:
            continue
        for j, y in b.items():
            if y == 0 or i + j > nt:
                continue
            c[i + j] = c.get(i + j, 0) + x*y
    return {k: v for k, v in c.items() if v != 0}

def sinv(a, nt=NS):
    """inverse of a scalar series with a[0] != 0 and no negative powers"""
    assert min(a) == 0
    b = {0: Fr(1)/a[0]}
    for n in range(1, nt + 1):
        s = sum(a.get(k, 0)*b[n - k] for k in range(1, n + 1))
        b[n] = -s/a[0]
    return b

E2 = {0: Fr(1)}; E2.update({n: Fr(-24*sigma(1, n)) for n in range(1, NS + 1)})
E4 = {0: Fr(1)}; E4.update({n: Fr(240*sigma(3, n)) for n in range(1, NS + 1)})
E6 = {0: Fr(1)}; E6.update({n: Fr(-504*sigma(5, n)) for n in range(1, NS + 1)})
# Delta/q = prod (1-q^n)^24 ; eta^{-8} q^{1/3} = prod (1-q^n)^{-8}
def eta_power_series(e, nt=NS):
    """prod_{n>=1} (1-q^n)^e as an exact series (e integer, may be negative)"""
    a = {0: Fr(1)}
    for n in range(1, nt + 1):
        # multiply by (1-q^n)^e = sum_k binom(e,k) (-q^n)^k
        f = {0: Fr(1)}
        k = 1; c = Fr(1)
        while n*k <= nt:
            c = c*(e - k + 1)/k
            f[n*k] = c*(-1)**k
            k += 1
        a = smul(a, f, nt)
    return a
P24 = eta_power_series(24)
Dl = {n + 1: v for n, v in P24.items() if n + 1 <= NS}          # Delta
invDlq = sinv(P24)                                              # q/Delta
E4c = smul(smul(E4, E4), E4)
Jser = {n - 1: v for n, v in smul(E4c, invDlq).items()}        # j = E4^3/Delta
Jser[0] -= 744                                                  # J = j - 744
Jser = {k: v for k, v in Jser.items() if v != 0 and k <= NS - 1}
assert Jser[-1] == 1 and 0 not in Jser and Jser[1] == 196884
DJ = {n: n*v for n, v in Jser.items() if n != 0}               # q dJ/dq = -E4^2 E6/Delta
Pm8 = eta_power_series(-8)                                      # eta^{-8} = q^{-1/3} Pm8

# ---------------- rho-forms ----------------
def zero():
    return [dict() for _ in range(d)]

def col(j, nt=NT + 4):
    """the column j of Xi: component z is sum_n Xi[n][z][j] q^{Lam_z+n}, n >= -1"""
    V = zero()
    for n in range(-1, nt + 1):
        for z in range(d):
            v = Xi[n][z][j]
            if v != 0:
                V[z][n] = Fr(v)
    return V

def vadd(A, B, ca=1, cb=1):
    C = zero()
    for z in range(d):
        for n in set(A[z]) | set(B[z]):
            v = ca*A[z].get(n, 0) + cb*B[z].get(n, 0)
            if v != 0:
                C[z][n] = v
    return C

def scal(s, V, nt=NT + 4):
    """scalar series s (integer exponents) times rho-form V"""
    C = zero()
    for z in range(d):
        for n, x in V[z].items():
            for k, y in s.items():
                if n + k > nt:
                    continue
                C[z][n + k] = C[z].get(n + k, 0) + x*y
        C[z] = {k: v for k, v in C[z].items() if v != 0}
    return C

def D(V):
    """q d/dq"""
    return [{n: (Lam[z] + n)*x for n, x in V[z].items() if (Lam[z] + n)*x != 0} for z in range(d)]

def serre(V, k):
    """Serre derivative of weight k: D - k/12 E2"""
    return vadd(D(V), scal(E2, V), 1, Fr(-k, 12))

def Jpow(m):
    s = {0: Fr(1)}
    for _ in range(m):
        s = smul(s, Jser)
    return s

# ---------------- the dual representation ----------------
# rho^* = (rho^T)^{-1} has S^* = S and T^* = conj(T_8).  With pi = (c1 c2)(d1 d5)(d2 d3)(d4 d6)
# (l -> 5l mod 13 on the labels d_l, folded by l ~ -l), the matrix S is pi-invariant and
# conj(T_8)_z = T_16(pi z), where T_16 = exp(-2 pi i 16/24) T_0.  Hence G_z := H_{pi z} is a form for
# rho^* whenever H is a form for rho_16 = rho_8 (x) chi, chi(T) = e^{-2 pi i/3}, chi(S) = 1, and
# H = eta^{-8} Y with Y a weight-four form for rho_8 is such a form of weight zero.
def _pi():
    p = list(range(d))
    def sw(a, b):
        ia, ib = names.index(a), names.index(b)
        p[ia], p[ib] = ib, ia
    sw('c1', 'c2')
    lab = {}
    for l in range(1, 7):
        m = (5*l) % 13
        m = min(m, 13 - m)
        lab[l] = m
    for l in range(1, 7):
        p[names.index('d%d' % l)] = names.index('d%d' % lab[l])
    return p
PI = _pi()

def check_pi_exact():
    """exact checks: (1) -Lambda_z = Lambda_{pi z} - 1/3 mod 1; (2) S is pi-invariant (symbolic, sympy)."""
    for z in range(d):
        x = -Lam[z] - (Lam[PI[z]] - Fr(1, 3))
        assert x.denominator == 1, (names[z], x)
    import sympy as sp
    s13 = sp.sqrt(13)
    x = (13 - 3*s13)/26; y = 3/s13
    def c(j):   # -2y cos(2 pi j/13) depends only on +-j mod 13: canonical representative
        j = j % 13
        return -2*y*sp.cos(2*sp.pi*min(j, 13 - j)/13)
    rows = [[x, 1 - x, 1, 1, 1, 1] + [y]*6, [1 - x, x, 1, 1, 1, 1] + [-y]*6,
            [1, 1, 2, -1, -1, -1] + [0]*6, [1, 1, -1, 2, -1, -1] + [0]*6,
            [1, 1, -1, -1, -1, 2] + [0]*6, [1, 1, -1, -1, 2, -1] + [0]*6]
    for l in range(1, 7):
        rows.append([y, -y, 0, 0, 0, 0] + [c(l*lp) for lp in range(1, 7)])
    bad = 0
    for i in range(d):
        for j in range(d):
            if sp.simplify(rows[PI[i]][PI[j]] - rows[i][j]) != 0:
                bad += 1
    assert bad == 0
    return True

def H_from_Y(Y):
    """H = eta^{-8} Y: returns list of dicts, H[w][n] = coefficient of q^{Lam_w - 1/3 + n}"""
    H = zero()
    for w in range(d):
        for n, x in Y[w].items():
            for k, y in Pm8.items():
                if n + k > NT + 4:
                    continue
                H[w][n + k] = H[w].get(n + k, 0) + x*y
        H[w] = {k: v for k, v in H[w].items() if v != 0}
    return H

def kz(z):
    """integer k_z with Lam_{pi z} - 1/3 = -Lam_z + k_z"""
    x = Lam[PI[z]] - Fr(1, 3) + Lam[z]
    assert x.denominator == 1
    return int(x)

# ---------------- exact linear algebra over Q ----------------
def rref(rows, ncols):
    R = [list(r) for r in rows]
    piv = []
    r = 0
    for c in range(ncols):
        p = None
        for i in range(r, len(R)):
            if R[i][c] != 0:
                p = i; break
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        pv = R[r][c]
        R[r] = [x/pv for x in R[r]]
        for i in range(len(R)):
            if i != r and R[i][c] != 0:
                f = R[i][c]
                R[i] = [a - f*b for a, b in zip(R[i], R[r])]
        piv.append(c); r += 1
        if r == len(R):
            break
    return R[:r], piv

def nullspace(rows, ncols):
    R, piv = rref(rows, ncols)
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for f in free:
        v = [Fr(0)]*ncols
        v[f] = Fr(1)
        for i, p in enumerate(piv):
            v[p] = -R[i][f]
        basis.append(v)
    return basis

def rank(rows, ncols):
    return len(rref(rows, ncols)[0])

# weakly holomorphic scalar forms of negative weight (exact): E4E6/Delta (-2), E4^2/Delta (-4), E6/Delta (-6), E4/Delta (-8)
_qinv = {n - 1: c for n, c in invDlq.items()}           # 1/Delta = q^{-1} * (q/Delta)
E4_over_Delta = smul(E4, _qinv)
E6_over_Delta = smul(E6, _qinv)
E4sq_over_Delta = smul(smul(E4, E4), _qinv)
E10_over_Delta = smul(smul(E4, E6), _qinv)
inv_Delta = dict(_qinv)
