"""a4_common.py -- construction of the DJX generators and of the lattice maps, shared by the a4_item*.py scripts."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import sys
import time
from gmpy2 import mpq
from a4_kfield import K, ONE, ZERO, R, XI, RXI
from a4_voa import Lattice, Vec, LatticeVOA, mat_inverse, mat_det, mat_mul, mat_vec, transpose
import a4_data as D


class Timer(object):
    def __init__(self):
        self.t0 = time.time()

    def __call__(self, msg=""):
        dt = time.time() - self.t0
        return "[%7.2f s] %s" % (dt, msg)


class Report(object):
    """collect (item, name, status, detail) and print; status in {VRAI, FAUX, INFO}"""

    def __init__(self, title):
        self.rows = []
        self.title = title
        self.nfalse = 0
        print("=" * 100)
        print(title)
        print("=" * 100)
        sys.stdout.flush()

    def check(self, name, ok, detail=""):
        status = "VRAI" if ok else "FAUX"
        if not ok:
            self.nfalse += 1
        self.rows.append((name, status, detail))
        line = "%-5s %s" % (status, name)
        if detail:
            line += "   [" + str(detail) + "]"
        print(line)
        sys.stdout.flush()
        return ok

    def info(self, name, detail=""):
        self.rows.append((name, "INFO", detail))
        line = "INFO  %s" % name
        if str(detail) != "":
            line += "   [" + str(detail) + "]"
        print(line)
        sys.stdout.flush()

    def summary(self):
        n = len([r for r in self.rows if r[1] != "INFO"])
        print("-" * 100)
        print("%s : %d checks, %d FAUX" % (self.title, n, self.nfalse))
        sys.stdout.flush()
        return self.nfalse == 0


# ---------------------------------------------------------------- lattice and VOA (one copy of E6+A2 with cocycle 4.2)
def make_voa():
    L = Lattice(D.G)
    return L, LatticeVOA(L)


def expsum(voa, terms):
    out = Vec()
    for (x, c) in terms:
        out += voa.exp(x, c)
    return out


def huang_generators(voa):
    """(4.6)-(4.8): HL, HR, EL, FL, ER, FR, u ; and omega"""
    g = {}
    g["HL"] = voa.heis(D.ELL)
    g["HR"] = voa.heis(D.HVEC)
    g["EL"] = voa.exp(D.ELL)
    g["FL"] = voa.exp(D.neg(D.ELL), -1)
    g["ER"] = expsum(voa, [(D.unit(i), 1) for i in D.I_SET])
    g["FR"] = expsum(voa, [(D.neg(D.unit(i)), -D.HVEC[i - 1]) for i in D.I_SET])
    g["u"] = expsum(voa, D.U_COEFFS)
    g["omega"] = voa.omega()
    return g


def jiao_generators(voa):
    """(4.11)-(4.14) and the table of u_{2,6}, u_{1,7}"""
    g = {}
    g["PE"] = expsum(voa, D.PE_TERMS)
    g["PF"] = expsum(voa, D.PF_TERMS)
    g["PH"] = voa.heis(D.P_VEC)
    g["BE"] = expsum(voa, D.BE_TERMS)
    g["BF"] = expsum(voa, D.BF_TERMS)
    g["BH"] = voa.heis(D.B_VEC)
    g["u26"] = expsum(voa, D.U26_TABLE)
    g["u17"] = expsum(voa, D.U17_TABLE)
    g["omega"] = voa.omega()
    return g


# ---------------------------------------------------------------- Section 5/6 lattice maps
def matrices():
    """return dict with M_BH (columns = B_H basis in alpha^J coords), R (columns r_i in alpha^J coords),
    D (swap 7,8), Rstar = R D, B (basis of L_0), Binv, Rstar_inv"""
    MBH = transpose([list(col) for col in D.BH_BASIS])           # 8x8, column j = BH_BASIS[j]
    RinBH = transpose([list(col) for col in D.R_IN_BH])          # column i = coordinates of r_i in B_H
    Rm = mat_mul(MBH, RinBH)                                     # column i = r_i in alpha^J coordinates
    Dm = [[mpq(1 if i == j else 0) for j in range(8)] for i in range(8)]
    Dm[6][6] = Dm[7][7] = mpq(0)
    Dm[6][7] = Dm[7][6] = mpq(1)
    Rstar = mat_mul(Rm, Dm)
    Bm = transpose([list(col) for col in D.B_BASIS])
    return {"MBH": MBH, "R": Rm, "D": Dm, "Rstar": Rstar, "Rstar_inv": mat_inverse(Rstar), "R_inv": mat_inverse(Rm),
            "B": Bm, "Binv": mat_inverse(Bm)}


def in_L0(x):
    """x in L_0 = {x in L_J : (lambda, x) in 2Z} = even coordinate sum (DJX l. 1639)"""
    return all(mpq(v).denominator == 1 for v in x) and sum(int(v) for v in x) % 2 == 0


def B_coords(mats, x):
    """n with x = B n (requires x in L_0)"""
    n = mat_vec(mats["Binv"], x)
    assert all(v.denominator == 1 for v in n), "x = %r is not in L_0 (B-coordinates %r)" % (x, n)
    return tuple(int(v) for v in n)


def q_star(n):
    """(6.11): q_*(n) = n2 n1 + n5(n1+n2+n3+n4) + n6(n2+n5) + n7(n2+n3+n4+n6) + n8 n1 (mod 2), n 1-based in DJX"""
    n1, n2, n3, n4, n5, n6, n7, n8 = n
    return (n2 * n1 + n5 * (n1 + n2 + n3 + n4) + n6 * (n2 + n5) + n7 * (n2 + n3 + n4 + n6) + n8 * n1) % 2


def gamma_star(mats, x):
    return -1 if q_star(B_coords(mats, x)) else 1


def gamma_from_61(L, mats, Ainv, x):
    """(6.1) with the lattice map A (A^{-1} given): gamma(x) = prod_{i>j} [eps_H(A^{-1}b_i, A^{-1}b_j) eps_J(b_i,b_j)]^{n_i n_j}"""
    n = B_coords(mats, x)
    cols = [tuple(int(v) for v in col) for col in D.B_BASIS]
    pre = []
    for col in cols:
        y = mat_vec(Ainv, col)
        assert all(v.denominator == 1 for v in y)
        pre.append(tuple(int(v) for v in y))
    val = 1
    for i in range(8):
        for j in range(i):
            if (n[i] * n[j]) % 2:
                val *= L.eps(pre[i], pre[j]) * L.eps(cols[i], cols[j])
    return val


def chi_star(mats, x):
    n = B_coords(mats, x)
    val = ONE
    for i in range(8):
        if n[i]:
            val = val * (D.CHI_STAR[i] ** n[i])
    return val


def Phi(voa, mats, v):
    """(6.13): Phi(P e^y) = (R_* P) gamma_*(R_* y) chi_*(R_* y) e^{R_* y}, Huang copy -> Jiao copy (on R_*^{-1} L_0)"""
    def factor(x):
        return K(gamma_star(mats, x)) * chi_star(mats, x)
    return voa.transport(v, mats["Rstar"], factor)


def Phi_inv(voa, mats, v):
    """(8.3): Phi^{-1}(Q e^x) = (R_*^{-1} Q) e^{R_*^{-1} x} / (gamma_*(x) chi_*(x))"""
    def factor(y):
        # y = R_*^{-1} x ; we need x = R_* y
        x = mat_vec(mats["Rstar"], y)
        x = tuple(int(v) for v in x)
        return (K(gamma_star(mats, x)) * chi_star(mats, x)).inverse()
    return voa.transport(v, mats["Rstar_inv"], factor)


def proportional(voa, v, w):
    """return c with v = c w if w != 0 and such c exists in K, else None"""
    if w.is_zero():
        return None
    key = next(iter(w.t))
    if key not in v.t:
        return None
    c = v.t[key] / w.t[key]
    if v == w.scale(c):
        return c
    return None
