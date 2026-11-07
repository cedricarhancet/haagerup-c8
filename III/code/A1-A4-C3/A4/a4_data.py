"""a4_data.py -- the finite input data of Dong-Jiao-Xu (arXiv:2611.04778v1), typed from the text
(the original verification files, with the formatted copy DJX26.txt used for the stacked fractions).

All lattice coordinates are in the ordered basis alpha_1..alpha_8 of DJX (4.1); Python tuples are 0-based
(position i-1 <-> alpha_i).  Field constants: r = sqrt13, xi = i sqrt3 (a4_kfield).
"""
from gmpy2 import mpq
from a4_kfield import K, R, XI, RXI, ONE

# ---------------------------------------------------------------- (4.1) Gram matrix of L = E6 + A2
P_EDGES = [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (7, 8)]
G = [[0] * 8 for _ in range(8)]
for i in range(8):
    G[i][i] = 2
for (i, j) in P_EDGES:
    G[i - 1][j - 1] = -1
    G[j - 1][i - 1] = -1

# ---------------------------------------------------------------- (4.4) glue vector
MU = (mpq(4, 3), mpq(5, 3), mpq(2), mpq(4, 3), mpq(2, 3), mpq(1), mpq(2, 3), mpq(1, 3))

# ---------------------------------------------------------------- Section 4.2: Huang
ELL = (-1, -2, -3, -2, -1, -2, 0, 0)
HVEC = (5, 8, 9, 8, 5, 0, 2, 2)
I_SET = [1, 2, 3, 4, 5, 7, 8]          # 1-based
A_VEC = (0, 0, 0, 0, 0, -1, 1, 1)      # a = alpha7 + alpha8 - alpha6
D_VEC = (0, 0, 0, 1, 1, -1, 0, 0)      # d = alpha4 + alpha5 - alpha6
C_VEC = (1, 1, 0, 0, 0, -1, 0, 0)      # c = alpha1 + alpha2 - alpha6
U_COEFFS = [(A_VEC, ONE), (D_VEC, (R + XI * 3) / 20), (C_VEC, (R - XI * 3) / 20)]   # (4.8)

# ---------------------------------------------------------------- Section 4.3: Jiao
RHO1 = (1, 1, 1, 0, 0, 0, 0, 0)
RHO2 = (0, 1, 1, 1, 0, 0, 0, 0)
RHO3 = (0, 0, 1, 1, 1, 0, 0, 0)
Q_VEC = (1, 2, 3, 2, 1, 2, 0, 0)
P_VEC = (1, 2, 3, 2, 1, 0, 0, 0)
B_VEC = (3, 4, 3, 4, 3, 2, 2, 2)


def unit(i):
    """alpha_i, 1-based"""
    return tuple(1 if k == i - 1 else 0 for k in range(8))


def neg(x):
    return tuple(-v for v in x)


# (4.12)-(4.13), formatted copy DJX26.txt l. 570-571:
#   BE = 2e^{a1} + e^{a2} + 2e^{a4} + e^{a5} + e^{q} + 2e^{a7} + e^{a8}
#   BF = -e^{-a1} - 2e^{-a2} - e^{-a4} - 2e^{-a5} - e^{-q} - e^{-a7} - 2e^{-a8}
BE_TERMS = [(unit(1), 2), (unit(2), 1), (unit(4), 2), (unit(5), 1), (Q_VEC, 1), (unit(7), 2), (unit(8), 1)]
BF_TERMS = [(neg(unit(1)), -1), (neg(unit(2)), -2), (neg(unit(4)), -1), (neg(unit(5)), -2), (neg(Q_VEC), -1),
            (neg(unit(7)), -1), (neg(unit(8)), -2)]
PE_TERMS = [(RHO1, 1), (RHO2, 1), (RHO3, 1)]
PF_TERMS = [(neg(RHO1), -1), (neg(RHO2), -1), (neg(RHO3), -1)]

# table of u_{2,6}, u_{1,7} (Section 4.3)
U26_TABLE = [
    ((0, 1, 1, 1, 1, 0, 1, 1), (R + XI) / 6),
    ((1, 1, 1, 1, 0, 0, 1, 1), (R - XI) / 3),
    ((1, 1, 1, 1, 1, 0, 0, 1), (ONE - XI) / 6),
    ((1, 1, 1, 1, 1, 0, 1, 0), (ONE + XI) / 3),
    ((1, 1, 1, 2, 1, 0, 0, 0), (ONE - RXI) / 6),
    ((1, 2, 1, 1, 1, 0, 0, 0), (ONE + RXI) / 12),
    ((2, 3, 4, 3, 2, 2, 0, 0), ONE),
]
U17_TABLE = [
    ((0, 0, -1, 0, 0, -1, 1, 1), XI),
    ((1, 1, 2, 2, 1, 1, 1, 1), R / 2 + XI / 6),
    ((1, 2, 2, 1, 1, 1, 1, 1), R / 4 - XI / 12),
    ((1, 2, 2, 2, 1, 1, 0, 1), ONE / 2 + XI / 6),
    ((1, 2, 2, 2, 1, 1, 1, 0), ONE - XI / 3),
    ((1, 2, 2, 2, 2, 1, 0, 0), ONE / 2 - RXI / 6),
    ((2, 2, 2, 2, 1, 1, 0, 0), ONE + RXI / 3),
]


# ---------------------------------------------------------------- Section 4.4: positive characters
def kappa_H(x):
    """kappa_H(x) = 5^{x1} 8^{x2} 9^{x3} 8^{x4} 5^{x5} 2^{x7} 2^{x8}"""
    return mpq(5) ** x[0] * mpq(8) ** x[1] * mpq(9) ** x[2] * mpq(8) ** x[3] * mpq(5) ** x[4] * mpq(2) ** x[6] * mpq(2) ** x[7]


def kappa_J(x):
    """kappa_J(x) = 2^{-x1 + x2 - x4 + x5 - x7 + x8}"""
    return mpq(2) ** (-x[0] + x[1] - x[3] + x[4] - x[6] + x[7])


# ---------------------------------------------------------------- Section 5
LAMBDA = tuple((pi + bi) * mpq(1, 2) for pi, bi in zip(P_VEC, B_VEC))        # lambda = (p+b)/2
S_VEC = tuple((pi + bi) * mpq(1, 4) for pi, bi in zip(P_VEC, B_VEC))         # s = (p+b)/4
T_VEC = tuple((13 * pi - 3 * bi) * mpq(1, 4) for pi, bi in zip(P_VEC, B_VEC))  # t = (13p-3b)/4


def _e(i):
    return tuple(mpq(v) for v in unit(i))


def _lin(*terms):
    out = [mpq(0)] * 8
    for c, v in terms:
        for k in range(8):
            out[k] += mpq(c) * mpq(v[k])
    return tuple(out)


# (5.4) ordered basis B_H = (2e1, e2-e1, e3+e1, e4-e1, e5-e1, s, e7-e1, e8-e1)  (columns, alpha^J coordinates)
BH_BASIS = [_lin((2, _e(1))), _lin((1, _e(2)), (-1, _e(1))), _lin((1, _e(3)), (1, _e(1))), _lin((1, _e(4)), (-1, _e(1))),
            _lin((1, _e(5)), (-1, _e(1))), S_VEC, _lin((1, _e(7)), (-1, _e(1))), _lin((1, _e(8)), (-1, _e(1)))]
# Section 6.1: basis B of L_0 = B_H with s replaced by 2s = lambda
B_BASIS = list(BH_BASIS)
B_BASIS[5] = LAMBDA

# (5.5) the r_i in the basis B_H
R_IN_BH = [
    (2, 2, 1, 1, 1, -1, 0, 0),
    (2, 1, 1, 1, 0, -1, 1, 1),
    (-6, -3, -2, -3, -2, 3, -2, -2),
    (2, 1, 1, 1, 1, -1, 1, 1),
    (2, 1, 1, 2, 1, -1, 0, 0),
    (3, 1, 0, 1, 1, -2, 1, 1),
    (2, 1, 1, 1, 1, -1, 1, 0),
    (2, 1, 1, 1, 1, -1, 0, 1),
]

# ---------------------------------------------------------------- Section 6.3: chi_* on the columns b_i of B (6.12)
CHI_STAR = [
    ONE / 720,
    (R - XI * 7) / 8,
    (R * (-537) + XI * 1377) / 128,
    (-R + XI * 3) / 2,
    -(ONE * 5 + RXI) / 16,
    -(ONE * 7 + RXI * 5) / 2388787200,
    (R - 3 + (ONE + R) * XI) / 4,
    (-R - 3 + (R - 1) * XI) / 8,
]
# Proposition 6.2 constants
C26 = (ONE * 7 + RXI * 5) / 92160
C17 = XI / 8640
CP = (ONE * 32765 + RXI * 71) / 32768

# ---------------------------------------------------------------- Appendix A.3: nu basis and Tables 2, 3 (Riesen's B_k)
NU = [unit(7), unit(8), unit(1), unit(2), unit(3), unit(4), unit(5)]   # nu_1..nu_7 (1-based in DJX)

# Table 2: quadratic coefficients q^{(k)}_{ij}, rows grouped; values (q0,q1,q2,q3)
TABLE2 = [
    ([(1, 1), (1, 2), (2, 2)], (0, 0, mpq(35, 3), mpq(-35, 6))),
    ([(1, 3), (2, 7)], (0, 0, mpq(-10, 3), mpq(5, 6))),
    ([(1, 4), (2, 6)], (0, 0, mpq(-20, 3), mpq(8, 3))),
    ([(1, 5), (2, 5)], (0, 0, mpq(-9), mpq(9, 2))),
    ([(1, 6), (2, 4)], (0, 0, mpq(-28, 3), mpq(16, 3))),
    ([(1, 7), (2, 3)], (0, 0, mpq(-20, 3), mpq(25, 6))),
    ([(3, 3), (7, 7)], (mpq(1, 3), mpq(4, 9), mpq(5, 3), mpq(-5, 6))),
    ([(3, 4), (4, 4), (6, 6), (6, 7)], (mpq(1, 3), mpq(2, 9), mpq(8, 3), mpq(-4, 3))),
    ([(3, 5), (5, 7)], (0, mpq(-4, 9), mpq(2), mpq(-1))),
    ([(3, 6), (4, 7)], (mpq(-1, 3), mpq(-5, 9), mpq(4, 3), mpq(-2, 3))),
    ([(3, 7)], (mpq(-2, 3), mpq(-1, 9), mpq(2, 3), mpq(-1, 3))),
    ([(4, 5), (5, 6)], (0, mpq(1, 9), mpq(4), mpq(-2))),
    ([(4, 6)], (mpq(-2, 3), mpq(-2, 9), mpq(8, 3), mpq(-4, 3))),
    ([(5, 5)], (0, mpq(1, 3), mpq(3), mpq(-3, 2))),
]


def _nu_comb(*terms):
    out = [0] * 8
    for c, i in terms:
        for k in range(8):
            out[k] += c * NU[i - 1][k]
    return tuple(out)


# Table 3: exponential coefficients g^{(k)}_Gamma, Gamma = set of lattice vectors (in nu's)
TABLE3 = [
    ([_nu_comb((-1, 1), (1, 7)), _nu_comb((-1, 2), (1, 3))], (0, 0, 4, -3)),
    ([_nu_comb((-1, 2), (1, 7)), _nu_comb((-1, 1), (1, 3))], (0, 0, 0, 1)),
    ([_nu_comb((-1, 1), (1, 6)), _nu_comb((-1, 2), (1, 4))], (0, 0, 3, -2)),
    ([_nu_comb((-1, 2), (1, 6)), _nu_comb((-1, 1), (1, 4))], (0, 0, 1, 0)),
    ([_nu_comb((-1, 1), (-1, 2), (1, i), (1, i + 1)) for i in range(3, 7)], (0, 0, -1, 1)),
    ([_nu_comb((1, 5), (-1, 7)), _nu_comb((1, 4), (-1, 7))], (0, mpq(5, 9), 0, 0)),
    ([_nu_comb((-1, 1), (1, 5)), _nu_comb((-1, 2), (1, 5))], (0, 0, 2, -1)),
    ([_nu_comb((1, 4), (-1, 6))], (1, 0, 0, 0)),
    ([_nu_comb((1, 3), (-1, 5))], (0, 1, 0, 0)),
    ([_nu_comb((1, 3), (-1, 6))], (0, mpq(8, 9), 0, 0)),
    ([_nu_comb((1, 3), (-1, 7))], (1, mpq(-1, 3), 0, 0)),
    ([_nu_comb((1, 3), (1, 4), (-1, 6), (-1, 7))], (-1, mpq(8, 9), 0, 0)),
]
