"""Modular data of Z(H3) (Drinfeld centre of the Haagerup fusion category H3).

Source: Evans--Gannon, CMP 307 (2011), arXiv:1006.1326v1, equations (2.5a)-(2.5b), p. 20.
Order of primaries: 0, b, a, c0, c1, c2, d1, ..., d6.
S is real symmetric; T is given at central charge c = 0 mod 24.
At central charge c we use T_c = exp(-2 pi i c/24) T.
"""
import mpmath as mp

mp.mp.dps = 60

def zh3_S():
    s13 = mp.sqrt(13)
    x = (13 - 3*s13)/26
    y = 3/s13
    c = lambda j: -2*y*mp.cos(2*mp.pi*j/13)
    M = [[mp.mpf(0)]*12 for _ in range(12)]
    row0 = [x, 1-x, 1, 1, 1, 1] + [y]*6
    rowb = [1-x, x, 1, 1, 1, 1] + [-y]*6
    rowa = [1, 1, 2, -1, -1, -1] + [0]*6
    rowc0 = [1, 1, -1, 2, -1, -1] + [0]*6
    rowc1 = [1, 1, -1, -1, -1, 2] + [0]*6
    rowc2 = [1, 1, -1, -1, 2, -1] + [0]*6
    rows = [row0, rowb, rowa, rowc0, rowc1, rowc2]
    for l in range(1, 7):
        rows.append([y, -y, 0, 0, 0, 0] + [c(l*lp) for lp in range(1, 7)])
    S = mp.matrix(12, 12)
    for i in range(12):
        for j in range(12):
            S[i, j] = mp.mpf(rows[i][j])/3
    return S

def zh3_T(c=0):
    xi3 = mp.exp(2j*mp.pi/3)
    xi13 = mp.exp(2j*mp.pi/13)
    diag = [1, 1, 1, 1, xi3, xi3**-1] + [xi13**(6*l*l) for l in range(1, 7)]
    ph = mp.exp(-2j*mp.pi*mp.mpf(c)/24)
    T = mp.matrix(12, 12)
    for i in range(12):
        T[i, i] = ph*diag[i]
    return T

def check(c=8):
    S = zh3_S(); T = zh3_T(c)
    I = mp.eye(12)
    S2 = S*S
    ST3 = (S*T)**3
    e1 = mp.mnorm(S2 - I, 1)
    # (ST)^3 = S^2 up to the phase coming from c: for a modular category, (ST)^3 = p S^2 with p = exp(2 pi i c_top/8) e^{-2pi i c/8}
    # We report (ST)^3 S^{-2}
    R = ST3*mp.inverse(S2)
    ph = R[0, 0]
    e2 = mp.mnorm(R - ph*I, 1)
    unit = mp.mnorm(S*S.H - I, 1)
    dims = [S[0, j]/S[0, 0] for j in range(12)]
    return e1, e2, ph, unit, dims

if __name__ == '__main__':
    for c in (0, 8):
        e1, e2, ph, unit, dims = check(c)
        print('c =', c, ' |S^2-I| =', mp.nstr(e1, 5), ' |(ST)^3 S^-2 - p I| =', mp.nstr(e2, 5),
              ' p =', mp.nstr(ph, 12), ' |SS*-I| =', mp.nstr(unit, 5))
    print('quantum dimensions:', [mp.nstr(d, 12) for d in dims])
    d = (3+mp.sqrt(13))/2
    print('d =', mp.nstr(d, 12), ' 3d =', mp.nstr(3*d, 12), ' 3d+1 =', mp.nstr(3*d+1, 12), ' 3d+2=', mp.nstr(3*d+2, 12))
