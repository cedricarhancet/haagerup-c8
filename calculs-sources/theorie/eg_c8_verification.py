"""Independent numerical verification of the Evans--Gannon classification of
character vectors for the modular data of Z(H3) at central charge 8.

Input: Lambda and chi of Evans--Gannon (arXiv:1006.1326v1, equation (2.7), p. 21-22),
transcribed below.  We
 (1) build Xi(tau) = q^Lambda (q^{-1} I + chi + sum_{n>=1} Xi_n q^n) from the
     Bantay--Gannon differential equation  D Xi = Xi [ (Lambda - I) F + (chi + [Lambda,chi]) G ],
     F = (J-240) Delta/E10,  G = Delta/E10,  D = q d/dq  (recursion in exact rationals);
 (2) test Xi(-1/tau) = S Xi(tau) at sample points (columns are vector-valued modular functions);
 (3) [the fundamental-matrix property is checked separately, with the exact
     Bantay--Gannon trace criterion, in bg_trace_criterion.py; an earlier determinant
     test was removed because det Xi * Delta^7 has weight 84 and is not constant];
 (4) derive the admissible character vectors (positivity: leading exponent of component z != 0
     strictly larger than -c/24 = -1/3) and print the first coefficients.
"""
from fractions import Fraction as Fr
import mpmath as mp
from zh3_modular_data import zh3_S, zh3_T

mp.mp.dps = 50
NMAX = 60

Lam = [Fr(2,3)]*4 + [Fr(0), Fr(1,3), Fr(5,39), Fr(20,39), Fr(32,39), Fr(2,39), Fr(8,39), Fr(11,39)]
chi = [
 [6, 80, 81, 81, 8748, 1215, 3549, 273, 13, 5538, 2275, 1378],
 [80, 6, 81, 81, 8748, 1215, -3549, -273, -13, -5538, -2275, -1378],
 [81, 81, 167, -81, -8748, -1215, 0, 0, 0, 0, 0, 0],
 [81, 81, -81, 167, -8748, -1215, 0, 0, 0, 0, 0, 0],
 [3, 3, -3, -3, -12, 18, 0, 0, 0, 0, 0, 0],
 [27, 27, -27, -27, 1458, -152, 0, 0, 0, 0, 0, 0],
 [7, -7, 0, 0, 0, 0, -88, -14, -1, 50, 63, 64],
 [42, -42, 0, 0, 0, 0, -1484, 92, 16, 2940, -192, -1041],
 [119, -119, 0, 0, 0, 0, -2142, 987, 11, -24990, -6035, 4641],
 [5, -5, 0, 0, 0, 0, 17, 13, -3, -2, 35, -14],
 [13, -13, 0, 0, 0, 0, 174, -1, -5, 294, -147, 51],
 [14, -14, 0, 0, 0, 0, 448, -77, 7, -343, 125, -24]]
names = ['0','b','a','c0','c1','c2','d1','d2','d3','d4','d5','d6']
d = 12

# ---------- q-series for E4, E6, Delta, J (exact integers) ----------
def sigma(k, n):
    return sum(m**k for m in range(1, n+1) if n % m == 0)
def series_mul(a, b, N):
    c = [0]*(N+1)
    for i, ai in enumerate(a[:N+1]):
        if ai == 0: continue
        for j, bj in enumerate(b[:N+1-i]):
            c[i+j] += ai*bj
    return c
def series_inv(a, N):  # a[0] != 0, exact Fractions
    b = [Fr(0)]*(N+1)
    b[0] = Fr(1, a[0])
    for n in range(1, N+1):
        s = sum(a[k]*b[n-k] for k in range(1, n+1))
        b[n] = -s/a[0]
    return b
N = NMAX + 5
E4 = [1] + [240*sigma(3, n) for n in range(1, N+1)]
E6 = [1] + [-504*sigma(5, n) for n in range(1, N+1)]
# Delta = q prod (1-q^n)^24, stored as series in q starting at q^1: Delta = q * Dl
Dl = [1] + [0]*N
for n in range(1, N+1):
    for _ in range(24):
        new = Dl[:]
        for k in range(n, N+1):
            new[k] -= Dl[k-n]
        Dl = new
E10 = series_mul(E4, E6, N)
invE10 = series_inv([Fr(x) for x in E10], N)
# G = Delta/E10 = q * (Dl * invE10)
GG = series_mul([Fr(x) for x in Dl], invE10, N)
g = [Fr(0)] + GG[:N]           # g[k] = coefficient of q^k, g[0] = 0
# J = E4^3/Delta - 744 ; J*Delta/E10 = E4^3/E10 - 744 G ; F = J G - 240 G = E4^3/E10 - 984 G
E4c = series_mul(series_mul(E4, E4, N), E4, N)
F1 = series_mul([Fr(x) for x in E4c], invE10, N)
f = [F1[k] - 984*g[k] for k in range(N+1)]
assert f[0] == 1 and f[1] == 0, (f[0], f[1])

# ---------- recursion ----------
LmI = [[(Lam[i]-1 if i == j else Fr(0)) for j in range(d)] for i in range(d)]
chiF = [[Fr(x) for x in r] for r in chi]
comm = [[(Lam[i]-Lam[j])*chiF[i][j] for j in range(d)] for i in range(d)]  # [Lam,chi]_{ij}
B = [[chiF[i][j] + comm[i][j] for j in range(d)] for i in range(d)]
def matmul(A, Bm):
    return [[sum(A[i][k]*Bm[k][j] for k in range(d) if A[i][k] != 0) for j in range(d)] for i in range(d)]
I = [[Fr(1) if i == j else Fr(0) for j in range(d)] for i in range(d)]
Xi = {-1: I, 0: chiF}
for n in range(1, NMAX+1):
    R = [[Fr(0)]*d for _ in range(d)]
    for l in range(-1, n):
        k = n - l
        M = [[f[k]*LmI[i][j] + g[k]*B[i][j] for j in range(d)] for i in range(d)]
        P = matmul(Xi[l], M)
        for i in range(d):
            for j in range(d):
                R[i][j] += P[i][j]
    Xi[n] = [[R[i][j]/(Lam[i]-Lam[j]+n+1) for j in range(d)] for i in range(d)]
    # consistency at n = 0 is automatic; check integrality
nonint = sum(1 for n in range(1, NMAX+1) for i in range(d) for j in range(d) if Xi[n][i][j].denominator != 1)
print('non-integral entries among Xi_1..Xi_%d:' % NMAX, nonint)

def Xi_num(tau):
    q = mp.exp(2j*mp.pi*tau)
    M = mp.matrix(d, d)
    for i in range(d):
        qi = mp.exp(2j*mp.pi*tau*mp.mpf(Lam[i].numerator)/Lam[i].denominator)
        for j in range(d):
            s = mp.mpc(0)
            qn = 1/q
            for n in range(-1, NMAX+1):
                v = Xi[n][i][j]
                if v != 0:
                    s += mp.mpf(v.numerator)/v.denominator*qn
                qn *= q
            M[i, j] = qi*s
    return M

S = zh3_S()
worst = 0
for tau in [mp.mpc(0.05, 1.02), mp.mpc(-0.3, 0.98), mp.mpc(0.2, 1.3), mp.mpc(0.45, 0.9), mp.mpc(0, 1)]:
    X1 = Xi_num(-1/tau); X2 = S*Xi_num(tau)
    err = mp.mnorm(X1 - X2, 1)/mp.mnorm(X2, 1)
    worst = max(worst, err)
    print('tau =', mp.nstr(tau, 4), ' relative |Xi(-1/tau) - S Xi(tau)| =', mp.nstr(err, 5))
print('worst S-relation error:', mp.nstr(worst, 5))

trL = sum(Lam)
print('tr Lambda =', trL, '(compare bg_trace_criterion.py)')

# admissible character vectors: X = Xi v with v_0 = 1 and v_z free iff Lam_z - 1 > -1/3 (z != 0)
free = [z for z in range(1, d) if Lam[z] - 1 > Fr(-1, 3)]
print('free columns (positivity allows a pole column):', [names[z] for z in free])
for z in range(d):
    coeffs = []
    for n in range(-1, 4):
        a0 = Xi[n][z][0]
        a1 = [Xi[n][z][j] for j in free]
        coeffs.append((int(a0), [int(x) for x in a1]))
    print(names[z].ljust(3), 'q^%s *' % str(Lam[z]), coeffs)
