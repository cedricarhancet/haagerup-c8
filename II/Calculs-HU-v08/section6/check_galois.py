import cmath, math
from math import gcd
from fractions import Fraction
import sympy as sp

# Haagerup: moments for m=1
z = [cmath.exp(2j*math.pi*6*l*l/13) for l in range(1,7)]
print("sum theta =", sum(z), " expected", -(1+math.sqrt(13))/2)
print("sum theta^2 =", sum(w*w for w in z), " expected eps=1:", (-1+(2-1)*math.sqrt(13))/2)

def mobius(n): return sp.mobius(n)
def phi(n): return sp.totient(n)

# <zeta>_10 computed by brute force: L = lcm(n,40); G_L = {u in (Z/L)^x : sigma_u(sqrt10)=sqrt10}
# sigma_u(sqrt10) = chi40(u) sqrt10, where chi40 = Kronecker symbol (40/u)? compute via explicit sqrt10 in Q(zeta40)
# sqrt10 = sqrt2*sqrt5, sqrt2 = z8+z8^-1, sqrt5 = sum_{k} (k/5) z5^k
def sigma_sqrt10_sign(u):
    z40 = cmath.exp(2j*math.pi/40)
    def s10(e):
        z8 = z40**(5*e); z5 = z40**(8*e)
        s2 = z8 + 1/z8
        s5 = z5 - z5**2 - z5**3 + z5**4
        return s2*s5
    v = s10(u % 40); w = s10(1)
    return round((v/w).real)
def avg(n):
    L = n*40//gcd(n,40)
    G = [u for u in range(1,L) if gcd(u,L)==1 and sigma_sqrt10_sign(u)==1]
    zeta = lambda u: cmath.exp(2j*math.pi*u*(L//n)/L)  # zeta = e^{2pi i/n}, sigma_u(zeta)= zeta^u
    s = sum(zeta(u) for u in G)/len(G)
    return s
bad = []
for n in range(1, 400):
    s = avg(n)
    A = s.real  # A + B sqrt10 is real
    assert abs(s.imag) < 1e-9
    if n % 40 != 0:
        exp = float(mobius(n))/float(phi(n))
        if abs(s.real-exp) > 1e-9: bad.append((n,'part1',s,exp))
    else:
        k = n//40
        if gcd(k,40) > 1:
            if abs(s) > 1e-9: bad.append((n,'part2',s))
        else:
            # A = 0, |B| = |r(k)|/(8 phi(k))
            B = abs(s.real)/math.sqrt(10)
            exp = abs(float(mobius(k)))/(8*float(phi(k)))
            if abs(B-exp) > 1e-9: bad.append((n,'part3',s,exp))
print("galois-average mismatches:", bad)
# check |A|+8|B| <=1, but A,B need separation: for case 3 A=0 so fine; case1 B=0.
# coset C_- for n=40
G40 = [u for u in range(1,40) if gcd(u,40)==1 and sigma_sqrt10_sign(u)==1]
nG40 = [u for u in range(1,40) if gcd(u,40)==1 and sigma_sqrt10_sign(u)==-1]
sG = sum(cmath.exp(2j*math.pi*u/40) for u in G40)
snG = sum(cmath.exp(2j*math.pi*u/40) for u in nG40)
print("G40 =", G40, "sum", sG)
print("other coset =", nG40, "sum", snG)
Cm = G40 if sG.real < 0 else nG40
print("C_- =", Cm, " mod 20:", sorted(u%20 for u in Cm))
print("21 in G40?", 21 in G40)
print("sum over C_- of zeta20^u:", sum(cmath.exp(2j*math.pi*u/20) for u in Cm))
# linear independence of the 8 roots of C_-: rank of their coordinates in power basis of Q(zeta40) (degree 16)
x = sp.symbols('x')
Phi40 = sp.cyclotomic_poly(40, x)
rows = []
for u in Cm:
    r = sp.Poly(sp.rem(x**u, Phi40, x), x).all_coeffs()[::-1]
    r = r + [0]*(16-len(r))
    rows.append(r)
print("rank of C_- roots:", sp.Matrix(rows).rank())
