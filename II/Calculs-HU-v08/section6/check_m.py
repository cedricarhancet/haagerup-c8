import sympy as sp, math
from sympy import factorint, totient, sqrt, Rational
# squarefree decomposition of 9m^2+4
viol = []; D10 = []
minratio = (10**9, None)
for m in range(1, 20001):
    Delta = 9*m*m+4
    f = factorint(Delta)
    D = 1; t = 1
    for p,e in f.items():
        if e % 2: D *= p
        t *= p**(e//2)
    assert t*t*D == Delta
    assert D % 3 == 1 and D > 1
    for p in factorint(D):
        assert p == 2 or p % 4 == 1
    if D == 10:
        D10.append((m,t)); continue
    r = float(totient(2*D))/math.sqrt(D)
    if r < minratio[0]: minratio = (r, (m, D))
    if r <= 12/math.sqrt(13): viol.append((m,D,r))
print("m with D=10 (m<=20000):", D10)
print("min phi(2D)/sqrt D for D!=10:", minratio, " 12/sqrt13 =", 12/math.sqrt(13))
print("violations of Huang's phi(2D)>=12 sqrt D/sqrt13:", viol[:5])
# D=10 bound
print("12/sqrt10 - 7/2 =", 12/math.sqrt(10)-3.5, " 2/(that) =", 2/(12/math.sqrt(10)-3.5))
# check inequality 4m^2+2 >= m^2/2 + 4mt for D=10 m values
for (m,t) in D10:
    print(m, t, 4*m*m+2, Rational(m*m,2)+4*m*t)
# small D counterexamples quoted in Remark 6.12(3)
for D in (2,5,7):
    print(D, totient(2*D), 12*math.sqrt(D)/math.sqrt(13))
# Huang identity
m = sp.symbols('m')
print(sp.factor(36*m**2*(9*m**2+4) - 13*(4*m**2+2)**2))
# values (p-1)/sqrt p
print(4/math.sqrt(5), 12/math.sqrt(13))
# m=2 second moment, eps=-1
print("eps=-1 rhs:", sp.nsimplify((-4+(2*(-1)-2)*sp.sqrt(40))/2))
