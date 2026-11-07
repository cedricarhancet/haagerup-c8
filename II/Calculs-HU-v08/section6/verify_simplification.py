# Check of the proposed shortcut for m = 2 in Proposition 6.6:
# squares of primitive 40th roots are primitive 20th roots, whose Galois average over Q(sqrt10) is 0.
import sympy as sp
x = sp.symbols('x')
print("mu(20)/phi(20) =", sp.Rational(sp.mobius(20), sp.totient(20)))
from math import gcd
sq = sorted({(2*u) % 40 for u in range(40) if gcd(u,40)==1})
print("2u mod 40 for u unit:", sq, " -> as 20th roots:", sorted({(u) % 20 for u in range(40) if gcd(u,40)==1}))
# order of zeta40^(2u) is 20 for u unit
print(all(sp.gcd(2*u,40)==2 for u in range(40) if gcd(u,40)==1))
for eps in (1,-1):
    val = sp.nsimplify((-4+(2*eps-2)*sp.sqrt(40))/2)
    print("eps",eps,"second moment =", val, " rational part", sp.nsimplify(val).as_coeff_add()[0])
