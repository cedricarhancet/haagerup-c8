"""Reduction modulo a prime p of the number field K = Q(sqrt13, sqrt-3) (note r = sqrt(6)sqrt(-7+sqrt13) = sqrt(-3)(sqrt13-1)), containing the
coefficients of the lines of case II (candidate 9).  p is chosen so that 13, -3 and 6(-7 + s13_p) are squares
mod p; an embedding K -> F_p is fixed by the choice of square roots (another choice gives a Galois conjugate
configuration, which is also a solution).  Ranks over F_p of p-integral vectors are lower bounds for ranks over K;
if a family of vectors has full rank r over F_p and over K, membership tests mod p are faithful."""
import sympy as sp
from gmpy2 import mpq
def find_prime(start=2**31 - 1):
    p = sp.prevprime(start)
    while True:
        if sp.legendre_symbol(13, p) == 1 and sp.legendre_symbol(p - 3, p) == 1:
            s13 = sp.sqrt_mod(13, p)
            s3 = sp.sqrt_mod(p - 3, p)
            # r = sqrt(6) sqrt(-7 + sqrt13) (principal branches) = sqrt(-3) (sqrt13 - 1): K = Q(sqrt13, sqrt-3)
            r = s3*(s13 - 1) % p
            assert (r*r - (-42 + 6*s13)) % p == 0
            return p, s13, s3, r
        p = sp.prevprime(p)
class Fp:
    def __init__(self, start=2**31 - 1):
        self.p, self.s13, self.s3, self.r = find_prime(start)
    def q(self, x):           # mpq or sympy Rational -> F_p
        n, d = (int(x.numerator), int(x.denominator)) if hasattr(x, 'numerator') else (int(x.p), int(x.q))
        return n*pow(d, -1, self.p) % self.p
    def alg(self, expr):
        """sympy expression in sqrt(13), sqrt(3)*I, sqrt(39)*I, sqrt(6)*sqrt(-7+sqrt(13)) -> F_p"""
        S13, S3, R = sp.symbols('S13 S3 R')
        e = sp.sympify(expr)
        e = e.subs(sp.sqrt(6)*sp.sqrt(-7 + sp.sqrt(13)), R)
        e = e.subs(sp.sqrt(-7 + sp.sqrt(13)), R/sp.sqrt(6))
        e = e.subs(sp.sqrt(39), S13*S3/sp.I).subs(sp.sqrt(3), S3/sp.I).subs(sp.sqrt(13), S13)
        e = sp.expand(e)
        assert not e.has(sp.I) and not e.has(sp.sqrt(6)) and not any(isinstance(a, sp.Pow) and a.exp == sp.Rational(1, 2) for a in sp.preorder_traversal(e)), e
        P = sp.Poly(e, S13, S3, R)
        tot = 0
        for (a, b, c), coef in P.terms():
            tot += self.q(sp.Rational(coef))*pow(self.s13, a, self.p)*pow(self.s3, b, self.p)*pow(self.r, c, self.p)
        return tot % self.p
    def vec(self, v, idx, n):
        out = [0]*n
        for k, c in v.items(): out[idx[k]] = (out[idx[k]] + self.q(c)) % self.p
        return out
