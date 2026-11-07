"""Fast exact arithmetic in K = Q(r, q), r^2 = 13, q^2 = -3, elements 4-tuples of gmpy2.mpq (basis 1, r, q, rq)."""
from gmpy2 import mpq
Z0 = mpq(0); ZERO = (Z0, Z0, Z0, Z0); ONE = (mpq(1), Z0, Z0, Z0)
def mul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return (a0*b0 + 13*a1*b1 - 3*a2*b2 - 39*a3*b3,
            a0*b1 + a1*b0 - 3*(a2*b3 + a3*b2),
            a0*b2 + a2*b0 + 13*(a1*b3 + a3*b1),
            a0*b3 + a3*b0 + a1*b2 + a2*b1)
def add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3])
def sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3])
def neg(a): return (-a[0], -a[1], -a[2], -a[3])
def iszero(a): return not (a[0] or a[1] or a[2] or a[3])
def inv(a):
    # a = x + y r with x = a0 + a2 q, y = a1 + a3 q in Q(q); 1/a = (x - y r)/(x^2 - 13 y^2)
    a0, a1, a2, a3 = a
    # x^2 = (a0^2 - 3 a2^2) + 2 a0 a2 q ; y^2 = (a1^2 - 3 a3^2) + 2 a1 a3 q
    u = a0*a0 - 3*a2*a2 - 13*(a1*a1 - 3*a3*a3)
    v = 2*a0*a2 - 26*a1*a3
    n = u*u + 3*v*v                       # (u + v q)^{-1} = (u - v q)/n
    iu, iv = u/n, -v/n
    # (x - y r)(iu + iv q): x - y r = (a0 + a2 q) - (a1 + a3 q) r
    return mul((a0, -a1, a2, -a3), (iu, Z0, iv, Z0))
