# exact arithmetic in K = Q(sqrt(-2), sqrt(3)) ; basis 1, s2=sqrt(-2), s3=sqrt(3), s6=s2*s3=sqrt(-6)
from fractions import Fraction as Fr
class Q4:
    __slots__ = ('a',)
    def __init__(self, a=0, b=0, c=0, d=0):
        self.a = (Fr(a), Fr(b), Fr(c), Fr(d))
    @staticmethod
    def lift(x):
        return x if isinstance(x, Q4) else Q4(x)
    def __add__(s, o):
        o = Q4.lift(o); return Q4(*[x+y for x, y in zip(s.a, o.a)])
    __radd__ = __add__
    def __neg__(s): return Q4(*[-x for x in s.a])
    def __sub__(s, o): return s + (-Q4.lift(o))
    def __rsub__(s, o): return Q4.lift(o) - s
    def __mul__(s, o):
        if not isinstance(o, Q4):
            o = Fr(o); return Q4(*[x*o for x in s.a])
        a1, b1, c1, d1 = s.a; a2, b2, c2, d2 = o.a
        # products: s2^2=-2, s3^2=3, s6^2=-6, s2 s3=s6, s2 s6=-2 s3, s3 s6=3 s2
        A = a1*a2 - 2*b1*b2 + 3*c1*c2 - 6*d1*d2
        B = a1*b2 + b1*a2 + 3*(c1*d2 + d1*c2)
        C = a1*c2 + c1*a2 - 2*(b1*d2 + d1*b2)
        D = a1*d2 + d1*a2 + b1*c2 + c1*b2
        return Q4(A, B, C, D)
    __rmul__ = __mul__
    def conj(s):  # complex conjugation: s2 -> -s2, s3 -> s3, s6 -> -s6
        a, b, c, d = s.a; return Q4(a, -b, c, -d)
    def __truediv__(s, o):
        o = Fr(o); return Q4(*[x/o for x in s.a])
    def __bool__(s): return any(s.a)
    def __eq__(s, o):
        o = Q4.lift(o); return s.a == o.a
    def __hash__(s): return hash(s.a)
    def __repr__(s):
        names = ['', '*sqrt(-2)', '*sqrt(3)', '*sqrt(-6)']
        t = [('%s%s' % (x, n)) for x, n in zip(s.a, names) if x]
        return '(' + ' + '.join(t) + ')' if t else '0'
    def __complex__(s):
        import math
        a, b, c, d = map(float, s.a)
        return complex(a + c*math.sqrt(3), b*math.sqrt(2) + d*math.sqrt(6))
