"""a4_kfield.py -- exact arithmetic in K = Q(sqrt(13), sqrt(-3)).

An element is a + b*r + c*xi + d*r*xi with a,b,c,d in Q (gmpy2.mpq), where
r = sqrt(13) > 0, xi = i*sqrt(3), so r^2 = 13, xi^2 = -3, (r xi)^2 = -39.
Complex conjugation: r -> r, xi -> -xi.
Written for the independent verification A4 of Dong-Jiao-Xu (arXiv:2611.04778v1).
"""
from gmpy2 import mpq

__all__ = ["K", "ZERO", "ONE", "R", "XI", "RXI", "kq"]


def _q(x):
    if isinstance(x, mpq):
        return x
    return mpq(x)


class K(object):
    __slots__ = ("a", "b", "c", "d")

    def __init__(self, a=0, b=0, c=0, d=0):
        self.a = _q(a)
        self.b = _q(b)
        self.c = _q(c)
        self.d = _q(d)

    # ---- basic predicates
    def is_zero(self):
        return self.a == 0 and self.b == 0 and self.c == 0 and self.d == 0

    def is_rational(self):
        return self.b == 0 and self.c == 0 and self.d == 0

    def is_real(self):
        return self.c == 0 and self.d == 0

    def __bool__(self):
        return not self.is_zero()

    def __eq__(self, other):
        if not isinstance(other, K):
            try:
                other = K(other)
            except Exception:
                return False
        return self.a == other.a and self.b == other.b and self.c == other.c and self.d == other.d

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.a, self.b, self.c, self.d))

    # ---- arithmetic
    def __add__(self, other):
        if not isinstance(other, K):
            other = K(other)
        return K(self.a + other.a, self.b + other.b, self.c + other.c, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return K(-self.a, -self.b, -self.c, -self.d)

    def __sub__(self, other):
        if not isinstance(other, K):
            other = K(other)
        return K(self.a - other.a, self.b - other.b, self.c - other.c, self.d - other.d)

    def __rsub__(self, other):
        return K(other).__sub__(self)

    def __mul__(self, other):
        if not isinstance(other, K):
            o = _q(other)
            return K(self.a * o, self.b * o, self.c * o, self.d * o)
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = other.a, other.b, other.c, other.d
        return K(a1 * a2 + 13 * b1 * b2 - 3 * c1 * c2 - 39 * d1 * d2,
                 a1 * b2 + b1 * a2 - 3 * (c1 * d2 + d1 * c2),
                 a1 * c2 + c1 * a2 + 13 * (b1 * d2 + d1 * b2),
                 a1 * d2 + d1 * a2 + b1 * c2 + c1 * b2)

    __rmul__ = __mul__

    def conj(self):
        """complex conjugation xi -> -xi"""
        return K(self.a, self.b, -self.c, -self.d)

    def conj_r(self):
        """Galois conjugation r -> -r (auxiliary, used for the inverse)"""
        return K(self.a, -self.b, self.c, -self.d)

    def norm_sq(self):
        """|z|^2 = z * conj(z), an element of Q(r)"""
        return self * self.conj()

    def inverse(self):
        if self.is_zero():
            raise ZeroDivisionError("K inverse of zero")
        # z^{-1} = conj(z) * conj_r(z*conj(z)) / N where N = (z conj z)(conj_r(z conj z)) in Q
        zc = self.conj()
        w = self * zc            # in Q(r): c = d = 0
        assert w.c == 0 and w.d == 0
        wr = w.conj_r()
        n = w * wr               # rational
        assert n.b == 0 and n.c == 0 and n.d == 0 and n.a != 0
        num = zc * wr
        inv = 1 / n.a
        return K(num.a * inv, num.b * inv, num.c * inv, num.d * inv)

    def __truediv__(self, other):
        if not isinstance(other, K):
            o = _q(other)
            return K(self.a / o, self.b / o, self.c / o, self.d / o)
        return self * other.inverse()

    def __rtruediv__(self, other):
        return K(other) * self.inverse()

    def __pow__(self, n):
        n = int(n)
        if n < 0:
            return self.inverse() ** (-n)
        result = K(1)
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    # ---- printing
    def __repr__(self):
        parts = []
        for coef, name in ((self.a, ""), (self.b, "r"), (self.c, "xi"), (self.d, "r*xi")):
            if coef == 0:
                continue
            if name == "":
                parts.append(str(coef))
            elif coef == 1:
                parts.append(name)
            elif coef == -1:
                parts.append("-" + name)
            else:
                parts.append("(%s)*%s" % (coef, name))
        if not parts:
            return "0"
        return " + ".join(parts).replace("+ -", "- ")

    __str__ = __repr__

    def approx(self):
        """floating approximation, for display only (never used in a test)"""
        import math
        r = math.sqrt(13.0)
        return complex(float(self.a) + float(self.b) * r,
                       (float(self.c) + float(self.d) * r) * math.sqrt(3.0))


ZERO = K(0)
ONE = K(1)
R = K(0, 1, 0, 0)
XI = K(0, 0, 1, 0)
RXI = K(0, 0, 0, 1)


def kq(x):
    """rational -> K"""
    return K(mpq(x))


if __name__ == "__main__":
    # self-test
    assert R * R == K(13)
    assert XI * XI == K(-3)
    assert RXI == R * XI
    z = K(mpq(1, 2), mpq(-3, 7), mpq(5, 11), mpq(2, 3))
    assert z * z.inverse() == ONE
    assert (z / z) == ONE
    assert z ** 3 == z * z * z
    assert (z ** -2) * (z ** 2) == ONE
    assert z.conj().conj() == z
    print("a4_kfield self-test OK")
