"""Exact Lie bracket extracted verbatim from the original find_signs.py."""
from lat import ip, vadd, eps
from fractions import Fraction as Fr

def br(x, y):
    """Lie bracket x_(0)y on W_1; elements dict: ('h', i) or ('e', alpha) -> coeff"""
    out = {}
    def addto(k, c):
        out[k] = out.get(k, 0) + c
        if out[k] == 0: del out[k]
    for kx, cx in x.items():
        for ky, cy in y.items():
            if kx[0] == 'h' and ky[0] == 'h':
                continue
            if kx[0] == 'h' and ky[0] == 'e':
                v = ky[1][kx[1]]  # doubled coord
                if v: addto(ky, Fr(v, 2)*cx*cy)
                continue
            if kx[0] == 'e' and ky[0] == 'h':
                v = kx[1][ky[1]]
                if v: addto(kx, -Fr(v, 2)*cx*cy)
                continue
            a, b = kx[1], ky[1]
            p = ip(a, b)
            if p == -1:
                addto(('e', vadd(a, b)), eps(a, b)*cx*cy)
            elif p == -2:
                e_ = eps(a, b)
                for i in range(16):
                    if a[i]: addto(('h', i), e_*Fr(a[i], 2)*cx*cy)
    return out
