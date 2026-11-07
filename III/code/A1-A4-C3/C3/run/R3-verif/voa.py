# Minimal lattice VOA V_L (L = E8+E8 in R^16, orthonormal Heisenberg basis h_0..h_15).
# State: dict {(mono, lat): Fraction}; mono = sorted tuple of (n, i) meaning h_i(-n), n>=1;
# lat = doubled lattice vector (16-tuple of ints).
# Y(e^b, z) = E^-(-b,z) E^+(-b,z) e_b z^{b(0)},  e_b e^g = eps(b,g) e^{b+g}.
# Composite modes by the iterate (Borcherds) formula.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import comb
from lat import ip, eps, vadd, vneg, ZERO

def sadd(acc, st, c=1):
    for k, v in st.items():
        nv = acc.get(k, 0) + c*v
        if nv:
            acc[k] = nv
        else:
            acc.pop(k, None)
    return acc

def clean(st):
    return {k: v for k, v in st.items() if v != 0}

def mono_mul(m1, m2):
    return tuple(sorted(m1 + m2))

def mono_wt(m):
    return sum(n for n, i in m)

def state_wt(key):
    m, lat = key
    return mono_wt(m) + Fr(sum(x*x for x in lat), 8)

# ---- Heisenberg ----
def heis_mode(i, n, st):
    """h_i(n) acting on state"""
    out = {}
    for (m, lat), c in st.items():
        if n < 0:
            k = (mono_mul(m, ((-n, i),)), lat)
            out[k] = out.get(k, 0) + c
        elif n == 0:
            if lat[i]:
                k = (m, lat)
                out[k] = out.get(k, 0) + c*Fr(lat[i], 2)
        else:
            cnt = m.count((n, i))
            if cnt:
                lst = list(m); lst.remove((n, i))
                k = (tuple(lst), lat)
                out[k] = out.get(k, 0) + c*n*cnt
    return clean(out)

def vec_mode(hv, n, st):
    """(sum_i hv_i h_i)(n), hv a dict i->coeff"""
    out = {}
    for i, a in hv.items():
        sadd(out, heis_mode(i, n, st), a)
    return out

# ---- polynomials in creation operators: dict mono->coeff ----
def annih(beta, m, poly):
    """beta(m), m>0, acting on Fock polynomial (dict mono->coeff); beta doubled vector"""
    out = {}
    for mono, c in poly.items():
        seen = set()
        for (n, i) in mono:
            if n != m or beta[i] == 0 or (n, i) in seen:
                continue
            seen.add((n, i))
            cnt = mono.count((n, i))
            lst = list(mono); lst.remove((n, i))
            k = tuple(lst)
            out[k] = out.get(k, 0) + c*m*cnt*Fr(beta[i], 2)
    return clean(out)

def Eplus(beta, poly):
    """E^+(-beta,z) poly = sum_e z^{-e} f_e ; returns dict e -> poly"""
    # exp(A), A = - sum_{m>0} beta(m) z^{-m}/m
    res = {0: dict(poly)}
    term = {0: dict(poly)}
    k = 0
    while term:
        k += 1
        new = {}
        for e, p in term.items():
            for m in range(1, mono_maxn(p) + 1):
                q = annih(beta, m, p)
                if q:
                    d = new.setdefault(e + m, {})
                    sadd(d, q, Fr(-1, m*k))
        term = {e: clean(p) for e, p in new.items() if clean(p)}
        for e, p in term.items():
            sadd(res.setdefault(e, {}), p)
    return {e: p for e, p in res.items() if p}

def mono_maxn(p):
    mx = 0
    for mono in p:
        for n, i in mono:
            mx = max(mx, n)
    return mx

_Scache = {}
def Sd(beta, d):
    """coefficient of z^d in E^-(-beta,z) = exp(sum_{m>0} beta(-m) z^m/m), as dict mono->coeff"""
    key = (beta, d)
    if key in _Scache:
        return _Scache[key]
    if d == 0:
        r = {(): Fr(1)}
    else:
        # recursion: d S_d = sum_{m=1}^d beta(-m) S_{d-m}
        r = {}
        for m in range(1, d+1):
            prev = Sd(beta, d-m)
            for mono, c in prev.items():
                for i in range(16):
                    if beta[i]:
                        k = mono_mul(mono, ((m, i),))
                        r[k] = r.get(k, 0) + c*Fr(beta[i], 2)/d
        r = clean(r)
    _Scache[key] = r
    return r

def e_mode(beta, n, st):
    """(e^beta)_(n) acting on state"""
    out = {}
    for (m, g), c in st.items():
        p = ip(beta, g)
        s = eps(beta, g)
        newlat = vadd(beta, g)
        ep = Eplus(beta, {m: Fr(1)})
        for e, poly in ep.items():
            d = -n - 1 - p + e
            if d < 0:
                continue
            S = Sd(beta, d)
            for m1, c1 in poly.items():
                for m2, c2 in S.items():
                    k = (mono_mul(m1, m2), newlat)
                    out[k] = out.get(k, 0) + s*c*c1*c2
    return clean(out)

# ---- general modes via iterate formula ----
def maxwt(st):
    return max(state_wt(k) for k in st) if st else 0

def mode(akey, n, st):
    """a_(n) st for basis state akey=(mono, lat)"""
    m, lat = akey
    if not m:
        return e_mode(lat, n, st)
    # peel the last creation operator: a = h_i(-q) b
    q, i = m[-1]
    b = (m[:-1], lat)
    wb = state_wt(b)
    wv = maxwt(st)
    out = {}
    # (h_(-q) b)_(n) = sum_j binom(q+j-1,j) [ h(-q-j) b_(n+j) - (-1)^q b_(n-q-j) h(j) ]
    j = 0
    while True:
        cf = comb(q+j-1, j)
        active = False
        # first term: b_(n+j) v nonzero needs output weight wb+wv-(n+j)-1 >= 0
        if wb + wv - (n + j) - 1 >= 0:
            active = True
            t = mode(b, n+j, st)
            if t:
                sadd(out, heis_mode(i, -q-j, t), cf)
        # second term: h(j) v nonzero needs j <= wv (j=0 allowed)
        if j <= wv:
            active = True
            t = heis_mode(i, j, st)
            if t:
                sadd(out, mode(b, n-q-j, t), -cf*(-1)**q)
        if not active:
            break
        j += 1
    return out

def mode_vec(a, n, st):
    out = {}
    for k, c in a.items():
        sadd(out, mode(k, n, st), c)
    return out

VAC = {((), ZERO): Fr(1)}
def h_state(i):
    return {(((1, i),), ZERO): Fr(1)}
def e_state(beta):
    return {((), beta): Fr(1)}
