"""Exact eigenvalues of the central class Z = e[omega]e on the twelve blocks of Lambda (x) K.
Uses the exact presentation upper_exact.pkl (z4x.py): with the irreducible paths as K-basis (exact confluence),
computes the element P(Z) = prod_{h in H} (Z - h) applied to the unit 1 = sum_v e_v exactly and checks P(Z) = 0,
where H = {0, 1/3, 2/3, 1, 5/13, 6/13, 7/13, 8/13, 11/13, 15/13}.  Since Lambda (x) K is semisimple and Z central,
the eigenvalue of Z on each block is then an element of H, determined by its reduction modulo p."""
import pickle, json, time, sys
from fractions import Fraction as Fr
from functools import lru_cache
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu-exact')
import kfast as KF, kexact as KX
from gmpy2 import mpq
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
HERE = '/home/claude/hj/work/Calculs/zhu-exact'
U = pickle.load(open(HERE + '/upper_exact.pkl', 'rb'))
D = pickle.load(open(HERE + '/rows_exact.pkl', 'rb'))
AR = D['ARROWS']; CONS = D['CONS']; VERT = D['VERT']
rules, alllin, killed, active, verts = U['rules'], U['alllin'], set(U['killed']), U['active'], U['verts']
def src(a): return AR[a]['source']
def tgt(a): return AR[a]['target']
def add(out, row, c=KF.ONE):
    for p, v in row.items():
        z = KF.add(out.get(p, KF.ZERO), KF.mul(c, v))
        if KF.iszero(z): out.pop(p, None)
        else: out[p] = z
LENS = sorted({len(q) for q in rules}, reverse=True)
@lru_cache(None)
def nf(p):
    for L in LENS:
        for i in range(len(p) - L + 1):
            q = p[i:i + L]
            if q in rules:
                o = {}
                for w, c in rules[q].items(): add(o, nf(p[:i] + w + p[i + L:]), c)
                return o
    return {p: KF.ONE}
@lru_cache(None)
def fullnf(p):
    if any(src(a) in killed or tgt(a) in killed for a in p): return {}
    for i, a in enumerate(p):
        if a in alllin:
            o = {}
            for w, c in alllin[a].items(): add(o, fullnf(p[:i] + w + p[i + 1:]), c)
            return o
    return nf(p)
def hA(t): r, l = t; return Fr(r * (r + 2), 20) + Fr(l * (l + 2), 60)
def kq(x): x = Fr(x); return (mpq(x.numerator, x.denominator), mpq(0), mpq(0), mpq(0))
# elements: dict key -> K, key = ('v', v) for vertices or a path tuple (normal form)
eloop = {A_['src_c']: x for x, A_ in enumerate(AR) if A_['z'] == 'P:e' and A_['src_c'] == A_['tgt_c']}
Zel = {}
for v in verts:
    c0 = [c for c in range(len(CONS)) if VERT[c] == v][0]
    Zel[('v', v)] = kq(hA(CONS[c0]))
for c, x in eloop.items():
    if VERT[c] in verts:
        for w, cc in fullnf((x,)).items():
            k = ('v', tgt(x)) if w == () else w
            Zel[k] = KF.add(Zel.get(k, KF.ZERO), KF.mul(kq(2), cc))
Zel = {k: c for k, c in Zel.items() if not KF.iszero(c)}
log('Z has', len(Zel), 'terms')
def ends(k): return (k[1], k[1]) if k[0] == 'v' else (tgt(k[0]), src(k[-1]))
@lru_cache(None)
def prod_keys(k1, k2):
    e1, e2 = ends(k1), ends(k2)
    if e1[1] != e2[0]: return {}
    if k1[0] == 'v': return {k2: KF.ONE}
    if k2[0] == 'v': return {k1: KF.ONE}
    out = {}
    for w, c in nf(k1 + k2).items():
        out[('v', tgt(k1[0])) if w == () else w] = c
    return out
def mul(x, y):
    out = {}
    for k1, a in x.items():
        for k2, b in y.items():
            ab = KF.mul(a, b)
            for k, c in prod_keys(k1, k2).items():
                z = KF.add(out.get(k, KF.ZERO), KF.mul(ab, c))
                if KF.iszero(z): out.pop(k, None)
                else: out[k] = z
    return out
H = [Fr(0), Fr(1, 3), Fr(2, 3), Fr(1), Fr(5, 13), Fr(6, 13), Fr(7, 13), Fr(8, 13), Fr(11, 13), Fr(15, 13)]
one = {('v', v): KF.ONE for v in verts}
x = dict(one)
for h in H:
    zx = mul(Zel, x)
    for k, c in x.items():
        z = KF.add(zx.get(k, KF.ZERO), KF.mul(kq(-h), c))
        if KF.iszero(z): zx.pop(k, None)
        else: zx[k] = z
    x = zx
    log('after factor', h, 'terms', len(x))
ok = (len(x) == 0)
# also check that no proper sub-product vanishes (each eigenvalue occurs): P_h(Z) != 0 for each h
missing = []
for h0 in H:
    y = dict(one)
    for h in H:
        if h == h0: continue
        zy = mul(Zel, y)
        for k, c in y.items():
            z = KF.add(zy.get(k, KF.ZERO), KF.mul(kq(-h), c))
            if KF.iszero(z): zy.pop(k, None)
            else: zy[k] = z
        y = zy
    if not y: missing.append(str(h0))
log('prod (Z - h) = 0 exactly:', ok, '; eigenvalues absent:', missing)
json.dump({'H': [str(h) for h in H], 'product_vanishes': ok, 'absent': missing, 'Z_terms': len(Zel)}, open(HERE + '/exactZ.json', 'w'), indent=1)
log('done')
