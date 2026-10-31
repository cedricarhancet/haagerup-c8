"""Axis B1, explicit form: the element z = (1, 1, -1) of Z is monomial on V'_1, i.e. the lift of a lattice isometry g
of E6 + A2 with signs: e^lam -> c(lam) e^{g lam}, h(-n) -> (g h)(-n).  We
 (1) read g and c(alpha_i) from the matrix, extend c by c(lam + mu) = c(lam) c(mu) eps(lam, mu) eps(g lam, g mu)
     (eta-formula), and check it against all 78 roots (so the formula is the unique lift, a VOA automorphism);
 (2) check z(a_(n) b) = z(a)_(n) z(b) on all pairs of weight-one basis states, n = -1, 0, 1 (sanity);
 (3) test EXACTLY over K whether z maps U_2(9, choice 0) onto U_2(10, choice 0 or 1) = sigma(U_2(9, choice)),
     sigma: sqrt13 -> -sqrt13."""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath(HERE + '/../generation-forte'))
os.chdir(HERE + '/../generation-forte')
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V, exact_engine as E, Kfield as K
BITS = tuple(int(b) for b in (sys.argv[1] if len(sys.argv) > 1 else '110'))
Zm = [[sp.sympify(x) for x in r] for r in pickle.load(open(HERE + '/j3_z_%d%d%d.pkl' % BITS, 'rb'))]
B1 = V.basis(1); n = len(B1); idx = {k: i for i, k in enumerate(B1)}
cart = [i for i, k in enumerate(B1) if k[0] == V.ZERO]
# Cartan block: column for h_i = alpha_i(-1)
hkey = {k[1][0][1]: i for i, k in enumerate(B1) if k[0] == V.ZERO}
g = [[0]*8 for _ in range(8)]       # g[j][i] = coefficient of alpha_j in g(alpha_i)
for i in range(8):
    col = hkey[i]
    for j in range(8):
        x = Zm[hkey[j]][col]
        assert x == int(x); g[j][i] = int(x)
def gv(lam): return tuple(sum(g[j][i]*lam[i] for i in range(8)) for j in range(8))
c_root = {}
for jcol, k in enumerate(B1):
    if k[0] == V.ZERO: continue
    rows = [i for i in range(n) if Zm[i][jcol] != 0]; assert len(rows) == 1
    tgt = B1[rows[0]][0]; assert tgt == gv(k[0]), (k, tgt, gv(k[0]))
    c_root[k[0]] = Zm[rows[0]][jcol]
simple = [tuple(1 if t == i else 0 for t in range(8)) for i in range(8)]
cs = [c_root[a] for a in simple]
delta = [[V.eps(simple[i], simple[j])*V.eps(gv(simple[i]), gv(simple[j])) for j in range(8)] for i in range(8)]
def c(lam):
    v = sp.Integer(1)
    for i in range(8):
        v *= cs[i]**lam[i]
        if delta[i][i] == -1 and (lam[i]*(lam[i]-1)//2) % 2: v = -v
        for j in range(i+1, 8):
            if delta[i][j] == -1 and (lam[i]*lam[j]) % 2: v = -v
    return v
bad = [a for a in c_root if sp.simplify(c(a) - c_root[a]) != 0]
print('isometry g:', g); print('det g =', sp.Matrix(g).det(), '; g preserves the form:', all(V.ip(gv(a), gv(b)) == V.ip(a, b) for a in simple for b in simple))
print('c on simple roots:', cs, '; c-formula agrees on all roots:', not bad)
def zmono(m):
    """image of an oscillator monomial: product over (n, i) of sum_j g[j][i] alpha_j(-n)"""
    terms = {(): sp.Integer(1)}
    for (nn, i) in m:
        new = {}
        for mono, x in terms.items():
            for j in range(8):
                if g[j][i]:
                    mm = tuple(sorted(mono + ((nn, j),)))
                    new[mm] = new.get(mm, 0) + x*g[j][i]
        terms = new
    return terms
def zvec(v, ring='Q'):
    out = {}
    for (lam, m), x in v.items():
        cl = c(lam); gl = gv(lam)
        for mm, y in zmono(m).items():
            k = (gl, mm)
            if ring == 'Q':
                out[k] = out.get(k, 0) + x*mpq(int(sp.fraction(cl*y)[0]), int(sp.fraction(cl*y)[1]))
            else:
                s = mpq(int(sp.fraction(cl*y)[0]), int(sp.fraction(cl*y)[1]))
                out[k] = K.add(out[k], K.scal(s, x)) if k in out else K.scal(s, x)
    if ring == 'Q': return V.clean(out)
    return {k: a for k, a in out.items() if not K.iszero(a)}
# (2) sanity: products of weight-one states
fails = 0; cnt = 0
import random; random.seed(1)
pairs = random.sample([(a, b) for a in B1 for b in B1], 600)
for a, b in pairs:
    for nn in (-1, 0, 1):
        lhs = zvec(V.clean(V.mode_vec({a: mpq(1)}, nn, {b: mpq(1)})))
        rhs = V.clean(V.mode_vec(zvec({a: mpq(1)}), nn, zvec({b: mpq(1)})))
        cnt += 1
        if V.clean({k: lhs.get(k, 0) - rhs.get(k, 0) for k in set(lhs) | set(rhs)}): fails += 1
print('product checks z(a_n b) = z(a)_n z(b):', cnt, 'failures:', fails)
# U_1 fixed?
T = E.setup(0)['T']
print('z fixes U_1:', all(zvec(T[s]) == V.clean(T[s]) for s in ('e1', 'f1', 'H1', 'e2', 'f2', 'H2')))
# (3) exact membership
def sigma(v): return {k: (a[0], -a[1], a[2], -a[3]) for k, a in v.items()}
res = {}
for src in (0, 1):
    S = E.setup(src); B2 = E.Blocks(S['T'], 2)
    imgs = [zvec(v, ring='K') for v in S['U2span']]
    for tgt in (0, 1):
        St = E.setup(tgt)
        for twist in (0, 1):
            span = [sigma(v) if twist else v for v in St['U2span']]
            piv = {}
            for v in span: piv.setdefault(B2.wt(v), []).append(v)
            basis = {mu: [vs[i] for i in E.select_pivots(vs, B2.pos[mu], len(B2.blk[mu]))] for mu, vs in piv.items()}
            ok_all = True
            bymu = {}
            for w in imgs: bymu.setdefault(B2.wt(w), []).append(w)
            for mu, ws in bymu.items():
                ok, badv = E.solve_block(basis.get(mu, []), ws, B2.pos[mu], len(B2.blk[mu]))
                if not ok or badv: ok_all = False; break
            name = 'z(U_2(9,%d)) in U_2(%d,%d)' % (src, 10 if twist else 9, tgt)
            res[name] = ok_all; print(name, ok_all, flush=True)
json.dump({'bits': BITS, 'g': g, 'c_simple': [str(x) for x in cs], 'product_checks': [cnt, fails], 'membership': res},
          open(HERE + '/j5_explicit_%d%d%d.json' % BITS, 'w'), indent=1)
