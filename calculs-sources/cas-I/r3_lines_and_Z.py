"""Axis A1 (case I): does an element of the centralizer Z of U_1 (j6, monomial lifts) fixing the Griess subalgebra 5
exchange the two isotropic lines of type (1,13)?  We apply the lift of each z (lattice isometry g + signs c) to v
(Riesen's generator, transported, r1) and locate z(v) in E."""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../cas-II')
sys.path.insert(0, GF); os.chdir(GF)
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, isotropic_lines as IL
FULL = {n: {k: sp.sympify(c) for k, c in d.items()} for n, d in pickle.load(open(HERE + '/r1_riesen_transported.pkl', 'rb')).items()}
B1 = V.basis(1); n = len(B1)
hkey = {k[1][0][1]: i for i, k in enumerate(B1) if k[0] == V.ZERO}
simple = [tuple(1 if t == i else 0 for t in range(8)) for i in range(8)]
def lift(Zm):
    g = [[0]*8 for _ in range(8)]
    for i in range(8):
        for j in range(8):
            x = Zm[hkey[j]][hkey[i]]; assert x == int(x); g[j][i] = int(x)
    gv = lambda lam: tuple(sum(g[j][i]*lam[i] for i in range(8)) for j in range(8))
    croot = {}
    for jcol, k in enumerate(B1):
        if k[0] == V.ZERO: continue
        rows = [i for i in range(n) if Zm[i][jcol] != 0]; assert len(rows) == 1
        assert B1[rows[0]][0] == gv(k[0]); croot[k[0]] = Zm[rows[0]][jcol]
    cs = [croot[a] for a in simple]
    delta = [[V.eps(simple[i], simple[j])*V.eps(gv(simple[i]), gv(simple[j])) for j in range(8)] for i in range(8)]
    def c(lam):
        v = sp.Integer(1)
        for i in range(8):
            v *= cs[i]**lam[i]
            if delta[i][i] == -1 and (lam[i]*(lam[i]-1)//2) % 2: v = -v
            for j in range(i+1, 8):
                if delta[i][j] == -1 and (lam[i]*lam[j]) % 2: v = -v
        return v
    assert all(sp.simplify(c(a) - croot[a]) == 0 for a in croot)
    def zvec(vv):
        out = {}
        for (lam, m), x in vv.items():
            terms = {(): sp.Integer(1)}
            for (nn, i) in m:
                new = {}
                for mono, y in terms.items():
                    for j in range(8):
                        if g[j][i]:
                            mm = tuple(sorted(mono + ((nn, j),))); new[mm] = new.get(mm, 0) + y*g[j][i]
                terms = new
            for mm, y in terms.items():
                k = (gv(lam), mm); out[k] = sp.expand(out.get(k, 0) + c(lam)*y*x)
        return {k: v for k, v in out.items() if v != 0}
    return g, cs, zvec
d = pickle.load(open(HERE + '/r2_type113.pkl', 'rb'))
P = d['P']; Em = sp.Matrix([[sp.sympify(x) for x in r] for r in d['E']]); lines = [[sp.sympify(c) for c in l] for l in d['lines']]
def rat(x): return sp.Rational(int(x.numerator), int(x.denominator))
def coordsP(target):
    keys = sorted(set().union(*[set(p) for p in P]) | set(target))
    M = sp.Matrix([[rat(p[k]) if k in p else 0 for p in P] for k in keys]); b = sp.Matrix([target.get(k, 0) for k in keys])
    rows = next(r for r in itertools.combinations(range(len(keys)), len(P)) if M[list(r), :].det() != 0)
    x = (M[list(rows), :].inv()*b[list(rows), :]).applyfunc(sp.expand)
    assert (M*x - b).applyfunc(lambda z: sp.simplify(sp.expand(z))) == sp.zeros(len(keys), 1), 'not in Pi'
    return x
def which_line(target):
    x = coordsP(target)
    rows = next(r for r in itertools.combinations(range(Em.rows), 2) if sp.simplify(Em[list(r), :].det()) != 0)
    y = (Em[list(rows), :].inv()*x[list(rows), :]).applyfunc(lambda z: sp.expand(sp.radsimp(z)))
    assert (Em*y - x).applyfunc(lambda z: sp.simplify(sp.expand(z))) == sp.zeros(Em.rows, 1), 'not in E'
    return [j for j, l in enumerate(lines) if sp.simplify(sp.expand(y[0]*l[1] - y[1]*l[0])) == 0]
res = {}
v = FULL['v']
print('v on line', which_line(v))
for bits in ('111', '110', '001', '000', '011', '010', '101', '100'):
    Zm = [[sp.sympify(x) for x in r] for r in pickle.load(open(HERE + '/../jumeau/j6_zI_%s.pkl' % bits, 'rb'))]
    g, cs, zvec = lift(Zm)
    zv = zvec(v)
    try:
        wl = which_line(zv)
    except AssertionError as ex:
        wl = 'not in this E (%s)' % ex
    eps = tuple(2*int(b) - 1 for b in bits)
    res[str(eps)] = {'g': g, 'c': [str(x) for x in cs], 'z(v) line': wl}
    print(eps, 'c =', cs, 'z(v) on line', wl, flush=True)
json.dump(res, open(HERE + '/r3_lines_and_Z.json', 'w'), indent=1)
