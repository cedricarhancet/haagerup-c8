"""Axis A1 (case I, Griess subalgebra 2), CERTIFIED version of r5: test (b) for the finite-option types (1,13) and (1,1).  For u in the (1,13)
line (highest weight vector) and w in the (1,1) line, every component of u_(1) w' (w' in the module of w) must have a
U_1-type occurring in U_2; in (1,13) x (1,1) only (0,12) does.  We project with the Casimir operators of g' and g''
(exact rational vectors, bilinear in the line coefficients), and test whether the forbidden components vanish for each
of the 2 x 2 choices of lines.  Exact arithmetic (sympy over the quadratic fields of the lines)."""
import sys, os, json, pickle, itertools, time
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../cas-II')
sys.path.insert(0, GF); os.chdir(GF)
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, isotropic_lines as IL, closure_explore as CE
from flint import acb, arb, ctx
ctx.prec = 400
def to_acb(e):
    e = sp.sympify(e)
    if e.is_Rational: return acb(arb(int(e.p))/int(e.q))
    if e is sp.I: return acb(0, 1)
    if e.is_Add:
        r = acb(0)
        for a in e.args: r += to_acb(a)
        return r
    if e.is_Mul:
        r = acb(1)
        for a in e.args: r *= to_acb(a)
        return r
    if e.is_Pow and e.exp == sp.Rational(1, 2):
        b = to_acb(e.base)
        assert b.imag == 0 and b.real > 0, ('sqrt of non-positive', e.base)
        return acb(b.real.sqrt())
    if e.is_Pow and e.exp.is_Integer:
        return to_acb(e.base)**int(e.exp)
    if e.is_Pow and e.exp == sp.Rational(-1, 2):
        return 1/to_acb(e.base**sp.Rational(1, 2))
    raise ValueError('unsupported ' + str(e))
t0 = time.time()
CI = int(sys.argv[1]) if len(sys.argv) > 1 else 2
H113 = sys.argv[2] if len(sys.argv) > 2 else '0'
H11 = sys.argv[3] if len(sys.argv) > 3 else '2/5'
D, out = IL.admissible_spaces('I', CI); T = D['T']
S13 = sp.sqrt(13)
def rat(x): return sp.Rational(int(x.numerator), int(x.denominator))
def Ebasis(spc, P):
    """E columns -> list of (coefficient, rational dict) pairs: column k = sum_j coef_j * vec_j (coef in {1, sqrt13})"""
    res = []
    for k in range(spc['E'].cols):
        vp, vq = CE.split_vec(list(spc['E'][:, k]), P)
        res.append([(sp.Integer(1), vp), (S13, vq)])
    return res
def space(tau, h1):
    c = out[tau]['spaces']
    return c[0] if len(c) == 1 else [x for x in c if str(x['h1']) == h1][0]
s113 = space((1, 13), H113); s11 = space((1, 1), H11)
L113 = [list(l) for l in IL.lines_in(s113)]; L11 = [list(l) for l in IL.lines_in(s11)]
E113 = Ebasis(s113, out[(1, 13)]['P']); E11 = Ebasis(s11, out[(1, 1)]['P'])
def M(x, n, v): return V.clean(V.mode_vec(x, n, v))
def casimir(which, v):
    e, f, H = (T['e1'], T['f1'], T['H1']) if which == 1 else (T['e2'], T['f2'], T['H2'])
    Hv = M(H, 0, v); HHv = M(H, 0, Hv); fev = M(f, 0, M(e, 0, v))
    out_ = {}
    V.vadd(out_, HHv, mpq(1, 2)); V.vadd(out_, Hv, mpq(1)); V.vadd(out_, fev, mpq(2))
    return V.clean(out_)
def proj_forbidden(z):
    """(C1 - c1(0)) z and (C2 - c2(12)) z must both vanish if z has only type (0,12); c(a) = a(a+2)/2"""
    a = casimir(1, z)                                   # c1(0) = 0
    b = casimir(2, z); V.vadd(b, z, mpq(-12*14, 2))     # c2(12) = 84
    return V.clean(a), V.clean(b)
mods11 = [[(c, CE.lowerings(T, (1, 1), v)) for c, v in col] for col in E11]
nw = len(mods11[0][0][1])
R = {}
for ia, cola in enumerate(E113):
    for ja, (ca, a) in enumerate(cola):
        if not a: continue
        for ib, colb in enumerate(mods11):
            for jb, (cb, lows) in enumerate(colb):
                for w in range(nw):
                    if not lows[w]: continue
                    R[(ia, ja, ib, jb, w)] = (ca*cb, proj_forbidden(M(a, 1, lows[w])))
print('bilinear pieces computed (%.0fs)' % (time.time() - t0), flush=True)
res = {}; CERT = {}
for i, l1 in enumerate(L113):
    for j, l2 in enumerate(L11):
        fails = False
        for w in range(nw):
            for part in (0, 1):
                acc = {}
                for (ia, ja, ib, jb, ww), (c, pr) in R.items():
                    if ww != w: continue
                    coef = sp.expand(l1[ia]*l2[ib]*c)
                    if coef == 0: continue
                    for k, x in pr[part].items(): acc[k] = acc.get(k, 0) + coef*rat(x)
                # certified nonvanishing: an acb ball excluding 0 for some coordinate
                cert = None
                for k, x in acc.items():
                    x = sp.expand(x)
                    if x == 0: continue
                    b = to_acb(x)
                    if not b.contains(0):
                        cert = (str(k) + ' [part ' + ('C1 z' if part == 0 else '(C2-84) z') + ', w index %d]' % w, str(b)); break
                if cert:
                    fails = True; print('   certified nonzero forbidden component at', cert[0], '=', cert[1][:60]); CERT[(i, j)] = cert; break
            if fails: break
        res['(1,13) line %d, (1,1) line %d' % (i, j)] = ('EXCLUDED (certified): ' + CERT[(i, j)][0] + ' -> ' + CERT[(i, j)][1]) if fails else 'not excluded'
        print(i, j, res['(1,13) line %d, (1,1) line %d' % (i, j)], flush=True)
json.dump({'lines_113': [[str(c) for c in l] for l in L113], 'lines_11': [[str(c) for c in l] for l in L11],
           'test_b_casimir': res, 'time': time.time() - t0}, open(HERE + '/r6_certify_sub%d.json' % CI, 'w'), indent=1)
