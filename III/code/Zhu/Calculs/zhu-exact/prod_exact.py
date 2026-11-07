"""Exact products for the exact Zhu relations of case II (independent rewrite of the product stage).
For a <= b in the lowering basis of U_2 (120 vectors z_b, exact over K):
  p2 = a(2)b in V'_1, p1 = a(1)b in V'_2,
  B(a,b) = a(0)b - T(a(1)b)/2 + T^2(a(2)b)/4                      in V'_3,
  N(a,b) = a(-1)b - T(a(0)b)/2 + T^2(a(1)b)/4 - T^3(a(2)b)/12      in V'_4 (only if wt(a)+wt(b) is dominant).
Usage: python3 prod_exact.py k m   (handles a with a % m == k). Output: prodx/a%03d.pkl (exact K-vectors)."""
import sys, os, time, pickle, gc
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu')
from zhu_common import *
k, m = int(sys.argv[1]), int(sys.argv[2])
OUT = '/home/claude/hj/work/Calculs/zhu-exact/prodx'; os.makedirs(OUT, exist_ok=True)
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']
d = pickle.load(open('/home/claude/hj/work/Calculs/zhu/U2_lowering_basis.pkl', 'rb'))
HV = d['hv']; LB = [(a, i, j) for a, i, j, _ in d['basis']]; TYPE = [w for _, w in HV]; Z = [b[3] for b in d['basis']]
def wt(b): a, i, j = LB[b]; r, l = TYPE[a]; return (r - 2 * i, l - 2 * j)
fnT = OUT + '/Tops.pkl'
if os.path.exists(fnT): TO = pickle.load(open(fnT, 'rb'))
else:
    TO = {w: E.build_op([(K.ONE, V.OMEGA)], 0, w) for w in (1, 2, 3)}
    if k == 0: pickle.dump(TO, open(fnT + '.tmp', 'wb')); os.replace(fnT + '.tmp', fnT)
log('T operators ready')
def Tapp(v, w, n=1):
    for _ in range(n): v = E.apply_op(TO[w], v); w += 1
    return v
def lc(*pairs):
    out = {}
    for c, v in pairs:
        if v: E.kv_add(out, v, K.scal(mpq(c), K.ONE))
    return E.kv_clean(out)
for a in range(120):
    if a % m != k: continue
    fn = OUT + '/a%03d.pkl' % a
    if os.path.exists(fn): continue
    need_m1 = any(wt(a)[0] + wt(b)[0] >= 0 and wt(a)[1] + wt(b)[1] >= 0 for b in range(a, 120))
    ops = {n: E.build_op(E.kv_terms(Z[a]), n, 2) for n in ((2, 1, 0, -1) if need_m1 else (2, 1, 0))}
    res = {}
    for b in range(a, 120):
        p2 = E.apply_op(ops[2], Z[b]); p1 = E.apply_op(ops[1], Z[b]); p0 = E.apply_op(ops[0], Z[b])
        Tp2 = Tapp(p2, 1); T2p2 = Tapp(Tp2, 2); Tp1 = Tapp(p1, 2)
        B = lc((1, p0), (mpq(-1, 2), Tp1), (mpq(1, 4), T2p2))
        N = None
        if wt(a)[0] + wt(b)[0] >= 0 and wt(a)[1] + wt(b)[1] >= 0:
            pm1 = E.apply_op(ops[-1], Z[b])
            N = lc((1, pm1), (mpq(-1, 2), Tapp(p0, 3)), (mpq(1, 4), Tapp(Tp1, 3)), (mpq(-1, 12), Tapp(T2p2, 3)))
        res[b] = {'p2': p2, 'p1': p1, 'Tp2': Tp2, 'B': B, 'N': N}
    del ops; V._MODE.clear(); V._EMODE.clear(); gc.collect()
    pickle.dump(res, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
    del res; gc.collect()
    log('a =', a, 'done')
log('finished')
