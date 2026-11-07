"""R3a step 3c: the full space S_2 x 1 (exact over Q(sqrt(-6))): primaries (+ lowering descendants) and the 27 affine
descendants x(-2)1, x(-1)y(-1)1 (x, y in S_1 x 1).  Checks: dimension 120 and sl2^2-decomposition = "U2" of
t12_table.json; commutant: (1 x x)_(n) v = 0 (n = 0, 1, 2), (omega_{1 x C})_(n) v = 0 (n = 1, 2, 3) for all v;
closure: a_(1) b in S_2 x 1, a_(2) b in S_1 x 1, a_(3) b in C 1 for a = primary hw vectors (and omega's), b in the basis.
Usage: python3 -I r3a_S2_full.py CI [full].  Output: S2full_CI.pkl, S2full_CI.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time, itertools
from collections import Counter
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import e8e8_voa as V
import common as CM
import qs6 as Q
t0 = time.time()
CI = int(sys.argv[1]) if len(sys.argv) > 1 else 0
FULL = len(sys.argv) > 2 and sys.argv[2] == 'full'
DG = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
OC = pickle.load(open(os.path.join(CM.HERE, 'omegaC1.pkl'), 'rb'))
S2 = pickle.load(open(os.path.join(CM.HERE, 'S2_%d.pkl' % CI), 'rb'))
u, wC1 = OC['cands'][CI]
om1C = Q.vclean(Q.vadd(Q.vreal(DG['omCC']), wC1, (mpq(-1), mpq(0))))
Ep, Hp, Fp = CM.SL['g1p']; Epp, Hpp, Fpp = CM.SL['g1pp']
S1 = [Ep, Hp, Fp, Epp, Hpp, Fpp]
from span import Span, smul
# ------------------------------------------------ build S_2 x 1
prim_hw = []
vecs = []
for tau, vv in sorted(S2['prims'].items()):
    for p in vv:
        prim_hw.append((tau, p))
        a, b = tau
        x = p
        for i in range(a + 1):
            y = x
            for j in range(b + 1):
                vecs.append(y); y = Q.mode_rq(Fpp, 0, y)
            x = Q.mode_rq(Fp, 0, x)
aff = []
for x in S1: aff.append(Q.vreal(V.clean(V.mode_vec(x, -2, V.one()))))
for x, y in itertools.combinations_with_replacement(S1, 2): aff.append(Q.vreal(V.clean(V.mode_vec(x, -1, y))))
allv = vecs + aff
SP = Span(allv)
print('dim S_2 x 1 = %d (expected 120) ; %d primary multiplets vectors + %d affine (%.0fs)' % (SP.rank, len(vecs), len(aff), time.time() - t0), flush=True)
# sl2^2 decomposition
def wt12(key): return CM.WW.weight_of(key)[:2]
ch = Counter()
for v in SP.B:
    ws = set(wt12(k) for part in v for k in part)
    assert len(ws) == 1; ch[ws.pop()] += 1
from weights_W import decomp
dec = decomp(ch)
U2 = {tuple(map(int, k.split(','))): m for k, m in json.load(open(os.path.join(CM.HERE, '..', 'data', 't12_table.json')))['U2'].items()}
print('sl2^2 decomposition:', sorted(dec.items()), '; equal to U2:', dict(dec) == U2, flush=True)
# commutant
bad = 0
for v in SP.B:
    for x in (CM.SL['g2p'] + CM.SL['g2pp']):
        for n in (0, 1, 2):
            if not Q.iszero_v(Q.mode_rq(x, n, v)): bad += 1
    for n in (1, 2, 3):
        if not Q.iszero_v(Q.mode(om1C, n, v)): bad += 1
print('commutant conditions (1 x S_1 modes 0,1,2 ; omega_{1 x C} modes 1,2,3) on the 120 basis vectors: failures', bad, '(%.0fs)' % (time.time() - t0), flush=True)
# closure
S1sp = Span([Q.vreal(x) for x in S1])
gens = [p for tau, p in prim_hw]
targets = SP.B if FULL else SP.B[::7]
fails = Counter(); tests = 0
for a in gens:
    for b in targets:
        r1 = Q.mode(a, 1, b); tests += 1
        if SP.coords(r1) is None: fails[1] += 1
        r2 = Q.mode(a, 2, b)
        if not Q.iszero_v(r2) and S1sp.coords(r2) is None: fails[2] += 1
        r3 = Q.mode(a, 3, b)
        rest = Q.vclean(Q.vadd(Q.vzero(), r3)); rest[0].pop(V.VAC, None); rest[1].pop(V.VAC, None)
        if not Q.iszero_v(rest): fails[3] += 1
    print('  closure: generator done (%.0fs)' % (time.time() - t0), flush=True)
print('closure tests a_(n) b (a: %d primary hw vectors, b: %d basis vectors, n = 1, 2, 3): %d pairs, failures %s' %
      (len(gens), len(targets), tests, dict(fails)), flush=True)
pickle.dump({'B': SP.B, 'prim_hw': prim_hw}, open(os.path.join(CM.HERE, 'S2full_%d.pkl' % CI), 'wb'))
json.dump({'dim': SP.rank, 'U2_match': dict(dec) == U2, 'commutant_failures': bad, 'closure_pairs': tests,
           'closure_failures': dict(fails), 'full': FULL}, open(os.path.join(CM.HERE, 'S2full_%d.json' % CI), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
