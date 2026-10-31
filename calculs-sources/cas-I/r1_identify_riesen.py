"""Axis A0: explicit identification of Riesen's case (I) data with ours.  Riesen's lattice: A2 on coordinates 0,1, E6
(Bourbaki) on 2..7, cocycle eps_R(a,b) = (-1)^{sum a_i b_i + sum_{i>j} a_i Q_ij b_j}; field K = Q(r,q), r^2 = 13,
q^2 = -3 (same basis as ours).  We compose: coordinate isometry phi (R 0,1 -> our 6,7; R 2..7 -> our 0..5), an
element w of W(E6) x W(A2) sending (H_s, H_t) to our (H1, H2), the cocycle correction eta, and a torus character t
matching E_s, E_t with our e1, e2.  Then we transport Riesen's e, f, v (his generators U = <A, e, f, v>) and test:
 e, f in our C_2 of the Griess subalgebra 5 or 6; v in the admissible space of type (1,13); which isotropic line."""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../cas-II')
sys.path.insert(0, GF); os.chdir(GF)
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, cases as C
RD = json.load(open(HERE + '/../../Sources/cas-I/riesen/shared/definition/certificates/definition/defining_states.json'))['states']
QR = json.load(open(HERE + '/../../Sources/cas-I/riesen/shared/definition/data/definition.json'))['Q']
PERM = [6, 7, 0, 1, 2, 3, 4, 5]
def phi(v): out = [0]*8; [out.__setitem__(PERM[i], v[i]) for i in range(8)]; return tuple(out)
assert all(V.G[PERM[i]][PERM[j]] == QR[i][j] for i in range(8) for j in range(8)), 'phi not an isometry'
def epsR(a, b):
    e = sum(a[i]*b[i] for i in range(8)) + sum(a[i]*QR[i][j]*b[j] for i in range(8) for j in range(i))
    return -1 if e % 2 else 1
T = C.build('I')
def cart(state):   # Cartan element -> root coordinate vector (Heisenberg in root basis)
    v = [0]*8
    for s in state:
        assert not any(s['momentum']) and len(s['oscillators']) == 1 and s['oscillators'][0][0] == 1
        v[s['oscillators'][0][1]] += sp.Rational(s['coefficient'][0])
    return tuple(v)
Hs, Ht = cart(RD['H_s']), cart(RD['H_t'])
def ours_cart(d):
    v = [0]*8
    for (lam, m), x in d.items():
        assert lam == V.ZERO and len(m) == 1; v[m[0][1]] += sp.Rational(int(x.numerator), int(x.denominator))
    return tuple(v)
H1, H2 = ours_cart(T['H1']), ours_cart(T['H2'])
print('Riesen H_s, H_t:', Hs, Ht, ' ours H1, H2:', H1, H2)
def refl(v, i):
    c = sum(v[j]*V.G[j][i] for j in range(8)); w = list(v); w[i] -= c; return tuple(w)
start = (phi(Hs), phi(Ht)); target = (H1, H2)
seen = {start: []}; frontier = [start]; word = None
while frontier and word is None:
    nxt = []
    for pair in frontier:
        for i in range(8):
            q = (refl(pair[0], i), refl(pair[1], i))
            if q not in seen:
                seen[q] = seen[pair] + [i]; nxt.append(q)
                if q == target: word = seen[q]; break
        if word is not None: break
    frontier = nxt
print('orbit explored', len(seen), 'word (0-based, ours):', word)
def psi(v):
    v = phi(v)
    for i in word: v = refl(v, i)
    return v
BAS = [tuple(1 if j == i else 0 for j in range(8)) for i in range(8)]
RIJ = [[V.eps(psi(BAS[i]), psi(BAS[j]))*epsR(BAS[i], BAS[j]) for j in range(8)] for i in range(8)]
assert all(RIJ[i][j] == RIJ[j][i] for i in range(8) for j in range(8))
def eta(p):
    s = 1
    for i in range(8):
        if RIJ[i][i] == -1 and (p[i]*(p[i]-1)//2) % 2: s = -s
        for j in range(i+1, 8):
            if RIJ[i][j] == -1 and (p[i]*p[j]) % 2: s = -s
    return s
s13, s3 = sp.sqrt(13), sp.sqrt(-3)
def kval(c): return sp.Rational(c[0]) + sp.Rational(c[1])*s13 + sp.Rational(c[2])*s3 + sp.Rational(c[3])*s13*s3
def transport(name):
    out = {}
    for s in RD[name]:
        p = tuple(s['momentum']); osc = tuple(sorted((m, PERM[i]) for m, i in s['oscillators']))
        if osc:   # Heisenberg oscillators: alpha_i -> psi(alpha_i) expanded
            terms = {(): sp.Integer(1)}
            for (m, i) in s['oscillators']:
                img = psi(BAS[i]); new = {}
                for mono, x in terms.items():
                    for j in range(8):
                        if img[j]:
                            mm = tuple(sorted(mono + ((m, j),))); new[mm] = new.get(mm, 0) + x*img[j]
                terms = new
            for mm, y in terms.items():
                k = (psi(p), mm); out[k] = out.get(k, 0) + eta(p)*y*kval(s['coefficient'])
        else:
            k = (psi(p), ()); out[k] = out.get(k, 0) + eta(p)*kval(s['coefficient'])
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}
TR = {n: transport(n) for n in RD}
def ourd(d): return {k: sp.Rational(int(x.numerator), int(x.denominator)) for k, x in d.items()}
ours = {n: ourd(T[n]) for n in ('e1', 'f1', 'H1', 'e2', 'f2', 'H2')}
# torus: t(p) TR[E](p) = ours[e](p)
eqs = []
for a, b in (('E_s', 'e1'), ('E_t', 'e2')):
    assert set(TR[a]) == set(ours[b]), (a, sorted(TR[a]), sorted(ours[b]))
    for k in TR[a]: eqs.append((k[0], sp.nsimplify(ours[b][k]/TR[a][k])))
print('ratios:', [(k, str(r)) for k, r in eqs])
tt = sp.symbols('t0:8')
sol = sp.solve([sp.Mul(*[tt[i]**k[i] for i in range(8)]) - r for k, r in eqs], tt, dict=True)
print('torus solutions:', len(sol)); sol = sol[0]; print(sol)
def tval(p): return sp.Mul(*[sol.get(tt[i], tt[i])**p[i] for i in range(8)])
def full(name): return {k: sp.expand(sp.radsimp(tval(k[0])*c)) for k, c in TR[name].items()}
res = {'word': word, 'torus': {str(k): str(v) for k, v in sol.items()}}
for a, b in (('E_s', 'e1'), ('F_s', 'f1'), ('H_s', 'H1'), ('E_t', 'e2'), ('F_t', 'f2'), ('H_t', 'H2')):
    m = full(a); o = ours[b]
    res[a + ' -> ' + b] = set(m) == set(o) and all(sp.simplify(m[k] - o[k]) == 0 for k in m)
    print(a, '->', b, res[a + ' -> ' + b], flush=True)
FULL = {n: full(n) for n in ('e', 'f', 'v', 'omega_C', 'omega_A', 'omega_W', 'b_0', 'b_1', 'b_2', 'b_3')}
pickle.dump({n: {k: str(c) for k, c in d.items()} for n, d in FULL.items()}, open(HERE + '/r1_riesen_transported.pkl', 'wb'))
json.dump(res, open(HERE + '/r1_identify_riesen.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
