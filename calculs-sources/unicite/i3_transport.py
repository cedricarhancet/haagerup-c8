"""Axis I: transport of Jiao's algebra into our lattice model.  psi = w o phi (phi: diagram isometry, w: Weyl element
of i2_weyl.json); e^p -> eta(p) e^{psi p}, h(-n) -> (psi h)(-n), with eta(p) = prod_{i<j} r_ij^{p_i p_j}
prod_i r_ii^{p_i(p_i-1)/2}, r(p,q) = eps_ours(psi p, psi q) / eps_Jiao(p, q) (bimultiplicative, symmetric).
Checks (exact): sl2 relations and levels of the transported triples; comparison with our currents."""
import sys, os, json, pickle
from fractions import Fraction as F
HERE = os.path.abspath('.'); sys.path.insert(0, HERE); sys.path.insert(0, os.path.abspath('../generation-forte'))
os.chdir('../generation-forte')
from gmpy2 import mpq
import lattice_voa as V, exact_engine as E, Kfield as K
import jiao_data as J
W = json.load(open(HERE + '/i2_weyl.json'))
PERM = [0, 2, 3, 4, 5, 1, 6, 7]
def phi(v): out = [0]*8; [out.__setitem__(PERM[i], v[i]) for i in range(8)]; return tuple(out)
def refl(v, i):
    c = sum(v[j]*V.G[j][i] for j in range(8)); w = list(v); w[i] -= c; return tuple(w)
def psi(v):
    v = phi(v)
    for i in W['word']: v = refl(v, i)
    return v
BAS = [tuple(1 if j == i else 0 for j in range(8)) for i in range(8)]
def r(p, q): return V.eps(psi(p), psi(q)) * J.epsJ(p, q)
RIJ = [[r(BAS[i], BAS[j]) for j in range(8)] for i in range(8)]
assert all(RIJ[i][j] == RIJ[j][i] for i in range(8) for j in range(8))
def eta(p):
    s = 1
    for i in range(8):
        if RIJ[i][i] == -1 and (p[i]*(p[i]-1)//2) % 2: s = -s
        for j in range(i+1, 8):
            if RIJ[i][j] == -1 and (p[i]*p[j]) % 2: s = -s
    return s
def tk(a): return tuple(mpq(x.numerator, x.denominator) for x in a)
def transport_e(d):
    return {(psi(p), ()): K.scal(mpq(eta(p)), tk(c)) for p, c in d.items()}
def heis(h):        # h in Jiao root coordinates -> sum_i (psi h)_i alpha_i(-1) 1 in our basis (Heisenberg in root basis)
    v = psi(h); return {(V.ZERO, ((1, i),)): (mpq(v[i]), mpq(0), mpq(0), mpq(0)) for i in range(8) if v[i]}
def kclean(v): return {k: a for k, a in v.items() if not K.iszero(a)}
TR = {'PE': transport_e(J.PE), 'PF': transport_e(J.PF), 'PH': heis(J.PH), 'BE': transport_e(J.BE), 'BF': transport_e(J.BF),
      'BH': heis(J.BH), 'U26': transport_e(J.U26), 'U17': transport_e(J.U17)}
def prod(u, n, v): return E.kprod(E.kv_terms(u), n, v)
def sub(a, b): out = dict(a); E.kv_add(out, b, K.scal(mpq(-1), K.ONE)); return kclean(out)
res = {}
for lab, (e, f, h, lev) in {'level3': ('PE', 'PF', 'PH', 3), 'level13': ('BE', 'BF', 'BH', 13)}.items():
    ef = prod(TR[e], 0, TR[f]); he = prod(TR[h], 0, TR[e]); e1f = prod(TR[e], 1, TR[f])
    ok1 = not sub(ef, TR[h]); ok2 = not sub(he, {k: K.scal(mpq(2), a) for k, a in TR[e].items()})
    ok3 = e1f == {V.VAC: (mpq(lev), mpq(0), mpq(0), mpq(0))}
    res[lab] = {'[E,F]=H': ok1, '[H,E]=2E': ok2, 'E(1)F=k': ok3, 'E(1)F': {str(k): str(a) for k, a in e1f.items()}}
S = E.setup(0); T = S['T']
ours = {'e1': E.kv_from_terms([(K.ONE, T['e1'])]), 'e2': E.kv_from_terms([(K.ONE, T['e2'])]),
        'H1': E.kv_from_terms([(K.ONE, T['H1'])]), 'H2': E.kv_from_terms([(K.ONE, T['H2'])])}
res['H1 = psi(P_H)'] = not sub(ours['H1'], TR['PH']); res['H2 = psi(B_H)'] = not sub(ours['H2'], TR['BH'])
res['support PE'] = sorted(str(k[0]) for k in TR['PE']); res['support e1'] = sorted(str(k[0]) for k in ours['e1'])
res['support BE'] = sorted(str(k[0]) for k in TR['BE']); res['support e2'] = sorted(str(k[0]) for k in ours['e2'])
res['coef PE'] = {str(k[0]): str(a) for k, a in TR['PE'].items()}; res['coef e1'] = {str(k[0]): str(a) for k, a in ours['e1'].items()}
res['coef BE'] = {str(k[0]): str(a) for k, a in TR['BE'].items()}; res['coef e2'] = {str(k[0]): str(a) for k, a in ours['e2'].items()}
print(json.dumps(res, indent=1))
pickle.dump(TR, open(HERE + '/jiao_transported.pkl', 'wb'))
