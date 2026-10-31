"""Axis E: exact test of phi_chi(U_2) = U_2 over K for the positive character chi of unitary_character.py
(kappa' = kappa'' = 1):  phi(h-monomial e^lam) = (-1)^{deg} chi(lam) conj(coef) h-monomial e^{-lam}.
If phi_chi(U_1) = U_1 and phi_chi(U_2) = U_2, then phi_chi preserves the generated algebra S, and the restriction of
the positive Hermitian form of V' (transported by chi) makes S unitary (Riesen, Proposition 3.1; FLM)."""
import sys, json, time
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V, cases as C, closure_explore as CE, griess_subalgebras as GS, Kfield as K
import exact_isotropy_s13 as ES, exact_weight2_K as EW
CHI = [mpq(2, 15), mpq(1), mpq(1), mpq(5, 3), mpq(9, 25), mpq(5, 24), mpq(2), mpq(2)]
def chi(lam):
    out = mpq(1)
    for c, x in zip(CHI, lam):
        out *= c**x if x >= 0 else (1/c)**(-x)
    return out
def phi(v):
    out = {}
    for (lam, mono), a in v.items():
        k2 = (V.neg(lam), mono)
        out[k2] = K.scal(chi(lam)*(-1)**len(mono), K.conj(a))
    return out
if __name__ == '__main__':
    chi_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    t0 = time.time()
    T, D, out, types, N1, N2, s13, U1, A2 = CE.setup('II', 9)
    # check phi on U_1
    Tc = C.build('II')
    for nm in ('H1', 'e1', 'f1', 'H2', 'e2', 'f2'):
        v = {k: (c, mpq(0), mpq(0), mpq(0)) for k, c in Tc[nm].items()}
        pv = phi(v)
        tgt = {'H1': ('H1', -1), 'e1': ('f1', -1), 'f1': ('e1', -1), 'H2': ('H2', -1), 'e2': ('f2', -1), 'f2': ('e2', -1)}[nm]
        want = {k: (c*tgt[1], mpq(0), mpq(0), mpq(0)) for k, c in Tc[tgt[0]].items()}
        print('phi(%s) = -%s :' % (nm, tgt[0]), pv == want)
    info = json.load(open('closure_explore_II_9.json')); Lj = json.load(open('isotropic_lines_II_9.json'))
    G = json.load(open('griess_subalgebras.json'))['II']['subalgebras'][9]
    _, basis, _, _ = GS.griess_data('II')
    tl = sorted(types); choice = info['combos'][chi_idx]
    lines = {tau: [[K.from_sympy(c) for c in sp.sympify(l)] for l in Lj[str(tau)][0]['lines']] for tau in tl}
    xs = [ES.split(c) for c in G['x']]
    xterms = []
    for i, (p_, q_) in enumerate(xs):
        if p_: xterms.append((K.scal(mpq(int(p_.p), int(p_.q)), K.ONE), basis[i + 1]))
        if q_: xterms.append((K.scal(mpq(int(q_.p), int(q_.q)), K.R), basis[i + 1]))
    U2K = [EW.kvec_from_parts([(K.ONE, a)]) for a in A2 if a] + [EW.kvec_from_parts([(K.ONE, basis[0])]), EW.kvec_from_parts(xterms)]
    for tau, l in zip(tl, choice):
        ty = types[tau]; lv = lines[tau][l]
        for w in range(len(ty['mods'][0][0])):
            terms = []
            for k in range(ty['d']):
                if ty['mods'][k][0][w]: terms.append((lv[k], ty['mods'][k][0][w]))
                if ty['mods'][k][1][w]: terms.append((K.mul(lv[k], K.R), ty['mods'][k][1][w]))
            U2K.append(EW.kvec_from_parts(terms))
    U2Q = [EW.weil(EW.kmul_vec(al, v)) for v in U2K for al in EW.BAS if v]
    byw = {}
    for rw in U2Q: byw.setdefault(C.weight_of_key(T, next(iter(rw))[0]), []).append(rw)
    def cols_of(rows): return sorted(set(c for rw in rows for c in rw))
    rk = {mu: EW.rank_rows(rows, cols_of(rows)) for mu, rows in byw.items()}
    fails = 0
    for v in U2K:
        if not v: continue
        pv = EW.weil(phi(v))
        mu = C.weight_of_key(T, next(iter(pv))[0])
        rows = byw.get(mu, [])
        cols = cols_of(rows + [pv])
        if EW.rank_rows(rows + [pv], cols) > rk.get(mu, 0): fails += 1
    print('choice', choice, ': vectors of U_2 whose image by phi_chi is NOT in U_2:', fails, ' (%.0fs)' % (time.time() - t0))
    json.dump({'choice': choice, 'chi': [str(c) for c in CHI], 'phi_U2_not_in_U2': fails}, open(f'unitarity_test_{chi_idx}.json', 'w'), indent=1)
