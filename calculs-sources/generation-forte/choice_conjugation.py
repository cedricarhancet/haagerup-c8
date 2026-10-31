"""Axis A3 (exact over K): relation between the two surviving choices of U_2 for the Griess subalgebra 9.
Test: complex conjugation (q -> -q) maps U_2^{(0)} onto U_2^{(1)}.  Since phi_chi = conj o theta_chi stabilizes each
U_2^{(i)} (unitarity_test.py), this is equivalent to theta_chi(U_2^{(0)}) = U_2^{(1)} for the LINEAR automorphism
theta_chi : e^lam -> (-1)^{#} chi(lam) e^{-lam} of V', which normalizes U_1 (theta_chi(e') = -f', ...)."""
import json, time
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V, cases as C, closure_explore as CE, griess_subalgebras as GS, Kfield as K
import exact_isotropy_s13 as ES, exact_weight2_K as EW
def build_U2K(chi_idx, T, types, A2, sub=9):
    info = json.load(open('closure_explore_II_9.json')); Lj = json.load(open('isotropic_lines_II_9.json'))
    G = json.load(open('griess_subalgebras.json'))['II']['subalgebras'][sub]
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
    return [v for v in U2K if v], choice
def span_data(T, vecs):
    rows = [EW.weil(EW.kmul_vec(al, v)) for v in vecs for al in EW.BAS]
    byw = {}
    for rw in rows: byw.setdefault(C.weight_of_key(T, next(iter(rw))[0]), []).append(rw)
    return byw
def contained(T, byw, vecs):
    bad = 0
    for v in vecs:
        z = EW.weil(v); mu = C.weight_of_key(T, next(iter(z))[0]); rows = byw.get(mu, [])
        cols = sorted(set(c for rw in rows + [z] for c in rw))
        if EW.rank_rows(rows + [z], cols) > EW.rank_rows(rows, cols): bad += 1
    return bad
if __name__ == '__main__':
    t0 = time.time()
    T, D, out, types, N1, N2, s13, U1, A2 = CE.setup('II', 9)
    U0, c0 = build_U2K(0, T, types, A2); U1_, c1 = build_U2K(1, T, types, A2)
    b1 = span_data(T, U1_); b0 = span_data(T, U0)
    conjU0 = [{k: K.conj(a) for k, a in v.items()} for v in U0]
    sig = lambda a: (a[0], -a[1], a[2], -a[3])          # r -> -r (sqrt13 -> -sqrt13), q fixed
    r1 = contained(T, b1, conjU0); r0 = contained(T, b0, conjU0)
    print('conj(U_2^(0)) not in U_2^(1):', r1, ' ; conj(U_2^(0)) not in U_2^(0):', r0, ' (%.0fs)' % (time.time() - t0))
    json.dump({'choice0': c0, 'choice1': c1, 'conjU0_not_in_U1': r1, 'conjU0_not_in_U0': r0}, open('choice_conjugation.json', 'w'), indent=1)
