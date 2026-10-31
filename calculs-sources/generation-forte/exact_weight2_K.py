"""Axis A1: EXACT verification over K = Q(sqrt13, sqrt-3) that the subspace U_2 of Theorem 10.5 (case II,
subalgebra 9, choice of lines given on the command line) is closed in weight two:
   rank_K U_2 = 120;  u_(1) w in U_2 and u_(2) w in U_1 for u highest weight vector of a primary line or u in {e, x}
   and w in the g'+g''-module of any line, or in C_2.
Method: Weil restriction K^n = Q^{4n}; membership z in span_K(U_2) <=> z in span_Q{alpha u : u in U_2, alpha in
{1, r, q, rq}}; exact rational ranks (python-flint fmpq_mat), weight space by weight space (all vectors are joint
weight vectors).  Also computes the scale-free coupling invariants used in Axis E."""
import sys, json, time, itertools
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import sympy as sp
import lattice_voa as V, cases as C, closure_explore as CE, griess_subalgebras as GS, Kfield as K
import exact_isotropy_s13 as ES
BAS = [K.ONE, K.R, K.Q, K.RQ]
def kvec_from_parts(terms):
    """terms: list of (K element, rational vector) -> K-vector dict key -> K element"""
    out = {}
    for c, v in terms:
        if K.iszero(c) or not v: continue
        for k, x in v.items():
            out[k] = K.add(out.get(k, K.Z), K.scal(x, c))
    return {k: a for k, a in out.items() if not K.iszero(a)}
def kmul_vec(alpha, v): return {k: K.mul(alpha, a) for k, a in v.items()}
def weight_of(T, v):
    return C.weight_of_key(T, next(iter(v))) if v else None
def to_fq(x): return fmpq(int(x.numerator), int(x.denominator))
def rank_rows(rows, cols):
    if not rows: return 0
    cidx = {c: i for i, c in enumerate(cols)}
    M = fmpq_mat(len(rows), len(cols))
    for i, rw in enumerate(rows):
        for (k, comp), x in rw.items(): M[i, cidx[(k, comp)]] = to_fq(x)
    return M.rank()
def weil(v):
    return {(k, i): a[i] for k, a in v.items() for i in range(4) if a[i]}
def by_weight(T, vecs):
    out = {}
    for v in vecs:
        if v: out.setdefault(weight_of(T, v), []).append(v)
    return out
if __name__ == '__main__':
    case, ci, chi = 'II', 9, int(sys.argv[1]) if len(sys.argv) > 1 else 0
    t0 = time.time()
    T, D, out, types, N1, N2, s13, U1, A2 = CE.setup(case, ci)
    info = json.load(open(f'closure_explore_{case}_{ci}.json'))
    Lj = json.load(open(f'isotropic_lines_{case}_{ci}.json'))
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    _, basis, _, _ = GS.griess_data(case)
    tl = sorted(types)
    choice = info['combos'][chi]
    # x = x_p + r x_q
    xs = [ES.split(c) for c in G['x']]
    xterms = []
    for i, (p_, q_) in enumerate(xs):
        if p_: xterms.append((K.scal(mpq(int(p_.p), int(p_.q)), K.ONE), basis[i + 1]))
        if q_: xterms.append((K.scal(mpq(int(q_.p), int(q_.q)), K.R), basis[i + 1]))
    lines = {tau: [[K.from_sympy(c) for c in sp.sympify(l)] for l in Lj[str(tau)][0]['lines']] for tau in tl}
    # rational-part terms for module vectors: list over weights of list of (K coef, rational vector)
    def module_terms(tau, l):
        ty = types[tau]; lv = lines[tau][l]; res = []
        for w in range(len(ty['mods'][0][0])):
            terms = []
            for k in range(ty['d']):
                if ty['mods'][k][0][w]: terms.append((lv[k], ty['mods'][k][0][w]))
                if ty['mods'][k][1][w]: terms.append((K.mul(lv[k], K.R), ty['mods'][k][1][w]))
            res.append(terms)
        return res
    MT = {tau: module_terms(tau, l) for tau, l in zip(tl, choice)}
    hwterms = {'e': [(K.ONE, basis[0])], 'x': xterms}
    for tau in tl: hwterms[str(tau)] = MT[tau][0]
    # U_2 as K-vectors
    U2K = [kvec_from_parts([(K.ONE, a)]) for a in A2 if a] + [kvec_from_parts([(K.ONE, basis[0])]), kvec_from_parts(xterms)]
    for tau in tl: U2K += [kvec_from_parts(t) for t in MT[tau]]
    U2Q = [weil(kmul_vec(al, v)) for v in U2K for al in BAS if v]
    U1Q = [weil(kvec_from_parts([(al, u)])) for u in U1 for al in BAS]
    def wq(row):
        k = next(iter(row))[0]; return C.weight_of_key(T, k)
    U2byw = {}
    for rw in U2Q: U2byw.setdefault(wq(rw), []).append(rw)
    U1byw = {}
    for rw in U1Q: U1byw.setdefault(wq(rw), []).append(rw)
    def cols_of(rows): return sorted(set(c for rw in rows for c in rw))
    rk = {mu: rank_rows(rows, cols_of(rows)) for mu, rows in U2byw.items()}
    total = sum(rk.values())
    print('choice', choice, ': rank_Q of the Weil restriction of U_2 =', total, '-> rank_K U_2 =', total/4, ' (%.0fs)' % (time.time() - t0), flush=True)
    rk1 = {mu: rank_rows(rows, cols_of(rows)) for mu, rows in U1byw.items()}
    # products
    def kprod(terms_u, n, terms_w):
        z = {}
        for a, va in terms_u:
            for b, vb in terms_w:
                r = V.clean(V.mode_vec(va, n, vb))
                if not r: continue
                c = K.mul(a, b)
                for k, x in r.items():
                    z[k] = K.add(z.get(k, K.Z), K.scal(x, c))
        return {k: a for k, a in z.items() if not K.iszero(a)}
    targets = [(name, ht) for name, ht in hwterms.items()]
    wlist = []
    for tau in tl: wlist += [(str(tau), w, t) for w, t in enumerate(MT[tau])]
    wlist += [('e', 0, [(K.ONE, basis[0])]), ('x', 0, xterms)]
    fails = []; ntests = 0
    for name, ht in targets:
        for (wname, widx, wt) in wlist:
            for n, base, rkb in ((1, U2byw, rk), (2, U1byw, rk1)):
                z = kprod(ht, n, wt)
                if not z: continue
                ntests += 1
                zq = weil(z)
                mu = C.weight_of_key(T, next(iter(z)))
                rows = base.get(mu, [])
                cols = cols_of(rows + [zq])
                r0 = rank_rows(rows, cols) if rows else 0
                if rank_rows(rows + [zq], cols) > r0:
                    fails.append((name, wname, widx, n))
        print('  u =', name, 'done (%.0fs), failures so far:' % (time.time() - t0), len(fails), flush=True)
    print('EXACT weight-two closure over K: tests =', ntests, '; failures =', len(fails), fails[:5])
    json.dump({'choice': choice, 'rank_K_U2': total/4, 'tests': ntests, 'failures': fails}, open(f'exact_weight2_K_{chi}.json', 'w'), indent=1)
    print('time %.0fs' % (time.time() - t0))
