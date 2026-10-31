"""Axis B4/C, exact verification modulo a prime p (embedding K -> F_p, see modp.py), for one choice of U_2:
 (W2) dim U_2 = 120 and full weight-two closure: u_1 w' in U_2, u_2 w' in U_1 for all pairs (rank tests mod p);
 (W3) S_3^{gen}: generated weight-three vectors (A_3, y(-1)U_2, L(-1)U_2, u_0 U_2 for u highest weight), closed
      under the zero modes of e', f', e'', f''; its rank mod p (a lower bound for dim S_3 over K; > 956 would be a
      contradiction);
 (C3) u_2 s in U_2 and u_3 s in U_1 for s in S_3^{gen} and u running over highest weight vectors of U_2."""
import sys, json, time
import numpy as np
import scipy.sparse as sps
import sympy as sp
from flint import nmod_mat
from gmpy2 import mpq
import lattice_voa as V, cases as C, closure_explore as CE, griess_subalgebras as GS, exact_isotropy_s13 as ES, modp
import os
F = modp.Fp(int(os.environ.get("PSTART", 2**27))); P = F.p
def rank_cols(cols, n):
    if not cols: return 0
    M = nmod_mat(len(cols), n, [int(x) for c in cols for x in c], P)
    return M.rank()
def sparse_op(vec, n, src_keys, out_idx, nout):
    rows, cols, vals = [], [], []
    for j, k in enumerate(src_keys):
        r = V.clean(V.mode_vec(vec, n, {k: mpq(1)}))
        for kk, c in r.items():
            rows.append(out_idx[kk]); cols.append(j); vals.append(F.q(c))
    return sps.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(nout, len(src_keys)))
def apply(op, v):            # exact (int64, p < 2^27, few terms per row)
    return np.asarray(op @ v) % P
if __name__ == '__main__':
    case, ci, chi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    do3 = len(sys.argv) > 4 and sys.argv[4] == 'w3'
    t0 = time.time()
    T, D, out, types, N1, N2, s13f, U1, A2 = CE.setup(case, ci)
    info = json.load(open(f'closure_explore_{case}_{ci}.json'))
    Lj = json.load(open(f'isotropic_lines_{case}_{ci}.json'))
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    _, basis, _, _ = GS.griess_data(case)
    tl = sorted(types)
    choice = info['combos'][chi] if chi >= 0 else [0]*len(tl)      # chi = -1: a rejected combination (sanity)
    xs = [ES.split(c) for c in G['x']]
    x_p = {}; x_q = {}
    for i, (p_, q_) in enumerate(xs):
        if p_: V.vadd(x_p, basis[i + 1], mpq(int(p_.p), int(p_.q)))
        if q_: V.vadd(x_q, basis[i + 1], mpq(int(q_.p), int(q_.q)))
    x_p, x_q = V.clean(x_p), V.clean(x_q)
    s13 = F.s13
    # U_2 as list of combinations [(coef mod p, rational vector)]
    lines = {tau: [[F.alg(c) for c in sp.sympify(l)] for l in Lj[str(tau)][0]['lines']] for tau in tl}
    U2 = [[(1, a)] for a in A2 if a] + [[(1, basis[0])], [(1, x_p), (s13, x_q)]]
    hw = [('e', U2[-2]), ('x', U2[-1])]
    mods = {}
    for tau, l in zip(tl, choice):
        ty = types[tau]; lv = lines[tau][l]; ml = []
        for w in range(len(ty['mods'][0][0])):
            comb = []
            for k in range(ty['d']):
                if ty['mods'][k][0][w]: comb.append((lv[k], ty['mods'][k][0][w]))
                if ty['mods'][k][1][w]: comb.append((lv[k]*s13 % P, ty['mods'][k][1][w]))
            ml.append(comb); U2.append(comb)
        mods[tau] = ml; hw.append((str(tau), ml[0]))
    keys1, keys2, keys3 = V.basis(1), V.basis(2), V.basis(3)
    idx1 = {k: i for i, k in enumerate(keys1)}; idx2 = {k: i for i, k in enumerate(keys2)}; idx3 = {k: i for i, k in enumerate(keys3)}
    def red(comb, idx, n):
        z = np.zeros(n, dtype=np.int64)
        for c, v in comb:
            z = (z + c*np.array(F.vec(v, idx, n), dtype=np.int64)) % P
        return z
    cache = {}
    def prod(cu, n, cv, idx, nn):
        z = np.zeros(nn, dtype=np.int64)
        for a, va in cu:
            for b, vb in cv:
                key = (id(va), n, id(vb))
                if key not in cache:
                    cache[key] = np.array(F.vec(V.clean(V.mode_vec(va, n, vb)), idx, nn), dtype=np.int64)
                z = (z + (a*b % P)*cache[key]) % P
        return z
    U2m = [red(c, idx2, len(keys2)) for c in U2]
    U1m = [np.array(F.vec(u, idx1, len(keys1)), dtype=np.int64) for u in U1]
    r2 = rank_cols(U2m, len(keys2)); r1 = rank_cols(U1m, len(keys1))
    print('p =', P, '; choice', choice, ': rank U_2 mod p =', r2, '; rank U_1 =', r1, flush=True)
    # (W2) closure: for all pairs (tau <= sigma), u = hw of tau, w' in the module of sigma
    worst = []
    for i, tau in enumerate(tl):
        for sigma in tl[i:]:
            bad1 = bad2 = 0
            for w in mods[sigma]:
                z1 = prod(mods[tau][0], 1, w, idx2, len(keys2))
                z2 = prod(mods[tau][0], 2, w, idx1, len(keys1))
                if np.any(z1) and rank_cols(U2m + [z1], len(keys2)) > r2: bad1 += 1
                if np.any(z2) and rank_cols(U1m + [z2], len(keys1)) > r1: bad2 += 1
            if bad1 or bad2: worst.append((str(tau), str(sigma), bad1, bad2))
    # also products with x (C_2) and among C_2: x_1 u in U_2, x_2 u in U_1
    for tau in tl:
        for w in mods[tau]:
            for n, idx, nn, base, rr in ((1, idx2, len(keys2), U2m, r2), (2, idx1, len(keys1), U1m, r1)):
                z = prod(U2[len(A2) + 1] if False else hw[1][1], n, w, idx, nn)
                if np.any(z) and rank_cols(base + [z], nn) > rr: worst.append(('x', str(tau), n))
    print('(W2) weight-two closure failures mod p:', worst if worst else 'none', ' (%.0fs)' % (time.time() - t0), flush=True)
    res = {'p': P, 'choice': choice, 'rankU2': r2, 'W2_failures': worst}
    if do3:
        n3 = len(keys3)
        gens = []
        for x in U1:
            gens.append(np.array(F.vec(V.clean(V.mode_vec(x, -3, {V.VAC: mpq(1)})), idx3, n3), dtype=np.int64))
            for y in U1: gens.append(np.array(F.vec(V.clean(V.mode_vec(x, -2, y)), idx3, n3), dtype=np.int64))
        for u in U2:
            gens.append(prod([(1, V.OMEGA)], 0, u, idx3, n3))
            for y in U1: gens.append(prod([(1, y)], -1, u, idx3, n3))
        for name, u in hw + [('A2', [(1, a)]) for a in A2 if a]:
            for v in U2: gens.append(prod(u, 0, v, idx3, n3))
        print('weight-three generators:', len(gens), '(%.0fs)' % (time.time() - t0), flush=True)
        cache.clear(); V._MODE.clear(); V._EMODE.clear()
        # sort by weight
        wts = [C.weight_of_key(T, k) for k in keys3]
        def wof(v):
            nz = np.nonzero(v)[0]
            return wts[nz[0]] if len(nz) else None
        byw = {}
        for g in gens:
            mu = wof(g)
            if mu is not None: byw.setdefault(mu, []).append(g)
        def reduce_basis(vs, mu):
            """row-reduced independent subset (exact mod p) restricted to the coordinates of weight mu"""
            cols = [i for i in range(n3) if wts[i] == mu]
            M = nmod_mat(len(vs), len(cols), [int(v[i]) for v in vs for i in cols], P)
            R, rk = M.rref()[0], M.rank()
            out_ = []
            for r_ in range(rk):
                z = np.zeros(n3, dtype=np.int64)
                for jj, i in enumerate(cols): z[i] = int(R[r_, jj])
                out_.append(z)
            return out_
        bases = {mu: reduce_basis(vs, mu) for mu, vs in byw.items()}
        print('rank before g-closure:', sum(len(b) for b in bases.values()), flush=True)
        ops = [sparse_op(T[k], 0, keys3, idx3, n3) for k in ('f1', 'f2', 'e1', 'e2')]
        for it in range(20):
            new = {}
            for mu, B in bases.items():
                for op in ops:
                    for v in B:
                        z = apply(op, v)
                        m2 = wof(z)
                        if m2 is not None: new.setdefault(m2, []).append(z)
            grown = False
            for mu, vs in new.items():
                old = bases.get(mu, [])
                nb = reduce_basis(old + vs, mu)
                if len(nb) > len(old): grown = True
                bases[mu] = nb
            tot = sum(len(b) for b in bases.values())
            print('  g-closure iteration', it, ': rank', tot, flush=True)
            if not grown: break
        dimS3 = sum(len(b) for b in bases.values())
        res['dimS3_gen_modp'] = dimS3
        S3 = [v for B in bases.values() for v in B]
        # module structure from dominant weights
        per = {mu: len(B) for mu, B in bases.items()}
        def rr(a, b): return per.get((a, b), 0)
        mult = {f'{a},{b}': rr(a, b) - rr(a + 2, b) - rr(a, b + 2) + rr(a + 2, b + 2) for (a, b) in per if a >= 0 and b >= 0}
        res['S3_structure'] = {k: v for k, v in mult.items() if v}
        print('S_3^{gen} (g-closed) rank mod p =', dimS3, ' (dim U_3 = 956) (%.0fs)' % (time.time() - t0), flush=True)
        # (C3) closure tests
        S3m = np.array(S3, dtype=np.int64).T
        fails = {}
        C3names = [] if os.environ.get('NOC3') else os.environ.get('C3NAMES', 'x;(3, 3);(1, 1)').split(';')
        for name, u in [h for h in hw if h[0] in C3names]:
            V._MODE.clear(); V._EMODE.clear()
            op2 = None; op3 = None
            for a, va in u:
                o2 = sparse_op(va, 2, keys3, idx2, len(keys2)); o3 = sparse_op(va, 3, keys3, idx1, len(keys1))
                Z2 = (a*(np.asarray(o2 @ S3m) % P)) % P; Z3 = (a*(np.asarray(o3 @ S3m) % P)) % P
                op2 = Z2 if op2 is None else (op2 + Z2) % P; op3 = Z3 if op3 is None else (op3 + Z3) % P
            ex2 = rank_cols(U2m + [op2[:, j] for j in range(op2.shape[1])], len(keys2)) - r2
            ex3 = rank_cols(U1m + [op3[:, j] for j in range(op3.shape[1])], len(keys1)) - r1
            fails[name] = (int(ex2), int(ex3))
            print('  u =', name, ': extra rank of u_2 S_3 over U_2 =', ex2, '; of u_3 S_3 over U_1 =', ex3, ' (%.0fs)' % (time.time() - t0), flush=True)
        if os.environ.get('C3U1'):
            for i_, y in enumerate(U1):
                V._MODE.clear(); V._EMODE.clear()
                o1 = sparse_op(y, 1, keys3, idx2, len(keys2)); o2 = sparse_op(y, 2, keys3, idx1, len(keys1))
                Z1 = np.asarray(o1 @ S3m) % P; Z2 = np.asarray(o2 @ S3m) % P
                ex1 = rank_cols(U2m + [Z1[:, j] for j in range(Z1.shape[1])], len(keys2)) - r2
                ex2 = rank_cols(U1m + [Z2[:, j] for j in range(Z2.shape[1])], len(keys1)) - r1
                fails['U1_%d' % i_] = (int(ex1), int(ex2))
                print('  y = U_1[%d] : extra rank of y_1 S_3 over U_2 =' % i_, ex1, '; of y_2 S_3 over U_1 =', ex2, flush=True)
        res['C3_extra_ranks'] = fails
    json.dump(res, open(f'closure_modp_{case}_{ci}_{chi}{"_w3" if do3 else ""}_p{P}.json', 'w'), indent=1)
    print('time %.0fs' % (time.time() - t0))
