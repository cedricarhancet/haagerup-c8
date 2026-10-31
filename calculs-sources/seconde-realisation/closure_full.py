"""Axis B4 (exploration): full weight-two closure for the choices surviving closure_explore.py:
   u_1 w' in U_2(choice) and u_2 w' in U_1 for all pairs, plus dim U_2(choice) = 120."""
import sys, json, time, itertools
import numpy as np
import closure_explore as CE, griess_subalgebras as GS
import sympy as sp
if __name__ == '__main__':
    case, ci = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    T, D, out, types, N1, N2, s13, U1, A2 = CE.setup(case, ci)
    info = json.load(open(f'closure_explore_{case}_{ci}.json'))
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    _, basis, _, _ = GS.griess_data(case)
    xv = [0] + [complex(sp.N(sp.sympify(c), 30)) for c in G['x']]
    e_num = N2.vec(basis[0]); x_num = sum((xv[i]*N2.vec(basis[i]) for i in range(1, len(basis))), np.zeros(len(N2.keys), dtype=complex))
    tl = sorted(types)
    U1m = np.array([N1.vec(u) for u in U1]).T
    A2m = np.array([N2.vec(a) for a in A2]).T
    def line_module(tau, l):
        ty = types[tau]; lv = ty['lines'][l]
        return [sum((lv[k]*(N2.vec(ty['mods'][k][0][w]) + s13*N2.vec(ty['mods'][k][1][w])) for k in range(ty['d'])), np.zeros(len(N2.keys), dtype=complex)) for w in range(len(ty['mods'][0][0]))]
    import lattice_voa as V
    for choice in info['combos']:
        mods = {tau: line_module(tau, l) for tau, l in zip(tl, choice)}
        U2 = np.hstack([A2m, e_num[:, None], x_num[:, None]] + [np.array(mods[t]).T for t in tl])
        r = np.linalg.matrix_rank(U2, tol=1e-9)
        worst_b, worst_a = 0, 0
        # products of highest weight vectors u with module vectors (numerically, from exact rational parts)
        for i, tau in enumerate(tl):
            ty = types[tau]
            for sigma in tl[i:]:
                tys = types[sigma]
                for w in range(len(tys['mods'][0][0])):
                    for n, Bm in ((1, U2), (2, U1m)):
                        z = 0
                        for k in range(ty['d']):
                            for k2 in range(tys['d']):
                                c = ty['lines'][choice[tl.index(tau)]][k]*tys['lines'][choice[tl.index(sigma)]][k2]
                                for s in (0, 1):
                                    for s2 in (0, 1):
                                        a = ty['parts'][k][s]; bb = tys['mods'][k2][s2][w]
                                        if not a or not bb: continue
                                        prod_ = V.clean(V.mode_vec(a, n, bb))
                                        f = (13.0 if (s, s2) == (1, 1) else (s13 if s + s2 == 1 else 1.0))
                                        z = z + c*f*(N2 if n == 1 else N1).vec(prod_)
                        res = CE.residual(Bm, z) if not isinstance(z, int) else 0.0
                        if n == 1: worst_b = max(worst_b, res)
                        else: worst_a = max(worst_a, res)
        # also the products with C_2: x_1 u in U_2 and x_2 u = 0 for u primary (sanity)
        print('choice', choice, ': dim U_2 =', r, '; worst residual u_1 w\' (mod U_2) = %.2e ; u_2 w\' (mod U_1) = %.2e' % (worst_b, worst_a), flush=True)
    print('time %.0fs' % (time.time() - t0))
