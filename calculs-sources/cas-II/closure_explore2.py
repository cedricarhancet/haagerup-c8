"""Axis B3/B4 generalized (exploration, floating point): each primary type tau carries a finite list of options,
an option being a list of highest weight vectors (one, or two for a type of multiplicity two) chosen among
isotropic lines of the admissible eigenspaces:  dim E = 2 -> common roots of binary forms; dim E = 3 -> common
zeros of the quadrics in P^2 (resultant); multiplicity two -> span of one line in each of two eigenspaces with the
cross condition (u, y_1 v) = 0.  Then the pairwise closure tests (a), (b) of closure_explore.py, and the full test
(c) for the surviving combinations.  Used for the control case (I)."""
import sys, json, time, itertools
import numpy as np
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, isotropic_lines as IL, closure_explore as CE, griess_subalgebras as GS, primaries as PR
S = sp.sqrt(13)
def num(z): return complex(sp.N(z, 30))
def lines_d3(space):
    """isotropic lines in a three-dimensional E (numerical, high precision): Q1 = Q2 = 0 in P^2, Q0 != 0"""
    s, t = sp.symbols('s t')
    Qs = [Q for Q in space['Qs'] if not Q.is_zero_matrix]
    out = []
    for chart in ((1, s, t), (0, 1, t), (0, 0, 1)):
        u = sp.Matrix(chart)
        eqs = [sp.expand((u.T*Q*u)[0, 0]) for Q in Qs]
        vars_ = [v for v in (s, t) if any(e.has(v) for e in eqs)]
        if not vars_:
            if all(e == 0 for e in eqs) and sp.simplify((u.T*space['Q0']*u)[0, 0]) != 0: out.append([num(c) for c in chart])
            continue
        if len(vars_) == 1:
            g = sp.Poly(eqs[0], vars_[0])
            for e in eqs[1:]: g = sp.gcd(g, sp.Poly(e, vars_[0]))
            for r in sp.Poly(g, vars_[0]).nroots(n=30):
                out.append([num(sp.sympify(c).subs(vars_[0], r)) for c in chart])
            continue
        res = sp.resultant(eqs[0], eqs[1], t)
        for r in sp.Poly(res, s).nroots(n=30):
            e0 = sp.Poly(sp.N(eqs[0].subs(s, r), 30), t)
            for r2 in e0.nroots(n=30):
                vals = [abs(complex(sp.N(e.subs({s: r, t: r2}), 20))) for e in eqs]
                if max(vals) < 1e-12:
                    out.append([num(sp.sympify(c).subs({s: r, t: r2})) for c in chart])
    # remove duplicates (projectively)
    uniq = []
    for v in out:
        v = np.array(v); v = v/v[np.argmax(np.abs(v))]
        if not any(np.allclose(v, w, atol=1e-10) for w in uniq): uniq.append(v)
    return [list(v) for v in uniq]
if __name__ == '__main__':
    case, ci = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    D, out = IL.admissible_spaces(case, ci)
    T = D['T']
    N1, N2 = CE.Num(1), CE.Num(2); s13 = float(np.sqrt(13))
    mult = PR.PRIMARY_TYPES[case]
    types = {}
    for tau, v in out.items():
        parts = []; offs = []; opts_spaces = []
        for spc in v['spaces']:
            E = spc['E']; d = E.shape[1]
            offs.append(len(parts))
            parts += [CE.split_vec(list(E[:, k]), v['P']) for k in range(d)]
            L = IL.lines_in(spc) if d <= 2 else None
            if d == 3: Ls = lines_d3(spc)
            elif isinstance(L, list): Ls = [[num(c) for c in l] for l in L]
            else: Ls = []
            opts_spaces.append((spc, Ls))
        nb = len(parts)
        def embed(si, coef):
            x = [0j]*nb
            for k, c in enumerate(coef): x[offs[si] + k] = c
            return x
        options = []
        if mult[tau] == 1:
            for si, (spc, Ls) in enumerate(opts_spaces):
                for l in Ls: options.append([embed(si, l)])
        else:
            # two lines from two different eigenspaces, with the cross isotropy condition
            for (si, (s1, L1)), (sj, (s2, L2)) in itertools.combinations(list(enumerate(opts_spaces)), 2):
                for l1 in L1:
                    for l2 in L2:
                        options.append([embed(si, l1), embed(sj, l2)])
        mods = [[CE.lowerings(T, tau, p) for p in pq] for pq in parts]
        types[tau] = {'parts': parts, 'mods': mods, 'options': options, 'nb': nb}
        print(tau, 'options:', len(options), flush=True)
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    _, basis, _, _ = GS.griess_data(case)
    import commutant_griess as CG
    U1 = CG.U1_basis(T); A2 = []
    for x in U1:
        A2.append(V.clean(V.mode_vec(x, -2, {V.VAC: mpq(1)})))
        for y in U1: A2.append(V.clean(V.mode_vec(x, -1, y)))
    xv = [0] + [num(sp.sympify(c)) for c in G['x']]
    e_num = N2.vec(basis[0]); x_num = sum((xv[i]*N2.vec(basis[i]) for i in range(1, len(basis))), np.zeros(len(N2.keys), dtype=complex))
    U1m = np.array([N1.vec(u) for u in U1]).T; A2m = np.array([N2.vec(a) for a in A2 if a]).T
    tl = sorted(types)
    cols = [A2m, e_num[:, None], x_num[:, None]]
    for tau in tl:
        ty = types[tau]
        for k in range(ty['nb']):
            for w in range(len(ty['mods'][k][0])):
                cols.append((N2.vec(ty['mods'][k][0][w]) + s13*N2.vec(ty['mods'][k][1][w]))[:, None])
    U2max = np.hstack(cols)
    cache = {}
    def pn(tau, k, s_, sigma, k2, s2, w, n):
        key = (tau, k, s_, sigma, k2, s2, w, n)
        if key not in cache:
            a = types[tau]['parts'][k][s_]; b = types[sigma]['mods'][k2][s2][w]
            cache[key] = (N1 if n == 2 else N2).vec(V.clean(V.mode_vec(a, n, b))) if a and b else None
        return cache[key]
    def prods(tau, u, sigma, v, n):
        nw = len(types[sigma]['mods'][0][0]); res = []
        for w in range(nw):
            z = np.zeros(len((N1 if n == 2 else N2).keys), dtype=complex)
            for k in range(types[tau]['nb']):
                if u[k] == 0: continue
                for k2 in range(types[sigma]['nb']):
                    if v[k2] == 0: continue
                    c = u[k]*v[k2]
                    for s_ in (0, 1):
                        for s2 in (0, 1):
                            r = pn(tau, k, s_, sigma, k2, s2, w, n)
                            if r is not None: z += c*(13.0 if (s_, s2) == (1, 1) else (s13 if s_ + s2 == 1 else 1.0))*r
            res.append(z)
        return res
    RES = {}
    def opt_ok(tau, o1, sigma, o2):
        worst = (0.0, 0.0)
        for u in o1:
            for v in o2:
                ra = max(CE.residual(U1m, z) for z in prods(tau, u, sigma, v, 2))
                rb = max(CE.residual(U2max, z) for z in prods(tau, u, sigma, v, 1))
                worst = (max(worst[0], ra), max(worst[1], rb))
        RES[(str(tau), str(sigma))] = RES.get((str(tau), str(sigma)), []) + [worst]
        return worst[0] < 1e-4 and worst[1] < 1e-4
    ok = {}
    for i, tau in enumerate(tl):
        for sigma in tl[i:]:
            for a, o1 in enumerate(types[tau]['options']):
                for b, o2 in enumerate(types[sigma]['options']):
                    if tau == sigma and b != a: continue
                    ok[(tau, a, sigma, b)] = opt_ok(tau, o1, sigma, o2)
            print('pair', tau, sigma, '(%.0fs)' % (time.time() - t0), flush=True)
    good = {t: [a for a in range(len(types[t]['options'])) if ok[(t, a, t, a)]] for t in tl}
    print('self residuals:', {k: v for k, v in RES.items() if k[0] == k[1]})
    print('options passing the self test:', {str(t): v for t, v in good.items()})
    combos = []
    def rec(i, cur):
        if i == len(tl): combos.append(list(cur)); return
        for a in good[tl[i]]:
            if all(ok[(tl[j], cur[j], tl[i], a)] for j in range(i)):
                rec(i + 1, cur + [a])
    rec(0, [])
    print('combinations passing all pairwise tests:', len(combos), combos[:10])
    json.dump({'good': {str(t): v for t, v in good.items()}, 'combos': combos,
               'noptions': {str(t): len(types[t]['options']) for t in tl}}, open(f'closure_explore2_{case}_{ci}.json', 'w'), indent=1)
    print('time %.0fs' % (time.time() - t0))
