"""Axis B4, rigorous part: certified exclusion of the line combinations rejected by closure_explore.py.
For a rejected pair (tau, l; sigma, m) we exhibit a rational linear form phi vanishing on U_1 (test a), or on the
rational span W_Q of A_2, e, x_p, x_q and of all module vectors of all E_rho (a superset of U2max, test b), such
that phi(z) != 0 for some product z = u_n w' (u in the line l of tau, w' in the module of the line m of sigma).
phi(z) = sum_j c_j phi(R_j) with R_j exact rational products and c_j algebraic coefficients evaluated as Arb
balls: the exclusion is certified when the ball does not contain 0."""
import sys, json, time
import sympy as sp
from flint import fmpq_mat, fmpq, acb, arb, ctx
from gmpy2 import mpq
import lattice_voa as V, closure_explore as CE, griess_subalgebras as GS, isotropic_lines as IL, linalg as LA
import exact_isotropy_s13 as ES, isotropy_test as IT
ctx.prec = 256
def annihilator(vecs, keys):
    """rational basis of linear forms (dicts key -> fmpq) vanishing on span(vecs)"""
    kidx = {k: i for i, k in enumerate(keys)}
    M = fmpq_mat(len(vecs), len(keys))
    for i, v in enumerate(vecs):
        for k, c in v.items(): M[i, kidx[k]] = LA.to_fq(c)
    ns = LA.nullspace(M)
    return [{keys[j]: x[j] for j in range(len(keys)) if x[j] != 0} for x in ns]
def phi_of(phi, v):
    s = fmpq(0)
    for k, c in v.items():
        f = phi.get(k)
        if f is not None: s += f*LA.to_fq(c)
    return s
if __name__ == '__main__':
    case, ci = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    T, D, out, types, N1, N2, s13, U1, A2 = CE.setup(case, ci)
    info = json.load(open(f'closure_explore_{case}_{ci}.json'))
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras'][ci]
    _, basis, _, _ = GS.griess_data(case)
    xs = [ES.split(c) for c in G['x']]
    x_p = {}; x_q = {}
    for i, (p, q) in enumerate(xs):
        if p: V.vadd(x_p, basis[i + 1], mpq(int(p.p), int(p.q)))
        if q: V.vadd(x_q, basis[i + 1], mpq(int(q.p), int(q.q)))
    tl = sorted(types)
    # exact line coefficients (sympy) -> acb
    L = json.load(open(f'isotropic_lines_{case}_{ci}.json'))
    def acb_of(expr):
        z = sp.N(sp.sympify(expr), 90)
        re, im = sp.re(z), sp.im(z)
        return acb(arb(str(sp.Float(re, 90))), arb(str(sp.Float(im, 90))))
    lines = {}
    for tau in tl:
        lines[tau] = [[acb_of(c) for c in sp.sympify(l.replace('I', '1j').replace('1j', 'I'))] if False else
                      [acb_of(c) for c in sp.sympify(l)] for l in L[str(tau)][0]['lines']]
    S13 = arb(13).sqrt()
    W = [a for a in A2 if a] + [basis[0], x_p] + ([x_q] if x_q else [])
    for tau in tl:
        for k in range(types[tau]['d']):
            for s in (0, 1):
                W += [v for v in types[tau]['mods'][k][s] if v]
    keys2 = V.basis(2); keys1 = V.basis(1)
    AnnW = annihilator(W, keys2); AnnU1 = annihilator(U1, keys1)
    print('dim of the rational span W_Q:', len(keys2) - len(AnnW), '; forms vanishing on W_Q:', len(AnnW), '; on U_1:', len(AnnU1), flush=True)
    rejected = [r for r in info['report'] if not (float(r[4]) < 1e-8 and float(r[5]) < 1e-8)]
    certified = 0; failed = []
    for (st, l, ss, m, ra, rb) in rejected:
        tau, sigma = eval(st), eval(ss)
        ok = False
        tests = []
        if float(ra) >= 1e-8: tests.append((2, AnnU1))
        if float(rb) >= 1e-8: tests.append((1, AnnW))
        for (n, Ann) in tests:
            nw = len(types[sigma]['mods'][0][0])
            for w in range(nw):
                # exact rational products R and their coefficients (balls)
                terms = []
                for k in range(types[tau]['d']):
                    for k2 in range(types[sigma]['d']):
                        c = lines[tau][l][k]*lines[sigma][m][k2]
                        for s in (0, 1):
                            for s2 in (0, 1):
                                a = types[tau]['parts'][k][s]; bb = types[sigma]['mods'][k2][s2][w]
                                if not a or not bb: continue
                                R = V.clean(V.mode_vec(a, n, bb))
                                if not R: continue
                                f = (acb(13) if (s, s2) == (1, 1) else (acb(S13) if s + s2 == 1 else acb(1)))
                                terms.append((c*f, R))
                if not terms: continue
                for phi in Ann[:60]:
                    val = acb(0)
                    for c, R in terms:
                        pv = phi_of(phi, R)
                        if pv != 0: val += c*acb(arb(int(pv.p))/int(pv.q))
                    if not val.contains(0):
                        ok = True; break
                if ok: break
            if ok: break
        if ok: certified += 1
        else: failed.append((st, l, ss, m))
    print('rejected pairs:', len(rejected), '; certified by an exact form + Arb ball:', certified, '; not certified:', failed)
    json.dump({'rejected': len(rejected), 'certified': certified, 'not_certified': failed}, open(f'closure_certify_{case}_{ci}.json', 'w'), indent=1)
    print('time %.0fs' % (time.time() - t0))
