"""Axis B2: for every candidate C_2 = span(e, x) (Axis A2) and every primary type tau of U_2, the admissible
lines p_tau = U_2 cap Pi_tau:
   (K) e_2 u = x_2 u = 0  (L_1 u and x_2 u are weight-one A-primaries of U, and U_1 has none),
   (E) u is an eigenvector of x_1 on Pi_tau.
Method: the eigenvectors of X = x_1 lying in K = ker(e_2) cap ker(x_2) are the eigenvectors of X restricted to
W = largest X-stable subspace of K = ker [K; K X; K X^2; ...]  (exact linear algebra).
Coefficients of x in Q or Q(sqrt 13): an element A + sqrt13 B acts on Q(sqrt13)^n = Q^{2n} by [[A,13B],[B,A]].
A candidate with dim W < (multiplicity of tau in U_2) for some tau is EXCLUDED (exact, no unitarity)."""
import sys, json, time, pickle
import sympy as sp
from flint import fmpq_mat, fmpq
s13 = sp.sqrt(13)
def h_tau(case, tau):
    k1, k2 = {'I': (1, 39), 'II': (3, 13)}[case]
    a, b = tau
    return sp.Rational(a*(a + 2), 4*(k1 + 2)) + sp.Rational(b*(b + 2), 4*(k2 + 2))
def load(case):
    d = pickle.load(open(f'primaries_{case}.pkl', 'rb'))
    out = {}
    for t, v in d['types'].items():
        out[eval(t)] = {'X1': [eval(m, sp.__dict__) for m in v['X1']], 'X2': [eval(m, sp.__dict__) for m in v['X2']],
                        'n': len(v['P'])}
    return out
def split(z):
    """z in Q(sqrt13) -> (p, q) rationals with z = p + q sqrt13"""
    z = sp.nsimplify(sp.expand(z))
    q = sp.Rational(sp.expand(z).coeff(s13)) if z.has(s13) else sp.Rational(0)
    p = sp.Rational(sp.expand(z - q*s13))
    return p, q
def fq(r): return fmpq(int(r.p), int(r.q))
def double(Mp, Mq):
    """flint matrix of A + sqrt13 B acting on Q^{2n}"""
    r, c = Mp.rows, Mp.cols
    M = fmpq_mat(2*r, 2*c)
    for i in range(r):
        for j in range(c):
            a, b = Mp[i, j], Mq[i, j]
            if a: M[i, j] = fq(a); M[i + r, j + c] = fq(a)
            if b: M[i, j + c] = fq(13*b); M[i + r, j] = fq(b)
    return M
def nullspace(M):
    Rr, rank = M.rref(); n = M.ncols(); piv = []
    for i in range(rank): piv.append(next(j for j in range(n) if Rr[i, j] != 0))
    free = [j for j in range(n) if j not in piv]; out = []
    for fj in free:
        x = [fmpq(0)]*n; x[fj] = fmpq(1)
        for i, pj in enumerate(piv): x[pj] = -Rr[i, fj]
        out.append(x)
    return out
def stable_subspace(X, K, n2):
    """largest X-stable subspace of ker K (flint matrices on Q^{n2})"""
    rows = []
    P = fmpq_mat(n2, n2, [fmpq(1) if i == j else fmpq(0) for i in range(n2) for j in range(n2)])
    blocks = []
    dim_prev = None
    for it in range(n2 + 1):
        blocks.append(K*P)
        S = fmpq_mat(sum(b.nrows() for b in blocks), n2)
        r0 = 0
        for b in blocks:
            for i in range(b.nrows()):
                for j in range(n2): S[r0 + i, j] = b[i, j]
            r0 += b.nrows()
        ns = nullspace(S)
        if dim_prev is not None and len(ns) == dim_prev:
            return ns
        dim_prev = len(ns)
        if not ns: return []
        P = X*P
    return ns
def analyse(case, cands, types, mult):
    res = []
    for ci, cand in enumerate(cands):
        xs = [split(sp.sympify(c)) for c in cand['x']]
        irr = any(q != 0 for p, q in xs)
        lam = sp.sympify(cand['lambda']); D = sp.sympify(cand['D'])
        row = {'c': cand.get('c'), 'x': cand['x'], 'field': 'Q(sqrt13)' if irr else 'Q', 'types': {}, 'excluded_by': []}
        for tau, d in types.items():
            n = d['n']
            Xp = sum((xs[i][0]*d['X1'][i + 1] for i in range(len(xs))), sp.zeros(n, n))
            Xq = sum((xs[i][1]*d['X1'][i + 1] for i in range(len(xs))), sp.zeros(n, n))
            r2 = d['X2'][0].rows
            X2p = sum((xs[i][0]*d['X2'][i + 1] for i in range(len(xs))), sp.zeros(r2, n))
            X2q = sum((xs[i][1]*d['X2'][i + 1] for i in range(len(xs))), sp.zeros(r2, n))
            Kp = d['X2'][0].col_join(X2p) if r2 else sp.zeros(1, n)
            Kq = sp.zeros(r2, n).col_join(X2q) if r2 else sp.zeros(1, n)
            X = double(Xp, Xq); K = double(Kp, Kq)
            W = stable_subspace(X, K, 2*n)
            dimW = len(W)//2          # dimension over Q(sqrt13)
            info = {'dimPi': n, 'dimW': dimW}
            if 0 < dimW <= 3:
                # eigenvalues of X on W: characteristic polynomial of the restriction (over Q, doubled)
                Wm = sp.Matrix([[sp.Rational(int(w[j].p), int(w[j].q)) for w in W] for j in range(2*n)])
                Xs = sp.Matrix(2*n, 2*n, lambda i, j: sp.Rational(int(X[i, j].p), int(X[i, j].q)))
                XW = Xs*Wm
                coef = Wm.solve_least_squares(XW) if False else (Wm.T*Wm).inv()*Wm.T*XW
                t = sp.Symbol('t')
                cp = sp.factor((coef - t*sp.eye(coef.rows)).det())
                info['charpoly_on_W_over_Q'] = str(cp)
            row['types'][str(tau)] = info
            if dimW < mult[tau]:
                row['excluded_by'].append(str(tau))
        res.append(row)
    return res
if __name__ == '__main__':
    G = json.load(open('griess_subalgebras.json'))
    import primaries as PR
    out = {}
    for case in sys.argv[1:] or ['I', 'II']:
        t0 = time.time()
        types = load(case)
        res = analyse(case, G[case]['subalgebras'], types, PR.PRIMARY_TYPES[case])
        out[case] = res
        print('case', case)
        for i, r in enumerate(res):
            tag = ('EXCLUDED: no admissible line for ' + ', '.join(r['excluded_by'])) if r['excluded_by'] else 'survives B2'
            print('  cand %2d  c = %-50s %-9s %s' % (i, r['c'], r['field'], tag))
            print('          dim W by type:', {t: v['dimW'] for t, v in r['types'].items()})
        print('  time %.1fs' % (time.time() - t0))
    json.dump(out, open('eigen_reduction.json', 'w'), indent=1)
