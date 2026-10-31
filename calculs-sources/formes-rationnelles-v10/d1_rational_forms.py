"""Axis D (question 16.2): Galois descent data for the two Haagerup algebras, K = Q(sqrt13, sqrt-3),
Gal(K/Q) = {1, sigma, tau, sigma tau}, sigma: sqrt13 -> -sqrt13, tau = complex conjugation.  V' has the rational form
V'_Q (states h_1(-n_1)...h_k(-n_k) e^lambda), and g~ = id (x) g acts semilinearly on V'_K.  We look for
psi_g = a_g o g~ with a_g in Z (automorphisms of V' fixing U_1 pointwise, Z = (Z/2)^3, all defined over Q) such that
psi_g(U) = U.  Since Z is elementary abelian and its elements commute with the g~, g -> psi_g is then a homomorphism,
i.e. a Galois descent datum, and U_Q = {u : psi_g(u) = u for all g} is a Q-form of U (Speiser's lemma in each weight).
 case II (python3 d1_rational_forms.py II): a_sigma = z (bits 110), a_tau = z z' (z' = bits 011); EXACT checks over K:
   z, z' rational, z^2 = z'^2 = 1, z z' = z' z on V'_1 (which determines the automorphisms), and
   psi_sigma(U_2) c U_2, psi_tau(U_2) c U_2 for U_2 = U_2(9,0) (exact membership, exact_engine.solve_block);
   explicit Q-form of U_2: Q-rank 120 of the traces sum_g psi_g(a u) (a in {1, r, q, rq}) and K-rank 120;
   psi_sigma, psi_tau fix the two idempotents f_1, f_2 of C_2.
 case I (python3 d1_rational_forms.py I): Riesen's generators e, f, v (transported, r1); a_tau = z_1 = (-1,-1,1),
   a_sigma in {z_1 z_A, z_A}, where z_A = (-1,1,1) is the lift of the diagram automorphism of A_2 with trivial signs
   (notation of Proposition prop:centralizer-one) and z_1 z_A = (1,-1,1) is the lift of the diagram automorphism of
   E_6 with signs; EXACT (sympy) checks that psi_g(e) = e, psi_g(f) = f, psi_g(v) in C v.  (Labels: eps = 2*bit - 1.)"""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); CASE = sys.argv[1]
GF = os.path.abspath(HERE + ('/../generation-forte' if CASE == 'II' else '/../cas-II'))
sys.path.insert(0, GF); os.chdir(GF)
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V
res = {'case': CASE}
def load_Z(path): return sp.Matrix([[sp.sympify(x) for x in r] for r in pickle.load(open(path, 'rb'))])
def lift_data(Zm):
    B1 = V.basis(1); n = len(B1)
    hkey = {k[1][0][1]: i for i, k in enumerate(B1) if k[0] == V.ZERO}
    g = [[int(Zm[hkey[j], hkey[i]]) for i in range(8)] for j in range(8)]
    gv = lambda lam: tuple(sum(g[j][i]*lam[i] for i in range(8)) for j in range(8))
    croot = {}
    for jcol, k in enumerate(B1):
        if k[0] == V.ZERO: continue
        rows = [i for i in range(n) if Zm[i, jcol] != 0]; assert len(rows) == 1
        assert B1[rows[0]][0] == gv(k[0]); croot[k[0]] = Zm[rows[0], jcol]
    simple = [tuple(1 if t == i else 0 for t in range(8)) for i in range(8)]
    cs = [croot[a] for a in simple]
    delta = [[V.eps(simple[i], simple[j])*V.eps(gv(simple[i]), gv(simple[j])) for j in range(8)] for i in range(8)]
    def c(lam):
        v = sp.Integer(1)
        for i in range(8):
            v *= cs[i]**lam[i]
            if delta[i][i] == -1 and (lam[i]*(lam[i] - 1)//2) % 2: v = -v
            for j in range(i + 1, 8):
                if delta[i][j] == -1 and (lam[i]*lam[j]) % 2: v = -v
        return v
    assert all(sp.simplify(c(a) - croot[a]) == 0 for a in croot)
    def zmono(m):
        terms = {(): sp.Integer(1)}
        for (nn, i) in m:
            new = {}
            for mono, x in terms.items():
                for j in range(8):
                    if g[j][i]:
                        mm = tuple(sorted(mono + ((nn, j),))); new[mm] = new.get(mm, 0) + x*g[j][i]
            terms = new
        return terms
    return g, cs, c, gv, zmono
if CASE == 'II':
    import exact_engine as E, Kfield as K
    Z1 = load_Z(HERE + '/../jumeau/j3_z_110.pkl'); Z2 = load_Z(HERE + '/../jumeau/j3_z_011.pkl')
    I86 = sp.eye(Z1.rows)
    res['z, zprime rational'] = all(x.is_Rational for x in list(Z1) + list(Z2))
    res['z^2 = 1'] = (Z1*Z1 == I86); res["z'^2 = 1"] = (Z2*Z2 == I86); res["z z' = z' z"] = (Z1*Z2 == Z2*Z1)
    print(res, flush=True)
    def make_zvecK(Zm):
        g, cs, c, gv, zmono = lift_data(Zm)
        def zvecK(v):
            out = {}
            for (lam, m), x in v.items():
                cl = c(lam); gl = gv(lam)
                for mm, y in zmono(m).items():
                    s = mpq(int(sp.fraction(cl*y)[0]), int(sp.fraction(cl*y)[1])); k = (gl, mm)
                    out[k] = K.add(out[k], K.scal(s, x)) if k in out else K.scal(s, x)
            return {k: a for k, a in out.items() if not K.iszero(a)}
        return zvecK, [str(x) for x in cs]
    zK, cz = make_zvecK(Z1); zpK, czp = make_zvecK(Z2)
    res['c(z)'] = cz; res["c(z')"] = czp
    sig = lambda v: {k: (a[0], -a[1], a[2], -a[3]) for k, a in v.items()}
    tau = lambda v: {k: (a[0], a[1], -a[2], -a[3]) for k, a in v.items()}
    PSI = {'sigma': lambda v: zK(sig(v)), 'tau': lambda v: zK(zpK(tau(v))), 'sigma tau': lambda v: zpK(sig(tau(v)))}
    S = E.setup(0); B2 = E.Blocks(S['T'], 2); span = S['U2span']
    piv = {}
    for v in span: piv.setdefault(B2.wt(v), []).append(v)
    basis = {mu: [vs[i] for i in E.select_pivots(vs, B2.pos[mu], len(B2.blk[mu]))] for mu, vs in piv.items()}
    res['dim_K U_2'] = sum(len(b) for b in basis.values())
    for name, psi in PSI.items():
        bymu = {}
        for v in span:
            w = psi(v); bymu.setdefault(B2.wt(w), []).append(w)
        ok_all = True
        for mu, ws in bymu.items():
            ok, bad = E.solve_block(basis.get(mu, []), ws, B2.pos[mu], len(B2.blk[mu]))
            if not ok or bad: ok_all = False
        res['psi_%s(U_2) c U_2 (exact)' % name] = ok_all; print(name, ok_all, flush=True)
    # explicit Q-form of U_2: traces T(a u) = sum_g psi_g(a u)
    BASK = [K.ONE, K.R, K.Q, K.RQ]
    def trace(v):
        out = dict(v)
        for psi in PSI.values(): E.kv_add(out, psi(v))
        return E.kv_clean(out)
    from flint import nmod_mat
    import numpy as np
    qrank = 0; krank = 0
    for mu, bs in basis.items():
        n_ = len(B2.blk[mu]); tr = [trace(E.kv_scale(a, b)) for b in bs for a in BASK]
        tr = [t for t in tr if t]
        # Q-rank: Weil coordinates modulo p (a lower bound for the Q-rank; = 4 dim_K would be impossible)
        cols = [E.weil_column(t, B2.pos[mu], n_)[0] for t in tr]
        qr = nmod_mat(len(cols), 4*n_, [x % E.P for col in cols for x in col], E.P).rank() if cols else 0
        kr = nmod_mat(len(tr), n_, [int(x) for t in tr for x in E.modp_vec(t, B2.pos[mu], n_)], E.P).rank() if tr else 0
        qrank += qr; krank += kr
        # fixed points: psi_g(T) = T exactly
        for t in tr[:2]:
            for name, psi in PSI.items():
                d = dict(psi(t)); E.kv_add(d, t, K.scal(mpq(-1), K.ONE)); assert not E.kv_clean(d), ('trace not fixed', name)
    res['Q-rank of the traces (mod p, lower bound)'] = qrank; res['K-rank of the traces'] = krank
    print('Q-rank', qrank, 'K-rank', krank, flush=True)
    # the idempotents f_1, f_2 of C_2 are fixed
    Bb = pickle.load(open('griess_B_II_basis.pkl', 'rb'))
    f1 = ['80/171', '-13/114 + sqrt(13)/114', '-25/228 - 5*sqrt(13)/228', '7/152 - sqrt(13)/152', '4*sqrt(13)/19']
    f2 = ['91/171', '13/114 - sqrt(13)/114', '5*sqrt(13)/228 + 25/228', '-7/152 + sqrt(13)/152', '-4*sqrt(13)/19']
    for nm, f in (('f_1', f1), ('f_2', f2)):
        vec = {}
        for m, t in enumerate(f):
            c = K.from_sympy(sp.sympify(t))
            for k, x in Bb[m].items(): vec[k] = K.add(vec.get(k, K.Z), K.scal(x, c))
        vec = E.kv_clean(vec)
        for name, psi in PSI.items():
            d = dict(psi(vec)); E.kv_add(d, vec, K.scal(mpq(-1), K.ONE))
            res['psi_%s(%s) = %s' % (name, nm, nm)] = not E.kv_clean(d)
    print(res, flush=True)
else:
    FULL = {n: {k: sp.sympify(c) for k, c in d.items()} for n, d in pickle.load(open(HERE + '/../cas-I/r1_riesen_transported.pkl', 'rb')).items()}
    def zvec_of(bits):
        Zm = load_Z(HERE + '/../jumeau/j6_zI_%s.pkl' % bits)
        g, cs, c, gv, zmono = lift_data(Zm)
        def zv(v):
            out = {}
            for (lam, m), x in v.items():
                for mm, y in zmono(m).items():
                    k = (gv(lam), mm); out[k] = sp.expand(out.get(k, 0) + c(lam)*y*x)
            return {k: w for k, w in out.items() if w != 0}
        return Zm, zv, [str(x) for x in cs]
    bits = {'z_1': '001', 'z_1 z_A': '101'}       # eps = 2*bit - 1:  001 -> (-1,-1,1) = z_1,  101 -> (1,-1,1) = z_1 z_A
    Zs = {nm: zvec_of(b) for nm, b in bits.items()}
    Mall = [load_Z(HERE + '/../jumeau/j6_zI_%s.pkl' % ''.join(b)) for b in itertools.product('01', repeat=3)]
    I86 = sp.eye(Mall[0].rows)
    res['Z_I rational'] = all(x.is_Rational for M in Mall for x in M)
    res['Z_I elementary abelian'] = all(M*M == I86 for M in Mall) and all(A*B == B*A for A in Mall for B in Mall)
    res['c(z_1)'] = Zs['z_1'][2]; res['c(z_1 z_A)'] = Zs['z_1 z_A'][2]
    R13, R3, II = sp.symbols('R13 R3 II')
    def galois(v, s13, s3, si):         # sqrt13 -> s13*sqrt13, sqrt3 -> s3*sqrt3, i -> si*i
        out = {}
        for k, x in v.items():
            e = sp.expand(x).xreplace({sp.sqrt(39): R13*R3, sp.sqrt(13): R13, sp.sqrt(3): R3, sp.I: II})
            out[k] = sp.expand(e.xreplace({R13: s13*sp.sqrt(13), R3: s3*sp.sqrt(3), II: si*sp.I}))
        return out
    def prop(a, b):                      # a = lambda b ?
        ks = set(a) | set(b)
        if set(a) != set(b): return False
        k0 = next(iter(b)); lam = sp.radsimp(a[k0]/b[k0])
        return all(sp.simplify(sp.expand(a[k] - lam*b[k])) == 0 for k in ks), lam
    def eq(a, b): return set(a) == set(b) and all(sp.simplify(sp.expand(a[k] - b[k])) == 0 for k in a)
    z1, z1zA = Zs['z_1'][1], Zs['z_1 z_A'][1]
    tau_ = lambda v: galois(v, 1, 1, -1)        # complex conjugation: sqrt-3 -> -sqrt-3
    sig_ = lambda v: galois(v, -1, 1, 1)        # an extension of sigma: sqrt13 -> -sqrt13, sqrt3, i fixed
    cand = {'tau': {'z_1': lambda v: z1(tau_(v))},
            'sigma': {'z_1 z_A (label (1,-1,1))': lambda v: z1zA(sig_(v)), 'z_A (label (-1,1,1)) = z_1 (z_1 z_A)': lambda v: z1(z1zA(sig_(v)))}}
    for g, opts in cand.items():
        for nm, psi in opts.items():
            out = {'e': eq(psi(FULL['e']), FULL['e']), 'f': eq(psi(FULL['f']), FULL['f'])}
            pv = prop(psi(FULL['v']), FULL['v'])
            out['v'] = (pv if isinstance(pv, bool) else (pv[0], str(pv[1])))
            res['psi_%s with a = %s' % (g, nm)] = out; print(g, nm, out, flush=True)
json.dump(res, open(HERE + '/d1_rational_forms_%s.json' % CASE, 'w'), indent=1, default=str)
print('done')
