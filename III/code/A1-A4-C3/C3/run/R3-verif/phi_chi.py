# Step 4: antilinear involution phi_chi ; character chi ; action on solutions ; sign-character symmetries
import sys, pickle, time, itertools, random
sys.path.insert(0, '.')
from voa import *
from lat import *
from q4 import Q4
import sympy as sp
t0 = time.time()
S = pickle.load(open('sl2_4.pkl', 'rb'))
G = pickle.load(open('G.pkl', 'rb'))
E, F, H = S['E'], S['F'], S['H']
Gb, omCC = G['Gb'], G['omCC']
sols = pickle.load(open('solutions.pkl', 'rb'))['sols']
# --- lattice basis coordinates (simple roots of each E8 copy) ---
def c16(v): return coords16(v)
supp = []
for X in range(4):
    for a, c in E[X].items():
        k = -F[X][((), vneg(a[1]))]*c   # d = -k/c  => k = -d*c
        supp.append((a[1], k, c))
Mrows = sp.Matrix([c16(a) for a, k, c in supp])
print('rank of span of the 24 support roots:', Mrows.rank())
# Smith form of the sublattice M in L (coordinates in simple-root basis)
from sympy.matrices.normalforms import smith_normal_form
snf = smith_normal_form(Mrows, domain=sp.ZZ)
print('Smith invariants of L/M:', [snf[i, i] for i in range(16)])
# chi(alpha) = k_alpha / c_alpha^2  -> log2 values ; solve <lam, alpha> = log2 for lam in Q^16 (coordinates: chi(beta)=2^{l . coords(beta)})
import math
vals = []
for a, k, c in supp:
    r = Fr(k)/(c*c)
    lg = math.log2(r); assert 2**round(lg) == r
    vals.append(round(lg))
lv = sp.symbols('l0:16')
eqs = [sum(ci*li for ci, li in zip(c16(a), lv)) - v for (a, k, c), v in zip(supp, vals)]
sol = sp.solve(eqs, lv, dict=True)
print('chi exponent solution (chi(beta) = 2^{sum l_i b_i}, b = coords in simple roots):', sol)
lam = [sol[0][l] for l in lv]
assert all(x.is_rational for x in lam)
def chi_exp(beta):
    return sum(sp.Rational(l)*b for l, b in zip(lam, c16(beta)))
def chi(beta):
    e = chi_exp(beta)
    assert e.is_integer, ('non-integer exponent', beta, e)
    return Fr(2)**int(e)
# also express chi in terms of orthonormal coordinates: chi(beta) = 2^{<mu,beta>}
mu = sp.symbols('m0:16')
eqs2 = [sum(sp.Rational(a[i], 2)*mu[i] for i in range(16)) - v for (a, k, c), v in zip(supp, vals)]
sol2 = sp.solve(eqs2, mu, dict=True)
print('chi(beta) = 2^{<mu,beta>} with mu =', [sol2[0].get(m, m) for m in mu])
def phi(st, chif=chi):
    """phi_chi(c f(h) e^beta) = conj(c) chi(beta) (-1)^{deg f} f(h) e^{-beta}"""
    out = {}
    for (m, b), c in st.items():
        cc = c.conj() if isinstance(c, Q4) else c
        sgn = (-1)**len(m)
        out[(m, vneg(b))] = cc*chif(b)*sgn
    return out
def eqst(a, b):
    keys = set(a) | set(b)
    return all(Q4.lift(a.get(k, 0)) == Q4.lift(b.get(k, 0)) for k in keys)
ok = all(eqst(phi(E[X]), {k: -v for k, v in F[X].items()}) and eqst(phi(F[X]), {k: -v for k, v in E[X].items()}) and eqst(phi(H[X]), {k: -v for k, v in H[X].items()}) for X in range(4))
print('phi_chi(E_X) = -F_X, phi_chi(F_X) = -E_X, phi_chi(H_X) = -H_X for X=1..4:', ok)
# phi_1 is an antilinear automorphism: check phi1(a_(n)b) = phi1(a)_(n) phi1(b) on random basis states, n=-2..3
one = lambda b: Fr(1)
random.seed(7)
pool = ROOTS + NORM4
ok = True
ntest = 0
for trial in range(300):
    a = random.choice(pool); b = random.choice(pool)
    if trial % 3 == 0:
        # force interaction
        cands = [x for x in ROOTS if ip(x, a) < 0]
        b = random.choice(cands)
    sa = {((), a): Fr(1)} if trial % 4 else {(((1, random.randrange(16)),), a): Fr(1)}
    sb = {((), b): Fr(1)} if trial % 5 else {(((1, random.randrange(16)),), b): Fr(1)}
    if max(state_wt(k) for k in sa) > 3 or max(state_wt(k) for k in sb) > 3:
        continue
    for n in range(-1, 3):
        l = phi(mode_vec(sa, n, sb), one)
        r = mode_vec(phi(sa, one), n, phi(sb, one))
        ntest += 1
        if clean(l) != clean(r):
            ok = False
print('phi_1 antilinear automorphism test on %d random mode products: %s' % (ntest, ok))
# same with chi (chi is a character, so phi_chi also an automorphism)
ok2 = True
for trial in range(100):
    a = random.choice(ROOTS); b = random.choice([x for x in pool if ip(x, a) < 0])
    for n in range(-1, 3):
        if clean(phi(mode_vec(e_state(a), n, e_state(b)))) != clean(mode_vec(phi(e_state(a)), n, phi(e_state(b)))):
            ok2 = False
print('phi_chi automorphism test (random):', ok2)
def tostate(xc):
    st = {}
    for c, g in zip(xc, Gb):
        if c:
            for k, v in g.items():
                nv = st.get(k, Q4(0)) + c*v
                if nv: st[k] = nv
                else: st.pop(k, None)
    return st
U = [tostate(x) for x in sols]
print('phi_chi(omega_CC) = omega_CC:', eqst(phi(omCC), omCC))
print('phi_chi(omega_W) = omega_W:', eqst(phi(G['omW']), G['omW']), ' phi_chi(omega_AxA) = omega_AxA:', eqst(phi(G['omAA']), G['omAA']))
for k, u in enumerate(U):
    pu = phi(u)
    img = [l for l in range(4) if eqst(pu, U[l])]
    print('phi_chi(u%d) = u%s' % (k, img))
# phi_chi on G in G-coordinates (matrix): phi(g_i) = sum Phi_ij g_j
piv = G['pivots']
Phi = [[phi(g).get(p, Fr(0)) for p in piv] for g in Gb]
print('phi_chi maps G to G:', all(eqst(phi(g), tostate([Q4(x) for x in row])) for g, row in zip(Gb, Phi)))
# ---- sign characters t: L -> {+-1} trivial on all support roots: centralize sl2^4, act on G ----
# solve over GF(2): t(beta) = (-1)^{s . coords(beta)}, s in GF(2)^16, s . coords(alpha) = 0 for support roots
import itertools
A2 = [[x % 2 for x in c16(a)] for a, k, c in supp]
# nullspace over GF(2)
def gf2_null(rows, n):
    rows = [r[:] for r in rows]
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x+y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = [0]*n; v[f] = 1
        for i, c in enumerate(piv):
            v[c] = rows[i][f]
        basis.append(v)
    return basis
NB = gf2_null(A2, 16)
print('sign characters of L trivial on the support roots: 2^%d' % len(NB))
def tact(svec, st):
    return {k: v*(-1 if sum(x*y for x, y in zip(svec, c16(k[1]))) % 2 else 1) for k, v in st.items()}
perms = set()
for bits in itertools.product((0, 1), repeat=len(NB)):
    svec = [sum(b*v[i] for b, v in zip(bits, NB)) % 2 for i in range(16)]
    okc = all(eqst(tact(svec, x), x) for x in E + F + H)
    perm = tuple(next((l for l in range(4) if eqst(tact(svec, U[k]), U[l])), None) for k in range(4))
    perms.add(perm)
    print('  t =', bits, ' centralizes sl2^4:', okc, ' fixes omega_CC:', eqst(tact(svec, omCC), omCC), ' permutation of (u0..u3):', perm)
print('time %.1fs' % (time.time()-t0))
