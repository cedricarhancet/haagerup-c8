"""Weight-two (Schellekens-type) relations for a VOA U with Rep(U) = Z(H3) at c = 8.

For h in a Cartan subalgebra of U_1, Zhu's theorem gives that the vector
    W_z(tau) = Z_{M_z}(h[-1]h, tau) = Tr_{M_z} h(0)^2 q^{L_0 - 1/3} - (<h,h>/12) E_2(tau) X_z(tau)
is a holomorphic vector-valued modular form of weight two for the same representation rho as the
characters X = X^(gamma), with W_z = O(q^{mu_z}), mu_z the leading exponent of X_z (with pole terms).
Then F = W * E_10/Delta is weakly holomorphic of weight zero, hence F = Xi (a + J b) with b supported
on the polar coordinates, and F must vanish at tau = i and tau = rho = exp(2 pi i/3).
We compute this finite-dimensional space and the linear relations it imposes on the low-order
coefficients of W.  Output: relations among
   w0  = coeff of q^{-1/3} in W_0,    w0p = coeff of q^{2/3} in W_0,
   w_z = leading coefficient of W_z (z != 0), at q^{mu_z}.
Usage: python3 weight2_relations.py gamma   (gamma = 0 or 1)
"""
import sys
from fractions import Fraction as Fr
import mpmath as mp
import xi_exact as X

mp.mp.dps = 60
gamma = int(sys.argv[1]) if len(sys.argv) > 1 else 1
d, Lam, Xi, names = X.d, X.Lam, X.Xi, X.names
N = X.NMAX
D3 = names.index('d3')
P = [0, D3] if gamma == 1 else [0]           # coordinates where a pole q^{Lambda_z - 1} is allowed

# q-series of J = E4^3/Delta - 744 (J[m] = coeff of q^m, m >= -1) and g = Delta/E10 (g[k], k >= 0)
Dl = [Fr(x) for x in X.Dl]
invDl = X.series_inv(Dl, N)
E4c = [Fr(x) for x in X.E4c]
JJ = X.series_mul(E4c, invDl, N)              # = q * E4^3/Delta, so J_m = JJ[m+1]
J = {m: JJ[m+1] for m in range(-1, N - 1)}
J[0] -= 744
g = X.g

# exact coefficient functionals of W(a,b) = (Xi a + J Xi b) g ; return dict (z,n) -> vector over (a, b_P)
nvar = d + len(P)
def F_coeff(z, n):
    """coefficient of q^{Lam_z + n} in F_z = (Xi a + J Xi b)_z, as a list of length nvar"""
    v = [Fr(0)]*nvar
    if n >= -1:
        for j in range(d):
            v[j] += Xi[n][z][j]
    for k, p in enumerate(P):
        s = Fr(0)
        for m in range(-1, n + 2):
            if (n - m) >= -1 and (n - m) in Xi:
                s += J[m]*Xi[n - m][z][p]
        v[d + k] += s
    return v
def W_coeff(z, n):
    v = [Fr(0)]*nvar
    for k in range(1, n + 3):                 # g_k, k >= 1; F index n - k >= -2
        if g[k] == 0:
            continue
        fc = F_coeff(z, n - k)
        for i in range(nvar):
            v[i] += g[k]*fc[i]
    return v

funcs = {}
funcs['w0'] = W_coeff(0, -1)
funcs['w0p'] = W_coeff(0, 0)
for z in range(1, d):
    nlead = -1 if z in P else 0
    funcs['w_' + names[z]] = W_coeff(z, nlead)
labels = list(funcs)

# constraints F(i) = 0 and F(rho) = 0
def Xi_num(tau):
    return X.Xi_num(tau)
def J_num(tau):
    q = mp.exp(2j*mp.pi*tau)
    return sum(mp.mpf(J[m].numerator)/J[m].denominator*q**m for m in J)
rows = []
for tau in (mp.mpc(0, 1), mp.exp(2j*mp.pi/3)):
    Xm = Xi_num(tau); Jv = J_num(tau)
    for z in range(d):
        row = [Xm[z, j] for j in range(d)] + [Jv*Xm[z, p] for p in P]
        rows.append(row)
A = mp.matrix(rows)
# null space via SVD
U_, Sv, Vh = mp.svd_c(A)
sv = [abs(Sv[i]) for i in range(min(A.rows, A.cols))]
tol = mp.mpf(10)**(-35)
rank = sum(1 for s in sv if s > tol)
print('gamma =', gamma, ' unknowns =', nvar, ' rank of vanishing conditions =', rank,
      ' => dim M_2(rho, mu) =', nvar - rank)
ss = sorted(sv, reverse=True)
print('singular values of the vanishing conditions (descending):', [mp.nstr(x, 4) for x in ss])
print('gap: smallest retained =', mp.nstr(ss[rank-1], 5), ' largest dropped =', mp.nstr(ss[rank], 5) if rank < len(ss) else 'none')
U_, Sv, Vh = mp.svd_c(A, full_matrices=True)
null = [[Vh[i, j].conjugate() for j in range(nvar)] for i in range(rank, nvar)]   # rows of Vh beyond rank
# image of the functionals on the null space
M = mp.matrix(len(labels), len(null))
for r, lab in enumerate(labels):
    fv = funcs[lab]
    for c, vec in enumerate(null):
        M[r, c] = sum(mp.mpf(fv[i].numerator)/fv[i].denominator*vec[i] for i in range(nvar))
# left kernel of M = relations among the functionals
U2, S2, V2 = mp.svd_c(M, full_matrices=True)
s2 = [abs(S2[i]) for i in range(min(M.cols, M.rows))]
r2 = sum(1 for s in s2 if s > mp.mpf(10)**(-30))
print('number of functionals =', len(labels), ' rank of image =', r2, ' => relations =', len(labels) - r2)
# left kernel: columns r2.. of U2 (conjugated), since U2^H M = S V^H
rels = [[U2[j, i].conjugate() for j in range(len(labels))] for i in range(r2, len(labels))]
# normalize each relation: make it real and with a convenient pivot; then try to identify rationals
out = []
for rel in rels:
    # rotate phase so that the largest entry is real positive
    k = max(range(len(rel)), key=lambda j: abs(rel[j]))
    ph = rel[k]/abs(rel[k])
    rel = [x/ph for x in rel]
    out.append(rel)
# row-reduce the relation space to a canonical (rref) basis
import itertools
R = mp.matrix(out)
def rref(R):
    R = R.copy(); m, n = R.rows, R.cols; piv = []; r = 0
    for c in range(n):
        if r >= m: break
        p = max(range(r, m), key=lambda i: abs(R[i, c]))
        if abs(R[p, c]) < mp.mpf(10)**(-25): continue
        for j in range(n): R[r, j], R[p, j] = R[p, j], R[r, j]
        pv = R[r, c]
        for j in range(n): R[r, j] /= pv
        for i in range(m):
            if i != r and abs(R[i, c]) > 0:
                f = R[i, c]
                for j in range(n): R[i, j] -= f*R[r, j]
        piv.append(c); r += 1
    return R, piv
Rr, piv = rref(R) if len(out) else (None, [])
print('relations (reduced row echelon form over the labels):')
print('labels:', labels)
for i in range(len(out)):
    terms = []
    for j in range(len(labels)):
        x = Rr[i, j]
        if abs(x) > mp.mpf(10)**(-20):
            re = mp.re(x); im = mp.im(x)
            idf = mp.identify(re, tol=mp.mpf(10)**(-20)) if abs(im) < mp.mpf(10)**(-20) else None
            terms.append('%s*%s' % (idf if idf else mp.nstr(x, 15), labels[j]))
    print('  0 = ' + ' + '.join(terms))

# ---- exact rational recognition of the reduced relations (checked to 1e-40) ----
from fractions import Fraction
exact = []
for i in range(len(out)):
    row = []
    for j in range(len(labels)):
        x = mp.re(Rr[i, j])
        assert abs(mp.im(Rr[i, j])) < mp.mpf(10)**(-40)
        fr = Fraction(str(mp.nstr(x, 45))).limit_denominator(10**8)
        assert abs(x - mp.mpf(fr.numerator)/fr.denominator) < mp.mpf(10)**(-40), (x, fr)
        row.append(fr)
    exact.append(row)
print('exact relations:')
for row in exact:
    print('  0 = ' + ' + '.join('(%s)*%s' % (c, labels[j]) for j, c in enumerate(row) if c != 0))
import json
json.dump({'gamma': gamma, 'labels': labels, 'relations': [[str(c) for c in row] for row in exact]},
          open('weight2_relations_gamma%d.json' % gamma, 'w'), indent=1)

# ---- independent checks ----
# (1) every relation must vanish exactly on the explicit weight-two forms D(Xi e_p), p in P (D = q d/dq)
def DXi_coeff(z, n, p):
    # coefficient of q^{Lam_z + n} in q d/dq of Xi[:, p] component z
    return (Lam[z] + n)*Xi[n][z][p] if n in Xi else Fr(0)
for p in P:
    vals = {}
    vals['w0'] = DXi_coeff(0, -1, p)
    vals['w0p'] = DXi_coeff(0, 0, p)
    for z in range(1, d):
        nlead = -1 if z in P else 0
        vals['w_' + names[z]] = DXi_coeff(z, nlead, p)
    for row in exact:
        assert sum(c*vals[lab] for c, lab in zip(row, labels)) == 0
print('check 1: all relations vanish exactly on D(Xi e_p) for p in', [names[p] for p in P])
# (2) the null-space elements are genuine weight-two forms: W(-1/tau) = tau^2 S W(tau) at a sample point
S = X.S if hasattr(X, 'S') else None
from zh3_modular_data import zh3_S
Sm = zh3_S()
def W_num(vec, tau):
    a = [vec[i] for i in range(d)]; b = [vec[d + k] for k in range(len(P))]
    Xm = Xi_num(tau); Jv = J_num(tau)
    q = mp.exp(2j*mp.pi*tau)
    E4 = 1 + 240*sum(X.sigma(3, n)*q**n for n in range(1, 60))
    E6 = 1 - 504*sum(X.sigma(5, n)*q**n for n in range(1, 60))
    Dl_ = q*mp.fprod([(1 - q**n)**24 for n in range(1, 200)])
    Fv = [sum(Xm[z, j]*a[j] for j in range(d)) + Jv*sum(Xm[z, p]*b[k] for k, p in enumerate(P)) for z in range(d)]
    return [f*Dl_/(E4*E6) for f in Fv]
worst = 0
for vec in null:
    for tau in (mp.mpc(0.13, 1.07), mp.mpc(-0.31, 0.93)):
        A1 = W_num(vec, -1/tau); A2 = W_num(vec, tau)
        lhs = mp.matrix(A1); rhs = tau**2*(Sm*mp.matrix(A2))
        worst = max(worst, mp.mnorm(lhs - rhs, 1)/max(mp.mnorm(rhs, 1), mp.mpf(10)**-30))
print('check 2: worst relative error of W(-1/tau) = tau^2 S W(tau) on the null space:', mp.nstr(worst, 5))
