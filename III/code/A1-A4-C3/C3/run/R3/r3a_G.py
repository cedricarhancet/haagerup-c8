"""R3a step 2: G = Com_W(A x A)_2 (exact over Q), Sugawara vectors, omega_CC = omega_W - omega_{AxA}, and the
determination of omega_{C x 1} in G by LINEAR conditions: o(omega_{C x 1}) acts on each sl2^4-isotypic component
(a,b,c,d) of W_1 lying in the sector R_z x R_{pi z} by h_z - h_A(a,b) (0 on (S x S)_1); we impose the traces of o(u)
on the highest weight spaces of all components (scalars on multiplicity-one components).
Output: G.pkl (basis of G, omega vectors)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, time, json
from fractions import Fraction as Fr
from gmpy2 import mpq
from flint import fmpq_mat, fmpq
import e8e8_voa as V
import common as CM
t0 = time.time()
G = CM.cond_kernel(2, (0, 0, 0, 0), CM.primary_ops())
print('dim G = dim Com_W(A x A)_2 =', len(G), '(%.0fs)' % (time.time() - t0), flush=True)
# Sugawara vectors
def sug(nm):
    E, H, F = CM.SL[nm]; k = CM.LEV[nm]
    v = {}
    V.vadd(v, V.mode_vec(E, -1, F)); V.vadd(v, V.mode_vec(F, -1, E)); V.vadd(v, V.mode_vec(H, -1, H), mpq(1, 2))
    return V.scale(v, mpq(1, 2*(k + 2)))
SUG = {nm: sug(nm) for nm in CM.NAMES}
for nm in CM.NAMES:
    c2 = V.mode_vec(SUG[nm], 3, SUG[nm]).get(V.VAC, 0)
    print('  Sugawara', nm, ': c/2 =', c2, ' expected', mpq(3*CM.LEV[nm], 2*(CM.LEV[nm] + 2)))
omW = V.omega()
omAA = {}
for nm in CM.NAMES: V.vadd(omAA, SUG[nm])
omCC = V.vadd(dict(omW), omAA, -1)
cc = CM.span_coords(G, omCC)
print('omega_CC in G:', cc is not None)
print('c(omega_CC)/2 =', V.mode_vec(omCC, 3, omCC).get(V.VAC), '(expected 36/10 = 18/5)')
# linear conditions for omega_{C x 1}
comp = CM.w1_sectors()
rows = []; rhs = []; rows1 = []; rhs1 = []
info = []
for kappa, secs in sorted(comp.items()):
    hw = CM.hw_W1(kappa)
    assert len(hw) == len(secs), (kappa, len(hw), secs)
    tr = [sum(CM.zero_mode_matrix(g, kappa)[i][i] for i in range(len(hw))) for g in G]
    target = sum(CM.omega_C1_eig(kappa, s) for s in secs)
    target1 = sum(CM.omega_1C_eig(kappa, s) for s in secs)
    rows.append(tr); rhs.append(target); rhs1.append(target1)
    info.append((kappa, secs, str(target), str(target1)))
print('trace conditions built (%.0fs)' % (time.time() - t0), flush=True)
def solve(rows, rhs):
    n = len(G)
    M = fmpq_mat(len(rows), n + 1)
    for i, (r, b) in enumerate(zip(rows, rhs)):
        for j in range(n): M[i, j] = r[j]
        M[i, n] = fmpq(b.numerator, b.denominator)
    Rr, rk = M.rref()
    A = fmpq_mat(len(rows), n)
    for i in range(len(rows)):
        for j in range(n): A[i, j] = M[i, j]
    rkA = A.rref()[1]
    sol = [fmpq(0)]*n; consistent = True
    for i in range(rk):
        j = next(j for j in range(n + 1) if Rr[i, j] != 0)
        if j == n: consistent = False; break
        sol[j] = Rr[i, n]
    return consistent, rkA, sol
ok, rk, sol = solve(rows, rhs)
print('omega_{C x 1}: system consistent', ok, '; rank', rk, 'of', len(G))
ok1, rk1, sol1 = solve(rows, rhs1)
print('omega_{1 x C}: system consistent', ok1, '; rank', rk1, 'of', len(G))
omC1 = V.clean(CM.combine([None]*0, []) or {})
def comb(coefs):
    v = {}
    for g, c in zip(G, coefs):
        if c != 0: V.vadd(v, g, CM.to_mq(c))
    return v
omC1 = comb(sol); om1C = comb(sol1)
print('omega_{C x 1} + omega_{1 x C} = omega_CC :', V.clean(V.vadd(dict(omC1), om1C)) == V.clean(omCC))
for nm, u, c in (('omega_{C x 1}', omC1, mpq(9, 5)), ('omega_{1 x C}', om1C, mpq(9, 5))):
    uu = V.mode_vec(u, 1, u)
    print(nm, ': u_(1) u = 2 u :', V.clean(uu) == V.scale(u, 2), '; u_(3) u =', V.mode_vec(u, 3, u).get(V.VAC), '(c/2 = 9/5)',
          '; support', len(u), flush=True)
print('omega_{C x 1}_(1) omega_{1 x C} = 0 :', V.clean(V.mode_vec(omC1, 1, om1C)) == {})
pickle.dump({'G': G, 'SUG': SUG, 'omCC': omCC, 'omC1': omC1, 'om1C': om1C, 'info': info, 'unique': rk == len(G)},
            open(os.path.join(CM.HERE, 'G.pkl'), 'wb'))
print('time %.0fs' % (time.time() - t0))
