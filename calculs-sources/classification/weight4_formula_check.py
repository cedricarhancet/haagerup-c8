"""Test of the weight-four one-point function formula (Axis B1).

Claim (Lemma, weight four): for h in V_1 with L(1)h = 0, h(0) semisimple, <h,h> = 1, and every module M,
   Z_M(h[-1]h[-1]h[-1]h, tau) = Tr_M h(0)^4 q^{L0-c/24} - (E2/2) Tr_M h(0)^2 q^{L0-c/24} + (E2^2/48) ch_M .
Proof idea: Heisenberg recursion Z(h[-1]b) = Z_{h(0)}(b) + Ehat_2 * Z(h[1]b) on Fock spaces, Ehat_2 = -E2/12.
Test: for lattice VOAs V_L (L = A2, E8, and the rank-one lattice sqrt(6)Z) and generic real h, the vector
   W_mu(tau) = eta^{-n} sum_{a in L+mu} ( <h,a>^4 - (E2/2) <h,h> <h,a>^2 + (E2^2/48) <h,h>^2 ) q^{<a,a>/2}
must satisfy W(-1/tau) = tau^4 S W(tau).  A wrong coefficient destroys this.  We also test that the
wrong coefficients (-E2/4, +E2^2/48) and (-E2/2, -E2^2/48) fail, to show the test is sensitive.
"""
import itertools
import mpmath as mp
mp.mp.dps = 30

def E2(tau):
    q = mp.exp(2j*mp.pi*tau)
    return 1 - 24*mp.nsum(lambda n: n*q**n/(1 - q**n), [1, mp.inf])
def eta(tau):
    q = mp.exp(2j*mp.pi*tau)
    return mp.exp(2j*mp.pi*tau/24)*mp.qp(q)

def lattice_W(G, B, cosets, h, tau, coef=(-0.5, 1/48.), R=7):
    """G: Gram matrix of basis B (rows = basis vectors in R^n); cosets: list of glue vectors (in R^n)"""
    n = len(B[0]); rk = len(B)
    e2 = E2(tau); q2 = mp.exp(mp.pi*1j*tau)   # q^{1/2}
    hh = sum(x*x for x in h)
    out = []
    for mu in cosets:
        s = mp.mpc(0)
        for c in itertools.product(range(-R, R + 1), repeat=len(B)):
            a = [mu[i] + sum(c[k]*B[k][i] for k in range(len(B))) for i in range(n)]
            nrm = sum(x*x for x in a)
            if nrm > 2*R*R/3:      # truncation (Im tau >= 0.8 makes this negligible)
                continue
            ha = sum(x*y for x, y in zip(h, a))
            s += (ha**4 + coef[0]*e2*hh*ha**2 + coef[1]*e2**2*hh**2)*q2**nrm
        out.append(s/eta(tau)**rk)
    return out

def check(name, B, cosets, Smat, h, coef):
    worst = 0
    for tau in (mp.mpc(0.11, 1.03), mp.mpc(-0.27, 0.97)):
        A = lattice_W(None, B, cosets, h, -1/tau, coef)
        Bv = lattice_W(None, B, cosets, h, tau, coef)
        for i in range(len(cosets)):
            rhs = tau**4*sum(Smat[i][j]*Bv[j] for j in range(len(cosets)))
            worst = max(worst, abs(A[i] - rhs)/max(abs(rhs), mp.mpf(10)**-20))
    print('%-6s coef=%s  worst relative error = %s' % (name, coef, mp.nstr(worst, 3)))
    return worst

# A2: basis in R^3 (sum zero), dual cosets
s = [1, -1, 0]; t = [0, 1, -1]
w1 = [mp.mpf(2)/3, -mp.mpf(1)/3, -mp.mpf(1)/3]
cos_A2 = [[0, 0, 0], w1, [-x for x in w1]]
om = mp.exp(2j*mp.pi/3)
S_A2 = [[1/mp.sqrt(3)*om**(-3*sum(x*y for x, y in zip(a, b))) for b in cos_A2] for a in cos_A2]
h_A2 = [mp.mpf('0.37'), mp.mpf('-0.81'), mp.mpf('0.44')]
# rank one: sqrt(6) Z, cosets k/sqrt(6), k = 0..5
r6 = mp.sqrt(6)
cos_1 = [[k/r6] for k in range(6)]
S_1 = [[mp.exp(-2j*mp.pi*j*k/6)/mp.sqrt(6) for k in range(6)] for j in range(6)]
res = []
for coef in ((-0.5, 1/48.), (-0.25, 1/48.), (-0.5, -1/48.)):
    res.append(check('A2', [s, t], cos_A2, S_A2, h_A2, coef))
    res.append(check('Z(s6)', [[r6]], cos_1, S_1, [mp.mpf('0.7')], coef))
# E8: Gram/basis via the even coordinate system D8+ ; small radius suffices for Im tau ~ 1
B_E8 = [[2, 0, 0, 0, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0], [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0], [0, 0, 0, 0, -1, 1, 0, 0], [0, 0, 0, 0, 0, -1, 1, 0],
        [mp.mpf(1)/2]*8]
def e8_vectors(maxnorm):
    out = []
    for v in itertools.product(range(-2, 3), repeat=8):
        if sum(v) % 2 == 0 and sum(x*x for x in v) <= maxnorm:
            out.append([mp.mpf(x) for x in v])
    for v in itertools.product([-3, -1, 1, 3], repeat=8):
        if (sum(v)//2) % 2 == 0 and sum(x*x for x in v) <= 4*maxnorm:
            out.append([mp.mpf(x)/2 for x in v])
    return out
V8 = e8_vectors(8)
def W_E8(h, tau, coef):
    e2 = E2(tau); q2 = mp.exp(mp.pi*1j*tau); hh = sum(x*x for x in h)
    s = sum((sum(x*y for x, y in zip(h, a))**4 + coef[0]*e2*hh*sum(x*y for x, y in zip(h, a))**2 + coef[1]*e2**2*hh**2)*q2**sum(x*x for x in a) for a in V8)
    return s/eta(tau)**8
h8 = [mp.mpf(x) for x in ('0.3', '-0.1', '0.25', '0.6', '-0.45', '0.05', '0.2', '-0.33')]
for coef in ((-0.5, 1/48.), (-0.25, 1/48.)):
    worst = 0
    for tau in (mp.mpc(0.05, 1.2),):
        A = W_E8(h8, -1/tau, coef); Bv = W_E8(h8, tau, coef)
        worst = max(worst, abs(A - tau**4*Bv)/abs(tau**4*Bv))
    print('E8     coef=%s  worst relative error = %s (norm <= 8 truncation)' % (coef, mp.nstr(worst, 3)))
