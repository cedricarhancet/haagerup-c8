"""Test of the one-point function formula in weight six (and four) on lattice vertex operator algebras.
Heisenberg recursion: Z(h[-1]^{n+1} 1) = h(0) Z(h[-1]^n 1) + n e <h,h> Z(h[-1]^{n-1} 1), e = -E2/12, gives
  n = 4 : h^4 + 6 e h^2 + 3 e^2      n = 6 : h^6 + 15 e h^4 + 45 e^2 h^2 + 15 e^3   (<h,h> = 1)
Check: W(-1/tau) = tau^n S W(tau) for V_{A2} and V_{sqrt6 Z}; a perturbed coefficient must fail.
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
def poly(n, lam, e, pert=0):
    if n == 4:
        return lam**4 + 6*e*lam**2 + 3*e**2
    return lam**6 + (15 + pert)*e*lam**4 + 45*e**2*lam**2 + 15*e**3
def W(B, cosets, h, tau, n, pert=0, R=7):
    e = -E2(tau)/12; q2 = mp.exp(mp.pi*1j*tau); hh = sum(x*x for x in h)
    out = []
    for mu in cosets:
        s = 0
        for c in itertools.product(range(-R, R + 1), repeat=len(B)):
            a = [mu[i] + sum(c[k]*B[k][i] for k in range(len(B))) for i in range(len(B[0]))]
            nrm = sum(x*x for x in a)
            if nrm > 30: continue
            lam = sum(x*y for x, y in zip(h, a))/mp.sqrt(hh)
            s += poly(n, lam, e, pert)*q2**nrm
        out.append(s/eta(tau)**len(B))
    return out
s_, t_ = [1, -1, 0], [0, 1, -1]
w1 = [mp.mpf(2)/3, -mp.mpf(1)/3, -mp.mpf(1)/3]
cos = [[0, 0, 0], w1, [-x for x in w1]]
om = mp.exp(2j*mp.pi/3)
S = [[om**(-3*sum(x*y for x, y in zip(a, b)))/mp.sqrt(3) for b in cos] for a in cos]
h = [mp.mpf('0.37'), mp.mpf('-0.81'), mp.mpf('0.44')]
r6 = mp.sqrt(6); cos1 = [[k/r6] for k in range(6)]
S1 = [[mp.exp(-2j*mp.pi*j*k/6)/r6 for k in range(6)] for j in range(6)]
for n, pert in ((4, 0), (6, 0), (6, 1)):
    for name, B, cs, Sm, hv in (('A2', [s_, t_], cos, S, h), ('Z(s6)', [[r6]], cos1, S1, [mp.mpf('0.7')])):
        worst = 0
        for tau in (mp.mpc(0.11, 1.03), mp.mpc(-0.27, 0.97)):
            A = W(B, cs, hv, -1/tau, n, pert); Bv = W(B, cs, hv, tau, n, pert)
            for i in range(len(cs)):
                rhs = tau**n*sum(Sm[i][j]*Bv[j] for j in range(len(cs)))
                worst = max(worst, abs(A[i] - rhs)/abs(rhs))
        print('weight %d %s %-6s worst relative error %s' % (n, '(perturbed)' if pert else '           ', name, mp.nstr(worst, 3)))
