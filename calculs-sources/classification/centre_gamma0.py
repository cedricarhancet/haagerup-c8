"""gamma = 0: the weight-one Lie algebra s = U_1 is semisimple (exact).
For a real central element h with <h,h> = 1 and p2(h) = 0 one has u = 1 and T = T4 = 0.  The weight-four
relations give Q5 = Tr_{R_{d5}} h(0)^4 = (12 T + 2 T4 - 27 u^2 - 90 u + 105)/26 = -6/13 < 0: impossible.
If dim z >= 2, ker(p2|z) != 0 since p2(z) lies in the centre of p2(s), which is 0 or a line."""
import sympy as sp
import weight4_solve as W
eqs, tsol, R4 = W.weight4_system(0)
sol = sp.solve(eqs, list(W.Q[:5]), dict=True)[0]
Q5 = sp.factor(sol[W.Q[4]])
print('Q5 =', Q5)
val = Q5.subs({W.u: 1, W.T: 0, W.T4: 0})
print('central h with p2(h) = 0 (u = 1, T = T4 = 0): Q5 =', val)
assert val < 0
# also: for central h, Q5 >= 0 forces u <= u_max
umax = max(sp.solve(sp.Eq(Q5.subs({W.T: 0, W.T4: 0}), 0), W.u))
print('central elements in general: Q5 >= 0 forces u <=', umax, '=', sp.N(umax, 6), '(and u >= 5/9 from weight two)')
