import sympy as sp
m = sp.symbols('m', positive=True, integer=True)
eps = sp.symbols('epsilon')
D = sp.sqrt(9*m**2+4)
d = (3*m+D)/2
# basic identities
print("d^2-1-3md:", sp.simplify(sp.expand(d**2-1-3*m*d)))
print("d-1/d-3m:", sp.simplify(d-1/d-3*m))
print("d+1/d-sqrtD:", sp.simplify(d+1/d-D))
dimF = 3+3*d**2
print("dimF-3(2+3md):", sp.simplify(sp.expand(dimF-3*(2+3*m*d))))
# center dims
dB = d**2; dX = 1+d**2; dZ = 2+3*m*d
print("dB=1+3md:", sp.simplify(sp.expand(dB-(1+3*m*d))))
print("dX=2+3md:", sp.simplify(sp.expand(dX-(2+3*m*d))))
# codegree of chi(rho)=-1/d
print("dimF/(3+3/d^2) - d^2:", sp.simplify(dimF/(3+3/d**2)-d**2))
# sum of squares
S = sp.symbols('S')
eq = sp.Eq(9*d**2*S, dimF**2 - 1 - dB**2 - dX**2 - 3*dZ**2)
sol = sp.solve(eq, S)[0]
print("sum n_j^2 =", sp.simplify(sp.radsimp(sp.expand(sol))))
rhs = 9*(1+d**2)**2-1-d**4-(1+d**2)**2-3*(1+d**2)**2
print("rhs - (4d^4+10d^2+4):", sp.simplify(sp.expand(rhs-(4*d**4+10*d**2+4))))
# first moment
T1 = sp.symbols('T1'); T2 = sp.symbols('T2')
s1 = sp.solve(sp.Eq(m*dB + m*dX + 3*d*T1, 0), T1)[0]
print("moment1 - (-(m^2+m sqrtD)/2):", sp.simplify(sp.radsimp(s1 + (m**2+m*D)/2)))
print("moment1 - (m^2 - m d):", sp.simplify(sp.radsimp(s1-(m**2-m*d))))
s2 = sp.solve(sp.Eq(m*dB + m*dX + 3*d*T2, eps*dimF), T2)[0]
print("moment2 - (-m^2+(2eps-m)sqrtD)/2:", sp.simplify(sp.radsimp(s2-(-m**2+(2*eps-m)*D)/2)))
# check numerically for m=1..5
for mm in range(1,6):
    print(mm, sp.N(sol.subs(m,mm)), 4*mm**2+2)
