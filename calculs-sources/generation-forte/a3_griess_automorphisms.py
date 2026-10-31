"""Axis A3 (exact): automorphisms of the Griess algebra B = Com_{V'}(A)_2 (case II, dim 5) fixing e.
An automorphism g permutes the 15 two-dimensional subalgebras containing e and preserves their central charges; the
lines with a charge pair occurring once (subalgebras 2, 3, 6, 7, 11, 12, 13) are fixed, with scaling 1 (lambda != 0).
They span the hyperplane H = {x_4 = 0} of e^perp, so g|_H = id and, g being orthogonal, g = id or the reflection s_H.
We test whether s_H is an algebra automorphism and how it acts on the 15 subalgebras (in particular 9 <-> 10?)."""
import json
import sympy as sp
d = json.load(open('griess_B_II.json'))
n = len(d['gram'])
prod = {tuple(map(int, k.split(','))): [sp.sympify(c) for c in v] for k, v in d['prod'].items()}
gram = sp.Matrix([[sp.sympify(c) for c in r] for r in d['gram']])
def mult(a, b):
    out = sp.zeros(n, 1)
    for i in range(n):
        for j in range(n):
            if a[i] != 0 and b[j] != 0: out += a[i]*b[j]*sp.Matrix(prod[(i, j)])
    return out.applyfunc(lambda z: sp.radsimp(sp.expand(z)))
G = json.load(open('griess_subalgebras.json'))['II']['subalgebras']
X = [sp.Matrix([0] + [sp.sympify(c) for c in s['x']]) for s in G]
# span check of the fixed lines
Hlines = [X[i] for i in (2, 3, 6, 7, 11, 12, 13)]
print('rank of fixed lines:', sp.Matrix.hstack(*Hlines).rank())
# reflection through H (and e): fix e and basis 1..3, normal vector n = gram^{-1}-orthogonal complement
E = [sp.Matrix([1 if i == j else 0 for i in range(n)]) for j in range(n)]
Hb = [E[0], E[1], E[2], E[3]]
# normal: solve <n, h> = 0 for h in Hb, n = E4 + sum c_i Hb_i
cs = sp.symbols('c0:4')
nv = E[4] + sum((c*h for c, h in zip(cs, Hb)), sp.zeros(n, 1))
sol = sp.solve([(h.T*gram*nv)[0] for h in Hb], cs, dict=True)[0]
nv = nv.subs(sol)
def s_H(v):
    return (v - 2*((nv.T*gram*v)[0]/(nv.T*gram*nv)[0])*nv).applyfunc(sp.radsimp)
# is s_H an algebra automorphism?
ok = all(sp.simplify(mult(s_H(E[i]), s_H(E[j])) - s_H(mult(E[i], E[j]))) == sp.zeros(n, 1) for i in range(n) for j in range(i, n))
print('s_H is an automorphism of B:', ok)
def same_line(u, v):
    M = sp.Matrix.hstack(u, v).applyfunc(sp.radsimp)
    return all(sp.radsimp(sp.expand(M[i, 0]*M[j, 1] - M[j, 0]*M[i, 1])) == 0 for i in range(n) for j in range(n))
perm = {}
for i in range(len(X)):
    img = s_H(X[i])
    # image line modulo e: subalgebra span{e, x}; compare e-perp parts (X already in e-perp)
    js = [j for j in range(len(X)) if same_line(img, X[j])]
    perm[i] = js
print('action of s_H on the 15 subalgebras:', perm)
json.dump({'automorphism': bool(ok), 'perm': {str(k): v for k, v in perm.items()}}, open('a3_griess_automorphisms.json', 'w'), indent=1)
