from fractions import Fraction as Fr
import itertools
# Bourbaki E8 simple roots (1-based e_1..e_8 -> index 0..7)
h = Fr(1,2)
al = [[h,-h,-h,-h,-h,-h,-h,h],[1,1,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],
      [0,0,-1,1,0,0,0,0],[0,0,0,-1,1,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,-1,1,0]]
al = [[Fr(x) for x in r] for r in al]
dot = lambda a,b: sum(x*y for x,y in zip(a,b))
C = [[dot(a,b) for b in al] for a in al]
import sympy as sp
Ci = sp.Matrix(C).inv()
print('inverse Cartan diag:', [Ci[i,i] for i in range(8)])
print('min entry:', min(Ci), ' all integer:', all(x.is_integer for x in Ci))
# fundamental weights in coordinates
A = sp.Matrix(al)  # rows simple roots
W = Ci*A  # rows = fundamental weights (since <w_i, a_j> = delta)
for i in range(8): print('varpi_%d ='%(i+1), list(W.row(i)), ' norm', dot(list(W.row(i)),list(W.row(i))))
theta = list(W.row(7))
print('theta(varpi_7) =', dot(theta, list(W.row(6))))
print('<alpha_8, e7+e8> =', dot(al[7],[0,0,0,0,0,0,1,1]), ' <alpha_i,e7+e8> =', [dot(a,[0,0,0,0,0,0,1,1]) for a in al])
# check smallest norms of nonzero dominant elements sum a_i varpi_i with a_i in {0,1,2}
vals = sorted(set(sum(a[i]*a[k]*Ci[i,k] for i in range(8) for k in range(8)) for a in itertools.product(range(3),repeat=8) if any(a)))
print('smallest norms of dominant integral coweights:', vals[:6])
