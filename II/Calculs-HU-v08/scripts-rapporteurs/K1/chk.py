import sympy as sp
# 1) L0 trick on coefficients: P_s Y(a,z)b = a_(k)b z^{-k-1}, s = ha+hb-k-1 ; (L_{-1}a)_(m) = -m a_(m-1)
ha,hb,k,z,lam=sp.symbols('h_a h_b k z lambda',positive=True)
s=ha+hb-k-1
lhs=s*z**(-k-1)                                  # L0 P_s Y(a,z)b  (coefficient of a_(k)b)
m=k+1                                            # (L_{-1}a)_(m) b has weight s  <=> m=k+1
rhs=ha*z**(-k-1)+z*(-m)*z**(-m-1)+hb*z**(-k-1)   # Y(L0 a + z L_{-1}a,z)b + Y(a,z)L0 b, s-component
print("L0 trick:",sp.simplify(lhs-rhs))
# 2) rescaling: P_s Y(lam^{ha+hb} a⊗b, lam z) = lam^s P_s Y(a,z)b
print("rescaling:",sp.simplify(lam**(ha+hb)*(lam*z)**(-k-1)-lam**s*z**(-k-1)))
# 3) R is swap-symmetric and (z,z) in R iff 2^{-1/2}<|z|<1  (R: 1>|z|>|conj(w)^{-1}-z|, 1>|w|>0)
import random,cmath
def inR(z,w): return 1>abs(z)>abs(1/w.conjugate()-z)>0 and 1>abs(w)>0
bad=0
for _ in range(200000):
    z=cmath.rect(random.random(),random.uniform(-3.2,3.2)); w=cmath.rect(random.random(),random.uniform(-3.2,3.2))
    if inR(z,w)!=inR(w,z): bad+=1
print("swap asym count:",bad)
print([inR(complex(r),complex(r)) for r in (0.70,0.71,0.72,0.99)])
