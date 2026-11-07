#!/usr/bin/env python3 -I
"""A1b -- exact verification of the additional lattice facts used in the identification U = S^natural
(Section 7.3 of Article 3).  Pure Fractions, no floating point.

Part 1 (coordinates (p, b) of h, Gram diag(6, 26)):
  (s, L_II') = Z/2 and (p, L_I') = Z/2 ; kernels of the characters x -> (s,x), x -> (p,x) on L_II, L_I are N ;
  B := {x in L_II' : (s,x) in Z} and B' := {x in L_I' : (p,x) in Z} satisfy B = B' = preimage of the 39-part P of N'/N ;
  L_II = N u (p+N), L_I = N u (s+N) ; 2s = lambda in N, 2p in N ; (s+p)^2 = 11 ;
  psi_N(class of gamma in L_II'/L_II) = class of gamma in L_I'/L_I for all gamma in B (definition of psi_N) ;
  isotropy: no nonzero isotropic element in L_II'/L_II nor in L_I'/L_I (q(k g) = 4k^2/39, 10k^2/39 not in Z for 39 !| k),
  so L_I and L_II have no proper even overlattice in h ;
  the sublattice {gamma in L_I' : (p,gamma) in Z} has index 2 in L_I' and does not contain s.
Part 2 (root coordinates alpha_1..alpha_8 of DJX (4.1), chain 1-2-3-4-5, 6 attached to 3, A2 = {7,8}):
  p, b, lambda in L = E6 + A2 ; s not in L ; L_0 = {x in L : (lambda, x) in 2Z} = {x : sum x_i even} ;
  lambda, 2p, (t-s)/2 in L_0 ; (s, L_0) in Z ; L_H = L_0 + Zs even ; (p, s) = 3/2, (p, L_0) in Z ;
  the vectors of h cap L (resp. h cap L_H, h cap L_0) on a box are exactly L_II (resp. L_I, N) ;
  E := L_0 cap h^perp = L cap h^perp = L_H cap h^perp on a box.
"""
from fractions import Fraction as F
from itertools import product

out = []
def rep(line):
    print(line); out.append(line)

# ------------------------------------------------------------------ Part 1: h-coordinates
G2 = [[F(6), F(0)], [F(0), F(26)]]
def ip2(u, v): return u[0]*v[0]*G2[0][0] + u[1]*v[1]*G2[1][1]
def add2(u, v): return (u[0]+v[0], u[1]+v[1])
def sc2(c, u): return (c*u[0], c*u[1])
p = (F(1), F(0)); b = (F(0), F(1))
lam = sc2(F(1, 2), add2(p, b)); s = sc2(F(1, 4), add2(p, b)); t = sc2(F(1, 4), add2(sc2(13, p), sc2(-3, b)))
u_I = sc2(F(1, 2), add2(s, t))
B_II = [p, lam]; B_I = [s, u_I]; B_N = [lam, sc2(2, p)]

def coords(v, B):
    (a, c), (d, e) = B[0], B[1]
    det = a*e - c*d
    return ((v[0]*e - v[1]*d)/det, (a*v[1] - c*v[0])/det)
def in_lat(v, B):
    x, y = coords(v, B); return x.denominator == 1 and y.denominator == 1
def dual_basis(B):
    # vectors d_i with (d_i, B_j) = delta_ij
    Gm = [[ip2(x, y) for y in B] for x in B]
    det = Gm[0][0]*Gm[1][1] - Gm[0][1]*Gm[1][0]
    inv = [[Gm[1][1]/det, -Gm[0][1]/det], [-Gm[1][0]/det, Gm[0][0]/det]]
    return [add2(sc2(inv[i][0], B[0]), sc2(inv[i][1], B[1])) for i in range(2)]
def box(B, R=6):
    return [add2(sc2(F(i), B[0]), sc2(F(j), B[1])) for i in range(-R, R+1) for j in range(-R, R+1)]

D_II = dual_basis(B_II); D_I = dual_basis(B_I); D_N = dual_basis(B_N)
gen_II = add2(sc2(F(1, 6), p), sc2(F(1, 26), b))      # (1,1) of Arh26 case (II)
gen_I = add2(sc2(F(1, 2), s), sc2(F(1, 78), t))       # (1,1) of Arh26 case (I)
assert in_lat(gen_II, D_II) and in_lat(gen_I, D_I)

# (s, L_II') = Z/2, (p, L_I') = Z/2
vals_s = {(ip2(s, x) * 2) % 2 for x in box(D_II)}
vals_p = {(ip2(p, x) * 2) % 2 for x in box(D_I)}
rep(f"(s, L_II') takes the values {sorted(vals_s)} mod 2 after doubling (expected {{0,1}}: (s,L_II') = Z/2)")
rep(f"(p, L_I')  takes the values {sorted(vals_p)} mod 2 after doubling (expected {{0,1}}: (p,L_I') = Z/2)")
assert vals_s == {0, 1} and vals_p == {0, 1}
assert all((2*ip2(s, x)).denominator == 1 for x in box(D_II)) and all((2*ip2(p, x)).denominator == 1 for x in box(D_I))
rep(f"(s,(1,1)_II) = {ip2(s, gen_II)} (expected 1/2), (p,(1,1)_I) = {ip2(p, gen_I)} (expected 1), (s,p) = {ip2(s,p)}, (p,(s+t)/2) = {ip2(p,u_I)}")
assert ip2(s, gen_II) == F(1, 2) and ip2(p, gen_I) == 1 and ip2(s, p) == F(3, 2) and ip2(p, u_I) == F(21, 2)

# kernels on L_II and L_I
kerS = [x for x in box(B_II) if ip2(s, x).denominator == 1]
kerP = [x for x in box(B_I) if ip2(p, x).denominator == 1]
assert all(in_lat(x, B_N) for x in kerS) and all(in_lat(x, B_N) for x in kerP)
assert all(ip2(s, x).denominator == 1 for x in box(B_N)) and all(ip2(p, x).denominator == 1 for x in box(B_N))
rep("ker(x -> (s,x) mod Z) on L_II = N and ker(x -> (p,x) mod Z) on L_I = N : OK")
# L_II = N u (p+N), L_I = N u (s+N)
assert all(in_lat(x, B_N) or in_lat(add2(x, sc2(-1, p)), B_N) for x in box(B_II))
assert all(in_lat(x, B_N) or in_lat(add2(x, sc2(-1, s)), B_N) for x in box(B_I))
assert not in_lat(p, B_N) and not in_lat(s, B_N) and in_lat(sc2(2, s), B_N) and in_lat(sc2(2, p), B_N)
rep("L_II = N u (p+N), L_I = N u (s+N), 2s = lambda in N, 2p in N, s,p not in N : OK")
rep(f"(s+p)^2 = {ip2(add2(s,p), add2(s,p))} (expected 11), so the weight of e^(s+p) is 11/2")
assert ip2(add2(s, p), add2(s, p)) == 11

# N'/N, its 39-part P, and the sets B, B'
def cls_N(x):  # canonical representative of x + N in a fundamental domain, as coords mod 1 in basis B_N
    c = coords(x, B_N); return (c[0] % 1, c[1] % 1)
Ndual_box = box(D_N, 80)
classes = {}
for x in Ndual_box:
    classes.setdefault(cls_N(x), x)
assert len(classes) == 156
P39 = {k for k, x in classes.items() if cls_N(sc2(39, x)) == (0, 0)}
rep(f"|N'/N| = {len(classes)}, 39-part P has {len(P39)} elements (expected 39)")
assert len(P39) == 39
Bset = {cls_N(x) for x in box(D_II, 40) if ip2(s, x).denominator == 1}
Bpset = {cls_N(x) for x in box(D_I, 40) if ip2(p, x).denominator == 1}
rep(f"B/N = {{x in L_II' : (s,x) in Z}}/N has {len(Bset)} classes, B'/N = {{x in L_I' : (p,x) in Z}}/N has {len(Bpset)} classes")
assert Bset == P39 and Bpset == P39
rep("B/N = B'/N = P (the 39-part of N'/N) : OK")
# s not in B' (index two), p not in B
assert cls_N(s) not in P39 and cls_N(p) not in P39
# psi_N : class in L_II'/L_II -> class in L_I'/L_I, through P ; compare with the known value u = 22 on (1,1)
def cls_lat(x, B):
    c = coords(x, B); return (c[0] % 1, c[1] % 1)
# the map P -> L_II'/L_II and P -> L_I'/L_I are bijections
imII = {cls_lat(classes[k], B_II) for k in P39}; imI = {cls_lat(classes[k], B_I) for k in P39}
assert len(imII) == 39 and len(imI) == 39
rep("P -> L_II'/L_II and P -> L_I'/L_I are bijections : OK")
# psi_N(gen_II) : find gamma in P with gamma = gen_II mod L_II, then its class mod L_I equals u*gen_I mod L_I
gam = next(classes[k] for k in P39 if cls_lat(classes[k], B_II) == cls_lat(gen_II, B_II))
u = next(k for k in range(39) if cls_lat(sc2(F(k), gen_I), B_I) == cls_lat(gam, B_I))
rep(f"psi_N((1,1)_II) = {u} (1,1)_I  (expected 22, as in a1_lattices.py)")
assert u == 22
# q-values and isotropy
def q(x): return (ip2(x, x)/2) % 1
rep(f"q((1,1)_II) = {q(gen_II)}, q((1,1)_I) = {q(gen_I)}")
iso_II = [k for k in range(1, 39) if q(sc2(F(k), gen_II)) == 0]
iso_I = [k for k in range(1, 39) if q(sc2(F(k), gen_I)) == 0]
rep(f"isotropic nonzero classes of L_II'/L_II: {iso_II}; of L_I'/L_I: {iso_I} (expected none): L_I, L_II have no proper even overlattice in h")
assert iso_II == [] and iso_I == []

# ------------------------------------------------------------------ Part 2: root coordinates of DJX
EDGES = [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (7, 8)]
G8 = [[F(0)]*8 for _ in range(8)]
for i in range(8): G8[i][i] = F(2)
for (i, j) in EDGES: G8[i-1][j-1] = G8[j-1][i-1] = F(-1)
def ip8(x, y): return sum(x[i]*G8[i][j]*y[j] for i in range(8) for j in range(8))
def add8(x, y): return tuple(x[i]+y[i] for i in range(8))
def sc8(c, x): return tuple(F(c)*x[i] for i in range(8))
P8 = tuple(F(v) for v in (1, 2, 3, 2, 1, 0, 0, 0)); B8 = tuple(F(v) for v in (3, 4, 3, 4, 3, 2, 2, 2))
LAM8 = sc8(F(1, 2), add8(P8, B8)); S8 = sc8(F(1, 4), add8(P8, B8)); T8 = sc8(F(1, 4), add8(sc8(13, P8), sc8(-3, B8)))
rep(f"root coordinates: p = {tuple(int(v) for v in P8)}, b = {tuple(int(v) for v in B8)}, lambda = {tuple(int(v) for v in LAM8)}")
rep(f"s = {tuple(str(v) for v in S8)}, t = {tuple(str(v) for v in T8)}")
assert ip8(P8, P8) == 6 and ip8(B8, B8) == 26 and ip8(P8, B8) == 0 and ip8(S8, S8) == 2 and ip8(T8, T8) == 78 and ip8(S8, T8) == 0
def integral(x): return all(v.denominator == 1 for v in x)
assert integral(P8) and integral(B8) and integral(LAM8) and not integral(S8) and not integral(T8)
rep("p, b, lambda in L ; s, t not in L : OK")
# L_0 = {x in L : (lambda,x) in 2Z} = {x : sum x_i even}: check (lambda, alpha_i) = +-1 for all i
Glam = [ip8(LAM8, tuple(F(1) if k == i else F(0) for k in range(8))) for i in range(8)]
rep(f"(lambda, alpha_i) = {[int(v) for v in Glam]} (all odd, so L_0 = {{x : sum x_i even}})")
assert all(v % 2 == 1 for v in Glam)
def in_L0(x): return integral(x) and (sum(x) % 2 == 0)
TmS2 = sc8(F(1, 2), add8(T8, sc8(-1, S8)))
rep(f"(t-s)/2 = {tuple(int(v) for v in TmS2)}, in L_0: {in_L0(TmS2)}; lambda in L_0: {in_L0(LAM8)}; 2p in L_0: {in_L0(sc8(2, P8))}; p in L_0: {in_L0(P8)}")
assert in_L0(TmS2) and in_L0(LAM8) and in_L0(sc8(2, P8)) and not in_L0(P8)
# (s, L_0) in Z, (p, L_0) in Z, (p,s) = 3/2 ; L_H = L_0 + Zs even
basis8 = [tuple(F(1) if k == i else F(0) for k in range(8)) for i in range(8)]
L0_gens = [add8(basis8[i], basis8[0]) for i in range(1, 8)] + [sc8(2, basis8[0])]   # generate L_0 (index 2)
assert all(ip8(S8, g).denominator == 1 for g in L0_gens) and all(ip8(P8, g).denominator == 1 for g in L0_gens)
assert ip8(P8, S8) == F(3, 2)
rep("(s, L_0) in Z, (p, L_0) in Z, (p,s) = 3/2 : OK ; hence sigma = exp(2 pi i s(0)) and tau' = exp(2 pi i p(0)) are trivial on V_{L_0}")
# L_H even: norms of g + eps s for generators and s itself even
assert all((ip8(add8(g, S8), add8(g, S8))) % 2 == 0 for g in L0_gens) and ip8(S8, S8) % 2 == 0
rep("L_H = L_0 + Zs is even : OK")
# h cap L, h cap L_H, h cap L_0 on a box: elements x p + y b ; (xp+yb, p) = 6x and (xp+yb, b) = 26y lie in Z/2 for
# vectors of L_H (which pair with s to Z/2 ... all vectors of L_H have half-integral root coordinates), so x in Z/12, y in Z/52
hits_L, hits_LH, hits_L0 = [], [], []
for i in range(-5*12, 5*12+1):
    for j in range(-5*52, 5*52+1):
        x = F(i, 12); y = F(j, 52)
        v8 = add8(sc8(x, P8), sc8(y, B8))
        if integral(v8):
            hits_L.append((x, y))
            if in_L0(v8): hits_L0.append((x, y))
        if in_L0(v8) or in_L0(add8(v8, sc8(-1, S8))):
            hits_LH.append((x, y))
assert all(in_lat(v, B_II) for v in hits_L) and all(in_lat(v, B_I) for v in hits_LH) and all(in_lat(v, B_N) for v in hits_L0)
# conversely the box vectors of L_II, L_I, N (coordinates bounded by 5) are hit
for v in box(B_II, 2): assert v in hits_L
for v in box(B_I, 2): assert v in hits_LH
for v in box(B_N, 2): assert v in hits_L0
rep(f"h cap L = L_II, h cap L_H = L_I, h cap L_0 = N on a box ({len(hits_L)}, {len(hits_LH)}, {len(hits_L0)} vectors) : OK")
# E = L cap h^perp = L_H cap h^perp = L_0 cap h^perp : enumerate integral vectors with |x_i| <= 2 orthogonal to p and b
E_L = [x for x in product(range(-2, 3), repeat=8) if ip8(tuple(F(v) for v in x), P8) == 0 and ip8(tuple(F(v) for v in x), B8) == 0]
assert all(sum(x) % 2 == 0 for x in E_L)
rep(f"integral vectors of L with |x_i|<=2 orthogonal to h: {len(E_L)}, all in L_0 : OK (indeed (lambda,x) = 0 for x in h^perp, so L cap h^perp = L_0 cap h^perp)")
# L_H cap h^perp subset L_0: an element x + s (x in L_0) orthogonal to p has (p, x+s) = (p,x) + 3/2 not in Z, impossible
rep("L_H cap h^perp = L_0 cap h^perp since (p, x+s) = (p,x) + 3/2 is not an integer for x in L_0 : OK (no enumeration needed)")
def rank(rows):
    M = [list(r) for r in rows]; rk = 0; cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((r for r in range(rk, len(M)) if M[r][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][c] != 0:
                f = M[r][c]/M[rk][c]; M[r] = [M[r][k] - f*M[rk][k] for k in range(cols)]
        rk += 1
    return rk
# rank of E = 6 : E is the kernel in L of the two linear forms x -> (p,x), x -> (b,x), which are independent
rows = [[ip8(P8, basis8[i]) for i in range(8)], [ip8(B8, basis8[i]) for i in range(8)]]
rep(f"(p, alpha_i) = {[int(v) for v in rows[0]]}, (b, alpha_i) = {[int(v) for v in rows[1]]}, rank {rank(rows)} (expected 2, so rank E = 6)")
assert rank(rows) == 2
# (lambda, x) = 0 for x in E, so E is contained in L_0 (consistent with the enumeration above)
rep("ALL CHECKS PASSED (a1b_extra)")
