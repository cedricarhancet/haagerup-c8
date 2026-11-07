#!/usr/bin/env python3 -I
"""A1 — exact integer/rational verification of the rank-two lattice facts used in Theorem A.

Ambient space: Q^2 with basis (p, b), Gram diag(6, 26)  (DJX: p^2 = 6, b^2 = 26, (p,b) = 0).
Vectors are pairs of Fractions (x, y) meaning x p + y b.
  lambda = (p+b)/2, s = (p+b)/4, t = (13p - 3b)/4.
  L_II = Z p + Z lambda            (Arh26 case (II), Gram [[6,3],[3,8]])
  L_I  = Z s + Z (s+t)/2           (Arh26 case (I),  Gram [[2,1],[1,20]])
  N    = Z lambda + Z 2p           (common index-two sublattice)
Checks: Gram matrices, L_I ∩ L_II = N, indices, the characters of sigma and tau, the discriminant form of N,
its 2-part = Z(Vec_{Z/2}) with a = s+N, b = p+N, ab fermion, the canonical isometry psi_N : L_II'/L_II -> L_I'/L_I and
its value u in {4, 17, 22, 35} on the generators (1,1) of [Arh26, Lemma 12.1].
"""
from fractions import Fraction as F
from itertools import product

G = [[F(6), F(0)], [F(0), F(26)]]

def ip(u, v):
    return u[0]*v[0]*G[0][0] + u[1]*v[1]*G[1][1]

def add(u, v): return (u[0]+v[0], u[1]+v[1])
def sc(c, u): return (c*u[0], c*u[1])

p = (F(1), F(0)); b = (F(0), F(1))
lam = sc(F(1, 2), add(p, b))
s = sc(F(1, 4), add(p, b))
t = sc(F(1, 4), add(sc(13, p), sc(-3, b)))
u_I = sc(F(1, 2), add(s, t))          # (s+t)/2 = (7p-b)/4

def gram(B):
    return [[ip(x, y) for y in B] for x in B]

out = []
def rep(line): print(line); out.append(line)

rep(f"lambda^2 = {ip(lam,lam)}, s^2 = {ip(s,s)}, t^2 = {ip(t,t)}, (s,t) = {ip(s,t)}, (s,p) = {ip(s,p)}, (p,lambda) = {ip(p,lam)}")
assert ip(lam, lam) == 8 and ip(s, s) == 2 and ip(t, t) == 78 and ip(s, t) == 0
rep(f"(s+t)/2 = {u_I}  (expected (7/4, -1/4) = (7p-b)/4), norm {ip(u_I,u_I)}")
assert u_I == (F(7, 4), F(-1, 4)) and ip(u_I, u_I) == 20

B_II = [p, lam]; B_I = [s, u_I]; B_N = [lam, sc(2, p)]
rep(f"Gram L_II = {gram(B_II)}  (expected [[6,3],[3,8]])")
rep(f"Gram L_I  = {gram(B_I)}  (expected [[2,1],[1,20]])")
rep(f"Gram N    = {gram(B_N)}  (expected [[8,6],[6,24]])")
assert gram(B_II) == [[6, 3], [3, 8]] and gram(B_I) == [[2, 1], [1, 20]] and gram(B_N) == [[8, 6], [6, 24]]

def det2(M): return M[0][0]*M[1][1] - M[0][1]*M[1][0]
rep(f"det L_I = {det2(gram(B_I))}, det L_II = {det2(gram(B_II))}, det N = {det2(gram(B_N))}")
assert det2(gram(B_I)) == 39 and det2(gram(B_II)) == 39 and det2(gram(B_N)) == 156

# --- membership tests via coordinates in a basis (solve 2x2 over Q) ---
def coords(v, B):
    """coordinates of v in basis B (vectors in Q^2 wrt (p,b))"""
    (a, c), (d, e) = B[0], B[1]
    det = a*e - c*d
    x = (v[0]*e - v[1]*d)/det
    y = (a*v[1] - c*v[0])/det
    return (x, y)
def in_lattice(v, B):
    x, y = coords(v, B)
    return x.denominator == 1 and y.denominator == 1

# L_I ∩ L_II = N : enumerate a box of L_I and test
box = range(-12, 13)
inter = set()
for m, n in product(box, repeat=2):
    v = add(sc(m, B_I[0]), sc(n, B_I[1]))
    if in_lattice(v, B_II):
        inter.add(v)
Nset = set(add(sc(m, B_N[0]), sc(n, B_N[1])) for m, n in product(range(-30, 31), repeat=2))
assert inter <= Nset and all(in_lattice(v, B_I) and in_lattice(v, B_II) for v in B_N)
rep(f"L_I ∩ L_II = N on a box of {len(inter)} vectors: OK (N ⊂ L_I ∩ L_II and L_I ∩ L_II ⊂ N)")
# indices
rep(f"[L_I : N] = {int(abs(det2(gram(B_N))/det2(gram(B_I)))**0.5)}, [L_II : N] = {int(abs(det2(gram(B_N))/det2(gram(B_II)))**0.5)}  (det ratio 156/39 = 4 = 2^2)")

# --- characters of sigma and tau ---
# sigma on L_II : alpha -> (-1)^{(alpha,lambda)} ;  tau on L_I : alpha -> (-1)^{(alpha,2p)} = (-1)^{(alpha,3s+t)/2}
assert sc(2, p) == sc(F(1, 2), add(sc(3, s), t))
for m, n in product(range(-6, 7), repeat=2):
    v = add(sc(m, B_II[0]), sc(n, B_II[1]))
    e = ip(v, lam); assert e.denominator == 1
    assert (e % 2 == 0) == in_lattice(v, B_N)
    w = add(sc(m, B_I[0]), sc(n, B_I[1]))
    f = ip(w, sc(2, p)); assert f.denominator == 1
    assert (f % 2 == 0) == in_lattice(w, B_N)
rep("ker((-1)^{(.,lambda)} on L_II) = N and ker((-1)^{(.,2p)} on L_I) = N : OK")
rep(f"(alpha, lambda) mod 2 on the basis (p, lambda) of L_II: {ip(p,lam)%2}, {ip(lam,lam)%2} ; (alpha,2p) on (s,(s+t)/2): {ip(s,sc(2,p))%2}, {ip(u_I,sc(2,p))%2}")

# --- dual lattices and discriminant forms ---
def dual_basis(B):
    Gm = gram(B); det = det2(Gm)
    inv = [[Gm[1][1]/det, -Gm[0][1]/det], [-Gm[1][0]/det, Gm[0][0]/det]]
    return [add(sc(inv[i][0], B[0]), sc(inv[i][1], B[1])) for i in range(2)]
Nd = dual_basis(B_N); Id = dual_basis(B_I); IId = dual_basis(B_II)
# N'/N as a finite abelian group: enumerate coset representatives
def disc_group(B):
    Bd = dual_basis(B)
    reps = {}
    for m, n in product(range(0, 200), repeat=2):
        v = add(sc(m, Bd[0]), sc(n, Bd[1]))
        key = tuple(x % 1 for x in coords(v, B))   # class mod lattice via fractional coordinates
        if key not in reps: reps[key] = v
    return reps
dN = disc_group(B_N); dI = disc_group(B_I); dII = disc_group(B_II)
rep(f"|N'/N| = {len(dN)}, |L_I'/L_I| = {len(dI)}, |L_II'/L_II| = {len(dII)}")
assert len(dN) == 156 and len(dI) == 39 and len(dII) == 39
def q(v): return (ip(v, v)/2) % 1
def order(v, B):
    for k in range(1, 400):
        if in_lattice(sc(k, v), B): return k
# 2-part of N'/N
two = [v for v in dN.values() if in_lattice(sc(2, v), B_N)]
rep("2-part of N'/N (q values): " + ", ".join(f"{coords(v,B_N)}: q={q(v)}" for v in two))
assert len(two) == 4
qs = sorted(q(v) for v in two)
assert qs == [0, 0, 0, F(1, 2)]
# identify a = s + N, b = p + N, ab = s + p + N
assert in_lattice(sc(2, s), B_N) and in_lattice(sc(2, p), B_N) and not in_lattice(s, B_N) and not in_lattice(p, B_N)
rep(f"q(s+N) = {q(s)}, q(p+N) = {q(p)}, q(s+p+N) = {q(add(s,p))}  -> a = s+N, b = p+N bosons, ab fermion: Z(Vec_Z2)")
assert q(s) == 0 and q(p) == 0 and q(add(s, p)) == F(1, 2)
# braiding bicharacter b(a,b) = q(a+b)-q(a)-q(b) = 1/2 : a, b not orthogonal
rep(f"bichar(a,b) = {(q(add(s,p)) - q(s) - q(p)) % 1} (expected 1/2)")
# 39-part
thirty9 = [v for v in dN.values() if in_lattice(sc(39, v), B_N)]
assert len(thirty9) == 39
rep(f"39-part of N'/N has {len(thirty9)} elements; it is cyclic: {any(order(v,B_N)==39 for v in thirty9)}")
# L_I'/N = e_I^perp : check that L_I' = {x in N' : (x, s) in Z} and L_II' = {x in N' : (x, p) in Z}
for v in dN.values():
    inLId = in_lattice(v, B_I) or any(in_lattice(add(v, sc(-1, w)), B_I) for w in [(F(0),F(0))])
# simpler: L_I' = dual of L_I ; test membership of v in L_I' by pairing with the basis of L_I
def in_dual(v, B): return all(ip(v, x).denominator == 1 for x in B)
cnt_I = sum(1 for v in dN.values() if in_dual(v, B_I)); cnt_II = sum(1 for v in dN.values() if in_dual(v, B_II))
rep(f"|L_I'/N| = {cnt_I} (expected 78), |L_II'/N| = {cnt_II} (expected 78)")
assert cnt_I == 78 and cnt_II == 78
assert all(in_dual(v, B_I) for v in thirty9) and all(in_dual(v, B_II) for v in thirty9)
rep("the 39-part of N'/N lies in L_I' ∩ L_II' : OK  (so L_I'/L_I ≅ 39-part ≅ L_II'/L_II canonically)")

# --- generators (1,1) of Arh26 Lemma 12.1 and the canonical isometry psi_N ---
# Arh26 case (II) coordinates (m',m'') -> (m'/6) p + (m''/26) b ; case (I): (m'/2) s + (m''/78) t
g_II = add(sc(F(1, 6), p), sc(F(1, 26), b))
g_I = add(sc(F(1, 2), s), sc(F(1, 78), t))
rep(f"q_II((1,1)) = {q(g_II)} (expected 4/39), q_I((1,1)) = {q(g_I)} (expected 10/39)")
assert q(g_II) == F(4, 39) and q(g_I) == F(10, 39)
assert in_dual(g_II, B_II) and order(g_II, B_II) == 39 and in_dual(g_I, B_I) and order(g_I, B_I) == 39
# psi_N: take y = g_II in N'; decompose y = y39 + y2 in N'/N; y39 in L_I'; find u with y39 ≡ u g_I mod L_I
def proj39(v):
    # 39-part component: 39-part = 4*(N'/N) ... use CRT: v = 40 v - 39 v ; 40v ≡ v mod ... simpler: y39 = 4*v*k where
    # multiplication by 4 kills the 2-part? no: 2-part has exponent 2, so 2v kills it; 39-part: mult by 39 kills it.
    # CRT: e = 39*k with 39k ≡ 1 mod 2 (k=1) -> 39v is in the 2-part; v - 39v = -38 v is in the 39-part? check orders.
    w = sc(F(-38), v)   # = v - 39 v
    assert in_lattice(sc(39, w), B_N)
    return w
y39 = proj39(g_II)
assert in_dual(y39, B_I)
u_found = None
for u in range(39):
    if in_lattice(add(y39, sc(-u, g_I)), B_I):
        u_found = u
rep(f"psi_N((1,1)_II) = u (1,1)_I with u = {u_found}  (Arh26 Lemma 12.1(3): u in {{4, 17, 22, 35}})")
assert u_found in (4, 17, 22, 35)
# isometry check for all elements
ok = True
for v in dII.values():
    w = proj39(v)
    if (q(w) - q(v)) % 1 != 0: ok = False
rep(f"psi_N preserves q on all 39 classes: {ok}")
assert ok
# Also: u^2 ≡ 16 mod 39 (Arh26 Lemma 12.1(3) : 10 u^2 ≡ 4)
rep(f"10*u^2 mod 39 = {10*u_found**2 % 39} (expected 4)")

# --- Weyl/orthogonal group check (Arh26 Lemma 12.1(4)): w' = diag(-1,1), w'' = diag(1,-1) in (p,b) coordinates act on L_II'/L_II by 14 and 25
def act(mat, v): return (mat[0]*v[0], mat[1]*v[1])
for name, mat, expect in [("w' (p -> -p)", (F(-1), F(1)), 14), ("w'' (b -> -b)", (F(1), F(-1)), 25)]:
    img = act(mat, g_II)
    k = [u for u in range(39) if in_lattice(add(img, sc(-u, g_II)), B_II)]
    rep(f"{name} acts on L_II'/L_II by {k} (expected {expect})")
    assert k == [expect]
# on L_I'/L_I: in (s,t) coordinates w' (s -> -s): acts by 1 ; w'' (t -> -t) by -1 = 38
for name, mat, expect in [("s -> -s", (F(-1), F(1)), 1), ("t -> -t", (F(1), F(-1)), 38)]:
    def act_st(v):
        # coordinates of v in (s,t): v = x s + y t with s,t orthogonal
        x = ip(v, s)/ip(s, s); y = ip(v, t)/ip(t, t)
        return add(sc(mat[0]*x, s), sc(mat[1]*y, t))
    img = act_st(g_I)
    k = [u for u in range(39) if in_lattice(add(img, sc(-u, g_I)), B_I)]
    rep(f"{name} acts on L_I'/L_I by {k} (expected {expect})")
    assert k == [expect]
rep("ALL CHECKS PASSED")
