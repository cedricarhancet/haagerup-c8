"""T12 -- independent check of E2 (bijection pi) and E3 (W_1, exclusion of so(32)).

Input : t12_table.json (produced by t12_table.py from Article I, read-only).
Output: t12_W1.json and stdout.  Exact integer arithmetic except the numerical
S-matrix check (floating point, tolerance 1e-9, used only as a control).
Usage: python3 -I t12_W1.py t12_table.json
"""
import sys, json, itertools, math, cmath
from fractions import Fraction as F
from collections import Counter

T = json.load(open(sys.argv[1]))
labs = ["0", "b", "a", "c0", "c1", "c2", "d1", "d2", "d3", "d4", "d5", "d6"]
h = {z: F(T[z]["h"]) for z in labs}
R = {z: {tuple(map(int, k.split(","))): m for k, m in T[z]["R"].items()} for z in labs}

# ---------------------------------------------------------------- Z(H3) data
d = (3 + math.sqrt(13)) / 2
FP = {"0": 1, "b": 3*d+1, "a": 3*d+2, "c0": 3*d+2, "c1": 3*d+2, "c2": 3*d+2}
for l in range(1, 7):
    FP[f"d{l}"] = 3*d
theta = {"0": F(0), "b": F(0), "a": F(0), "c0": F(0), "c1": F(1, 3), "c2": F(2, 3)}
for l in range(1, 7):
    theta[f"d{l}"] = F(6*l*l, 13) % 1           # Evans-Gannon (2.5a)
# (H_J) + pseudo-unitarity: exp(2 pi i h_z) = theta_z
assert all((h[z] - theta[z]) % 1 == 0 for z in labs), "twists do not match h"
print("twists exp(2 pi i h_z) = theta_z (EG (2.5a)) : OK")
print("sum FPdim^2 =", sum(FP[z]**2 for z in labs), " (9d+6)^2 =", (9*d+6)**2)

# EG S-matrix (2.5b), order labs, S = (1/3) * M
x = (13 - 3*math.sqrt(13)) / 26
y = 3 / math.sqrt(13)
c = lambda j: -2*y*math.cos(2*math.pi*j/13)
M = [[0.0]*12 for _ in range(12)]
row0 = [x, 1-x, 1, 1, 1, 1] + [y]*6
rowb = [1-x, x, 1, 1, 1, 1] + [-y]*6
M[0] = row0; M[1] = rowb
blk = [[2, -1, -1, -1], [-1, 2, -1, -1], [-1, -1, -1, 2], [-1, -1, 2, -1]]
for i in range(4):
    M[2+i] = [1, 1] + blk[i] + [0]*6
for l in range(1, 7):
    M[5+l] = [y, -y, 0, 0, 0, 0] + [c((l*m) % 13) for m in range(1, 7)]
S = [[M[i][j]/3 for j in range(12)] for i in range(12)]
# controls: symmetric, unitary (real orthogonal), S_{0z}/S_{00} = FPdim
assert all(abs(S[i][j]-S[j][i]) < 1e-12 for i in range(12) for j in range(12))
for i in range(12):
    for j in range(12):
        v = sum(S[i][k]*S[k][j] for k in range(12))
        assert abs(v - (1 if i == j else 0)) < 1e-9
for i, z in enumerate(labs):
    assert abs(S[0][i]/S[0][0] - FP[z]) < 1e-9
print("EG S-matrix: real symmetric, S^2 = 1 (all simples self-dual), S_0z/S_00 = FPdim : OK")
# Verlinde: integrality/nonnegativity of fusion coefficients (control)
for i in range(12):
    for j in range(12):
        for k in range(12):
            N = sum(S[i][m]*S[j][m]*S[k][m]/S[0][m] for m in range(12))
            assert abs(N - round(N)) < 1e-8 and round(N) >= 0
print("Verlinde fusion coefficients nonnegative integers : OK")

# all permutations pi with pi(0)=0, FP preserved, theta_pi(i) = -theta_i, S_{pi i, pi j} = conj S_ij = S_ij
cands = {z: [w for w in labs if abs(FP[w]-FP[z]) < 1e-9 and (theta[w] + theta[z]) % 1 == 0] for z in labs}
print("candidates for pi(z):", cands)
sols = []
def rec(k, cur, used):
    if k == 12:
        sols.append(dict(cur)); return
    z = labs[k]
    for w in cands[z]:
        if w in used:
            continue
        cur[z] = w
        okk = all(abs(S[labs.index(cur[u])][labs.index(w)] - S[labs.index(u)][k]) < 1e-9 for u in labs[:k+1])
        if okk:
            rec(k+1, cur, used | {w})
        del cur[z]
rec(0, {}, frozenset())
print("number of permutations of the modular data inverting T and fixing S:", len(sols))
for s in sols:
    print("   ", {z: s[z] for z in labs if s[z] != z})
pi = sols[0]

# ---------------------------------------------------------------- W = sum M_i x M_pi(i)
print("\nh_i + h_pi(i):", {z: str(h[z] + h[pi[z]]) for z in labs})
assert all((h[z] + h[pi[z]]).denominator == 1 for z in labs)
assert all(h[z] > 0 for z in labs if z != "0")
w1pairs = [z for z in labs if h[z] + h[pi[z]] == 1]
print("pairs (i,pi(i)) of weight 1:", [(z, pi[z]) for z in w1pairs])

def dimV(t):
    p = 1
    for a in t:
        p *= a + 1
    return p

W1 = Counter()
# (S x S)_1 = S_1 x 1 + 1 x S_1 ; S_1 = g' + g'' = V(2,0)+V(0,2)
for t in [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)]:
    W1[t] += 1
for z in w1pairs:
    for (a, b), m in R[z].items():
        for (c_, d_), n in R[pi[z]].items():
            W1[(a, b, c_, d_)] += m*n
dimW1 = sum(dimV(t)*m for t, m in W1.items())
print("dim W_1 =", dimW1)
print("trivial constituents of W_1:", W1.get((0, 0, 0, 0), 0))
print("number of constituents:", sum(W1.values()), " distinct types:", len(W1))

def ind1(n):          # Dynkin index of V(n) for sl2, normalised ind V(1) = 1
    return n*(n+1)*(n+2)//6

def index_vec(mod):
    out = [0, 0, 0, 0]
    for t, m in mod.items():
        for k in range(4):
            oth = dimV(t) // (t[k]+1)
            out[k] += m*oth*ind1(t[k])
    return out

iw = index_vec(W1)
print("Dynkin indices of W_1 per factor:", iw, " = 60 *", [v/60 for v in iw])

# centraliser of S_1 x 1 in W_1 (types (0,0,*,*))
cent = {t: m for t, m in W1.items() if t[0] == 0 and t[1] == 0}
print("centraliser of S_1 x 1 in W_1:", cent, " dim", sum(dimV(t)*m for t, m in cent.items()))

# single-factor restrictions
for k in range(4):
    res = Counter()
    for t, m in W1.items():
        res[t[k]] += m*dimV(t)//(t[k]+1)
    print(f"W_1 restricted to factor {k}:", dict(sorted(res.items())))

# ---------------------------------------------------------------- so(32) test
def weights(t):
    return list(itertools.product(*[range(-a, a+1, 2) for a in t]))

def char(mod):
    C = Counter()
    for t, m in mod.items():
        for w in weights(t):
            C[w] += m
    return C

charW1 = char(W1)

def wedge2_char(U):
    ws = []
    for t, m in U.items():
        ws += weights(t)*m
    C = Counter()
    for i in range(len(ws)):
        for j in range(i+1, len(ws)):
            C[tuple(ws[i][k]+ws[j][k] for k in range(4))] += 1
    return C

target = [6, 26, 6, 26]          # 2 * levels (index of the vector rep of so(N) is 2)
types = [t for t in itertools.product(range(3), range(5), range(3), range(5)) if t != (0, 0, 0, 0)]
types = [t for t in types if dimV(t) <= 32 and all(dimV(t)//(t[k]+1)*ind1(t[k]) <= target[k] for k in range(4))]
types.sort(key=lambda t: -dimV(t))
tind = {t: [dimV(t)//(t[k]+1)*ind1(t[k]) for k in range(4)] for t in types}

def enum(multfree):
    res = []
    def rec(i, budget, dimleft, cur):
        if budget == [0, 0, 0, 0]:
            U = dict(cur)
            if dimleft > 0:
                U[(0, 0, 0, 0)] = dimleft
            res.append(U); return
        if i == len(types):
            return
        t = types[i]
        dt = dimV(t); it = tind[t]
        maxm = 1 if multfree else 32
        if multfree and sum(t) % 2 == 1:
            maxm = 0             # symplectic type: multiplicity must be 0 when centraliser is 0
        m = 0
        while m <= maxm:
            nb = [budget[k]-m*it[k] for k in range(4)]
            if min(nb) < 0 or dimleft - m*dt < 0:
                break
            if multfree or sum(t) % 2 == 0 or m % 2 == 0:
                if m:
                    cur.append((t, m))
                rec(i+1, nb, dimleft-m*dt, cur)
                if m:
                    cur.pop()
            m += 1
    rec(0, list(target), 32, [])
    return res

allU = enum(False)
orth = [U for U in allU if all(sum(t) % 2 == 0 or m % 2 == 0 for t, m in U.items())]
mf = [U for U in orth if all(m == 1 for t, m in U.items() if t != (0, 0, 0, 0)) and U.get((0, 0, 0, 0), 0) <= 1
      and all(sum(t) % 2 == 0 for t in U)]
print("\n32-dim reps U of sl2^4 with so(32)-indices (3,13,3,13):", len(allU))
print("   of which orthogonal (symplectic types with even multiplicity):", len(orth))
print("   of which with zero centraliser in so(U) (mult-free, orthogonal types, <=1 trivial):", len(mf))
hits = [U for U in orth if wedge2_char(U) == charW1]
print("orthogonal U with Lambda^2 U = W_1 :", len(hits))
# finer: which invariants already fail
def summary(U):
    return {"index": index_vec(U), "dim": sum(dimV(t)*m for t, m in U.items())}
for U in mf:
    print("   mult-free candidate:", sorted(U.items()), summary(U))

json.dump({"pi": pi, "all_pi": sols, "W1": {",".join(map(str, t)): m for t, m in sorted(W1.items())},
           "dimW1": dimW1, "index_W1": iw, "n_U_all": len(allU), "n_U_orth": len(orth),
           "n_U_multfree": len(mf), "n_hits": len(hits),
           "multfree_U": [[[list(t), m] for t, m in sorted(U.items())] for U in mf]},
          open("t12_W1.json", "w"), indent=1)
