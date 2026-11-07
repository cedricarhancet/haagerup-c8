"""T12 -- explicit realisations of W_1 = e8 + e8 as an sl2^4-module ("next step").

Configuration T12 (the enumeration below is exhaustive for this configuration, as checked by the
referees R3P and F1, see ../../COMPUTATIONS.md and ../scripts-rapporteurs/F1/enum_sl2.py;
the only other configuration, "option B", is excluded by ../R3/t12b_optionB.py):
in each copy e8^(k), one level-3 factor g'_k is the
2A1-subalgebra (Dynkin index 2), whose centraliser is so(13) (type B6), with
     e8 | sl2 + so(13) = (V2,1) + (V2,13) + (V0,78) + (V1,64)          (*)
(checked below on so(16)+128).  The three other factors act on C^13 through a 13-dim
orthogonal representation U'_k, i.e. a homomorphism sl2^3 -> so(13); its Dynkin indices in e8
equal those in so(13) = ind(U'_k)/2.  Then e8^(k) restricted to sl2^4 is
     X_k = (V2 x 1) + V2 x U'_k + Lambda^2 U'_k + V1 x Spin(U'_k).
We search U'_1, U'_2 with X_1 + X_2 = W_1 (characters in Z^4).  Every solution is an explicit
injective homomorphism sl2^4 -> e8 + e8 realising the branching of W_1.
Factor order everywhere: (g'_1, g''_1, g'_2, g''_2).
Usage: python3 -I t12_e8split.py t12_W1.json
"""
import sys, json, itertools
from collections import Counter
from fractions import Fraction as F

D = json.load(open(sys.argv[1]))
W1 = {tuple(map(int, k.split(","))): m for k, m in D["W1"].items()}
ind1 = lambda n: n*(n+1)*(n+2)//6
def dimV(t):
    p = 1
    for a in t: p *= a+1
    return p
def weights(t):
    return list(itertools.product(*[range(-a, a+1, 2) for a in t]))
def char(mod):
    C = Counter()
    for t, m in mod.items():
        for w in weights(t): C[w] += m
    return C
charW1 = char(W1)

# ---- check (*) via e8 = so(16) + 128 under so(4) + so(12), diagonal sl2 in so(4) = sl2 + sl2
# so(16)|(sl2a,sl2b,so12) = (2,0,1)+(0,2,1)+(1,1,12)+(0,0,66); 128 = (1,0,32)+(0,1,32')
diag = Counter()
for (a, b, d), m in {(2, 0, 1): 1, (0, 2, 1): 1, (1, 1, 12): 1, (0, 0, 66): 1, (1, 0, 32): 1, (0, 1, 32): 1}.items():
    # V(a) x V(b) restricted to the diagonal = sum V(a+b-2j)
    for j in range(min(a, b)+1):
        diag[a+b-2*j] += m*d
print("e8 | diagonal sl2 (2A1):", dict(diag), " -> V2:", diag[2], "V1:", diag[1], "V0:", diag[0])
assert diag[2] == 14 and diag[1] == 64 and diag[0] == 78

# ---- 13-dim orthogonal reps of sl2^3 with prescribed so(13)-indices
def orth_reps(target, maxw=(4, 4, 4)):
    types = [t for t in itertools.product(*[range(m+1) for m in maxw]) if any(t)]
    tind = {t: [dimV(t)//(t[i]+1)*ind1(t[i]) for i in range(3)] for t in types}
    types = [t for t in types if dimV(t) <= 13 and all(tind[t][i] <= 2*target[i] for i in range(3))]
    out = []
    def rec(i, budget, dl, cur):
        if budget == [0, 0, 0]:
            U = dict(cur)
            if dl: U[(0, 0, 0)] = dl
            out.append(U); return
        if i == len(types): return
        t = types[i]; m = 0
        while True:
            nb = [budget[k]-m*tind[t][k] for k in range(3)]
            if min(nb) < 0 or dl-m*dimV(t) < 0: break
            if sum(t) % 2 == 0 or m % 2 == 0:
                if m: cur.append((t, m))
                rec(i+1, nb, dl-m*dimV(t), cur)
                if m: cur.pop()
            m += 1
    rec(0, [2*x for x in target], 13, [])
    return out

def spin_weights(U):
    ws = []
    for t, m in U.items(): ws += weights(t)*m
    ws = Counter(ws)
    z = (0, 0, 0)
    assert ws[z] % 2 == 1
    ws[z] -= 1
    mus = []
    for w in sorted(ws):
        while ws[w] > 0:
            neg = tuple(-x for x in w)
            if w == neg:
                assert ws[w] >= 2; ws[w] -= 2; mus.append(w)
            else:
                assert ws[neg] >= 1; ws[w] -= 1; ws[neg] -= 1; mus.append(w)
    assert len(mus) == 6
    C = Counter()
    for eps in itertools.product((1, -1), repeat=6):
        s = [F(0)]*3
        for e, mu in zip(eps, mus):
            for k in range(3): s[k] += F(e*mu[k], 2)
        assert all(x.denominator == 1 for x in s)
        C[tuple(int(x) for x in s)] += 1
    return C

def Xchar(U, anchor):
    """character in Z^4 of e8 restricted to sl2^4; anchor = index (0..3) of the 2A1 factor."""
    others = [k for k in range(4) if k != anchor]
    def emb(a, w3):
        v = [0]*4; v[anchor] = a
        for k, x in zip(others, w3): v[k] = x
        return tuple(v)
    C = Counter()
    for a in (-2, 0, 2): C[emb(a, (0, 0, 0))] += 1                      # (V2,1)
    wU = []
    for t, m in U.items(): wU += weights(t)*m
    for a in (-2, 0, 2):
        for w in wU: C[emb(a, w)] += 1                                     # (V2, U')
    for i in range(len(wU)):
        for j in range(i+1, len(wU)):
            C[emb(0, tuple(wU[i][k]+wU[j][k] for k in range(3)))] += 1     # Lambda^2 U'
    for w, m in spin_weights(U).items():
        for a in (-1, 1): C[emb(a, w)] += m                                # (V1, Spin)
    assert sum(C.values()) == 248
    return C

def decomp(C, r=4):
    C = Counter(C); out = Counter()
    while C:
        w = max(v for v in C if C[v] > 0)
        m = C[w]; out[w] += m
        for u in weights(w):
            C[u] -= m
            if C[u] == 0: del C[u]
            assert C.get(u, 0) >= 0
    return out

sols = []
# e8^(1): anchor g'_1 (factor 0); U'_1 rep of (g''_1, g'_2, g''_2) with indices (s, 1, t)
# e8^(2): anchor g'_2 (factor 2); U'_2 rep of (g'_1, g''_1, g''_2) with indices (1, 13-s, 13-t)
cache2 = {}
for s in range(14):
    for t in range(14):
        L1 = orth_reps((s, 1, t))
        if not L1: continue
        L2 = orth_reps((1, 13-s, 13-t))
        if not L2: continue
        X2 = {}
        for U2 in L2:
            key = frozenset(Xchar(U2, 2).items())
            X2.setdefault(key, []).append(U2)
        for U1 in L1:
            c1 = Xchar(U1, 0)
            need = Counter(charW1); need.subtract(c1)
            if min(need.values()) < 0: continue
            need = +need
            key = frozenset(need.items())
            if key in X2:
                for U2 in X2[key]:
                    sols.append((s, t, U1, U2))
print("number of solutions (s,t,U'_1,U'_2):", len(sols))
seen = set()
for s, t, U1, U2 in sols:
    print(f"  split e8^(1) indices (2,{s},1,{t}) ; e8^(2) indices (1,{13-s},2,{13-t})")
    print("     U'_1 [(g''_1,g'_2,g''_2)] =", sorted(U1.items()))
    print("     U'_2 [(g'_1,g''_1,g''_2)] =", sorted(U2.items()))
    print("     X_1 =", sorted(decomp(Xchar(U1, 0)).items()))
    print("     X_2 =", sorted(decomp(Xchar(U2, 2)).items()))
json.dump([{"s": s, "t": t, "U1": [[list(k), m] for k, m in sorted(U1.items())],
            "U2": [[list(k), m] for k, m in sorted(U2.items())]} for s, t, U1, U2 in sols],
          open("t12_e8split.json", "w"), indent=1)
