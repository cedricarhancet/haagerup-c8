"""Exhaustiveness of the T12 family (R3a): each level-3 factor g'_k has Dynkin index split (2,1) or (1,2) between the
two copies of e8 (index 3 in one copy forces 3A1, which contains V(3), excluded; index 2 forces 2A1, index 1 forces A1).
T12 treats g'_1 = (2,1), g'_2 = (1,2) (and the swap).  Remaining option B: g'_1 = g'_2 = (2,1).  We show that the
restriction of W_1 to g'_1 + g'_2 is then incompatible with the T12 data, by comparing characters in Z^2.
Copy with both index 2: g'_1 = 2A1 with centraliser so(13), g'_2 in so(13) with 13-dim orthogonal representation U' of
index 4 (2 x Dynkin index in e8): U' = V2 + 10 V0 or 4 V1 + 5 V0;  e8 | = (V2,1) + V2 x U' + Lambda^2 U' + V1 x Spin(U').
Copy with both index 1: two commuting A1 of index 1 = two orthogonal roots (A1 orbit), e8 | = (2,0)+(0,2)+12(1,1)+32(1,0)+32(0,1)+66(0,0)."""
import json, itertools, os
from collections import Counter
from fractions import Fraction as F
W1 = {tuple(map(int, k.split(','))): m for k, m in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 't12_W1.json')))['W1'].items()}
def weights(t): return list(itertools.product(*[range(-a, a+1, 2) for a in t]))
target = Counter()
for (a, b, c, d), m in W1.items():
    for wa in range(-a, a+1, 2):
        for wc in range(-c, c+1, 2):
            target[(wa, wc)] += m*(b+1)*(d+1)
def spin_weights_1(Uw):
    ws = Counter(Uw); ws[0] -= 1
    mus = []
    for w in sorted(ws):
        while ws[w] > 0:
            if w == 0: ws[w] -= 2; mus.append(0)
            else: ws[w] -= 1; ws[-w] -= 1; mus.append(w)
    C = Counter()
    for eps in itertools.product((1, -1), repeat=6):
        s = sum(F(e*mu, 2) for e, mu in zip(eps, mus)); assert s.denominator == 1; C[int(s)] += 1
    return C
results = []
for name, Uw in (('V2+10V0', [2, 0, -2] + [0]*10), ('4V1+5V0', [1, -1]*4 + [0]*5)):
    C1 = Counter()
    for a in (-2, 0, 2): C1[(a, 0)] += 1
    for a in (-2, 0, 2):
        for w in Uw: C1[(a, w)] += 1
    for i in range(13):
        for j in range(i+1, 13): C1[(0, Uw[i] + Uw[j])] += 1
    for w, m in spin_weights_1(Uw).items():
        for a in (-1, 1): C1[(a, w)] += m
    assert sum(C1.values()) == 248
    C2 = Counter()
    for t, m in {(2, 0): 1, (0, 2): 1, (1, 1): 12, (1, 0): 32, (0, 1): 32, (0, 0): 66}.items():
        for w in weights(t): C2[w] += m
    tot = C1 + C2
    ok = (tot == target)
    print('option B, U\'|g\'_2 =', name, ': character of W_1 on (g\'_1, g\'_2) equals the T12 data:', ok)
    results.append(ok)
print('T12 data restricted to (g\'_1, g\'_2):', sorted(target.items()))
print('option B excluded:', not any(results))
