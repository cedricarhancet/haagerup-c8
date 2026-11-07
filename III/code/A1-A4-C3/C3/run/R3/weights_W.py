"""sl2^4-weights of the basis states of W_0, W_1, W_2 and the character of W_1 (R3a step 1)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, itertools
from collections import Counter
from gmpy2 import mpq
import e8e8_voa as V
HERE = os.path.dirname(os.path.abspath(__file__))
SL = pickle.load(open(os.path.join(HERE, 'sl2x4.pkl'), 'rb'))
NAMES = ['g1p', 'g1pp', 'g2p', 'g2pp']
# H vectors as elements of h (coefficients on eps_i)
HV = []
for nm in NAMES:
    H = SL[nm][1]; h = [mpq(0)]*16
    for (b, m), c in H.items():
        assert b == V.Z_IDX and len(m) == 1 and m[0][0] == 1
        h[m[0][1]] = c
    HV.append(h)
def wt_lat(ib):
    y = V.VECS[ib]
    w = []
    for h in HV:
        s = sum(h[i]*y[i] for i in range(16))/2
        assert s.denominator == 1
        w.append(int(s))
    return tuple(w)
LW = [wt_lat(i) for i in range(len(V.VECS))]
def weight_of(key): return LW[key[0]]
def weights(t): return list(itertools.product(*[range(-a, a+1, 2) for a in t]))
def decomp(Cc):
    Cc = Counter(Cc); out = Counter()
    while Cc:
        w = max(v for v in Cc if Cc[v] > 0)
        mlt = Cc[w]; out[w] += mlt
        for u in weights(w):
            Cc[u] -= mlt
            if Cc[u] == 0: del Cc[u]
            assert Cc.get(u, 0) >= 0
    return out
if __name__ == '__main__':
    B1 = V.basis(1)
    ch = Counter(weight_of(k) for k in B1)
    dec = decomp(ch)
    t12 = json.load(open(os.path.join(HERE, '..', 'data', 't12_W1.json')))['W1']
    t12 = {tuple(map(int, k.split(','))): m for k, m in t12.items()}
    print('W_1 under sl2^4 (g1p, g1pp, g2p, g2pp):', sorted(dec.items()))
    print('equal to T12 W1 data:', dict(dec) == t12)
    B2 = V.basis(2)
    ch2 = Counter(weight_of(k) for k in B2)
    print('dim of weight (0,0,0,0) in W_2:', ch2[(0, 0, 0, 0)])
    for mu in [(0, 4, 0, 0), (1, 1, 0, 0), (1, 3, 0, 0), (1, 5, 0, 0), (1, 7, 0, 0), (2, 2, 0, 0), (2, 6, 0, 0), (3, 3, 0, 0), (2, 0, 0, 0), (0, 2, 0, 0), (4, 0, 0, 0)]:
        print('  dim weight', mu, ':', ch2[mu])
