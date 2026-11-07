"""Independent check of Lemma sl2-four: (a) option B exclusion; (b) uniqueness of (U1,U2).
Factors order: 0 = g'x1, 1 = g''x1, 2 = 1xg', 3 = 1xg''."""
import itertools
from collections import Counter
def wts(n): return list(range(-n, n+1, 2))
def ind(n): return n*(n+1)*(n+2)//6
def char_irr(hw):  # Counter over Z^k
    return Counter(itertools.product(*[wts(n) for n in hw]))
def add(c, d, m=1):
    for k, v in d.items(): c[k] += m*v
def tensor(c1, c2):
    out = Counter()
    for a, x in c1.items():
        for b, y in c2.items():
            out[tuple(a) + tuple(b)] += x*y
    return out
# ---------------- target character of W_1 (Table 2 case II, pi from Prop diagonal)
R = {'c1': {(0,2):1}, 'c2': {(0,0):1,(0,4):1,(1,1):1,(1,3):1,(2,2):1}, 'd1': {(0,2):1,(1,1):1},
     'd5': {(0,4):1,(1,3):1}, 'd4': {(0,0):1,(1,1):1}, 'd6': {(0,2):1,(1,3):1,(2,0):1}}
PI = {'c1':'c2','c2':'c1','d1':'d5','d5':'d1','d4':'d6','d6':'d4'}
target = Counter()
for hw in [(2,0,0,0),(0,2,0,0),(0,0,2,0),(0,0,0,2)]: add(target, char_irr(hw))
for z, pz in PI.items():
    for (a,b), m in R[z].items():
        for (c,d), n in R[pz].items():
            add(target, char_irr((a,b,c,d)), m*n)
print('dim target =', sum(target.values()))
def decompose(ch):
    ch = Counter(ch); out = Counter()
    while True:
        ch = Counter({k:v for k,v in ch.items() if v})
        if not ch: return out
        top = max(ch, key=lambda w: (sum(w), w))
        # highest weight: maximal in dominance; pick lexicographic max among dominant maximal
        cands = [w for w in ch if all(x >= 0 for x in w)]
        top = max(cands, key=lambda w: (sum(w), w))
        m = ch[top]; assert m > 0
        out[top] += m
        add(ch, char_irr(top), -m)
        assert all(v >= 0 for v in ch.values()), 'negative'
# ---------------- e8 under a + so(13) acting via U (rep of sl2^3 on C^13)
def spin_char(Uchar, k):
    """Uchar: Counter of weights (k-tuples) of an orthogonal rep of odd dim; spin rep character (weights doubled->halved)"""
    ws = []
    for w, m in Uchar.items(): ws += [w]*m
    pos = []; zero = 0
    rem = Counter(ws)
    for w in sorted(rem):
        pass
    # pair each nonzero weight with its negative
    cnt = Counter(ws)
    for w in list(cnt):
        if all(x == 0 for x in w): continue
        neg = tuple(-x for x in w)
        if w > neg:
            assert cnt[w] == cnt[neg]
            pos += [w]*cnt[w]
    nz0 = cnt[tuple([0]*k)]
    # zero weights: (nz0-1)/2 pairs of zero mu's
    assert nz0 % 2 == 1
    pos += [tuple([0]*k)]*((nz0-1)//2)
    assert len(pos) == 6, len(pos)
    ch = Counter({tuple([0]*k): 1})
    for mu in pos:
        # e^{mu/2} + e^{-mu/2}; mu has integer coords; weights of sl2 normalized as integers already (h-eigenvalues),
        # so mu/2 must be integral on sum: we keep h-eigenvalues: spin weights = (1/2) sum +-mu_i
        new = Counter()
        for w, m in ch.items():
            new[tuple(a + b for a, b in zip(w, mu))] += m
            new[tuple(a - b for a, b in zip(w, mu))] += m
        ch = new
    out = Counter()
    for w, m in ch.items():
        assert all(x % 2 == 0 for x in w)
        out[tuple(x//2 for x in w)] += m
    return out
def lam2(Uchar):
    ws = []
    for w, m in Uchar.items(): ws += [w]*m
    out = Counter()
    for i, j in itertools.combinations(range(len(ws)), 2):
        out[tuple(a+b for a, b in zip(ws[i], ws[j]))] += 1
    return out
def e8_char(Uchar, k):
    """character of e8 under sl2 (=a, first coordinate) x sl2^k (via U on E_3^perp)"""
    z = tuple([0]*k)
    out = Counter()
    V2 = char_irr((2,)); V1 = char_irr((1,))
    add(out, tensor(V2, Counter({z:1})))               # Lambda^2 E3 = adjoint of a
    add(out, tensor(V2, Uchar))                        # E3 x U
    add(out, tensor(Counter({(0,):1}), lam2(Uchar)))   # Lambda^2 U (concatenation: a-weight 0)
    add(out, tensor(V1, spin_char(Uchar, k)))          # Delta_3 x Delta_13
    assert sum(out.values()) == 248, sum(out.values())
    return out
# ---------------- (a) option B
def to4(ch, order):  # order: list mapping local coords -> global factor index; others weight 0
    out = Counter()
    for w, m in ch.items():
        g = [0]*4
        for loc, glob in enumerate(order): g[glob] = w[loc]
        out[tuple(g)] += m
    return out
# restricted target to factors 0 and 2
tgt02 = Counter()
for w, m in target.items(): tgt02[(w[0], w[2])] += m
for name, Uhw in (('V(2)+10V(0)', {(2,):1, (0,):10}), ('4V(1)+5V(0)', {(1,):4, (0,):5})):
    U = Counter()
    for hw, m in Uhw.items(): add(U, char_irr(hw), m)
    c1 = e8_char(U, 1)     # coords (g'x1 = a, 1xg' via U)
    # copy 2: two commuting index-1 sl2: (2,0)+(0,2)+12(1,1)+32(1,0)+32(0,1)+66(0,0)
    c2 = Counter()
    for hw, m in {(2,0):1,(0,2):1,(1,1):12,(1,0):32,(0,1):32,(0,0):66}.items(): add(c2, char_irr(hw), m)
    tot = Counter(); add(tot, c1); add(tot, c2)
    print('option B,', name, ': equal to target restricted to (g\'x1, 1xg\') ?', tot == tgt02)
    print('   decomposition option B:', dict(decompose(tot)))
print('   target (g\'x1,1xg\'):', dict(decompose(tgt02)))
# ---------------- (b) enumeration of (U1, U2)
irrs = [hw for hw in itertools.product(range(5), repeat=3) if (hw[0]+1)*(hw[1]+1)*(hw[2]+1) <= 13]
def dimh(hw): return (hw[0]+1)*(hw[1]+1)*(hw[2]+1)
def reps13():
    out = []
    def rec(i, rem, cur):
        if rem == 0: out.append(dict(cur)); return
        if i == len(irrs): return
        hw = irrs[i]; d = dimh(hw)
        for m in range(rem//d, -1, -1):
            if m: cur[hw] = m
            rec(i+1, rem - m*d, cur)
            if m: del cur[hw]
    rec(0, 13, {})
    return out
allreps = reps13()
def orthogonal(r): return all(m % 2 == 0 for hw, m in r.items() if sum(hw) % 2 == 1)
def indices(r):
    out = [0,0,0]
    for hw, m in r.items():
        for f in range(3):
            o = [hw[g]+1 for g in range(3) if g != f]
            out[f] += m*ind(hw[f])*o[0]*o[1]
    return out
orth = [r for r in allreps if orthogonal(r)]
print('13-dim reps of sl2^3:', len(allreps), '; orthogonal:', len(orth))
# U1: factors (g''x1, 1xg', 1xg'') ; rep-indices (2 j1, 2, 2 k1) ; U2: factors (g'x1, g''x1, 1xg'') : (2, 2(13-j1), 2(13-k1))
U1s = [r for r in orth if indices(r)[1] == 2 and indices(r)[0] % 2 == 0 and indices(r)[2] % 2 == 0]
U2s = [r for r in orth if indices(r)[0] == 2 and indices(r)[1] % 2 == 0 and indices(r)[2] % 2 == 0]
print('candidates U1:', len(U1s), ' U2:', len(U2s))
def uchar(r):
    c = Counter()
    for hw, m in r.items(): add(c, char_irr(hw), m)
    return c
cache1 = {}
sols = []
for r1 in U1s:
    i1 = indices(r1)
    for r2 in U2s:
        i2 = indices(r2)
        if i1[0] + i2[1] != 26 or i1[2] + i2[2] != 26: continue
        k1 = id(r1)
        if k1 not in cache1: cache1[k1] = to4(e8_char(uchar(r1), 3), [0, 1, 2, 3])
        c1 = cache1[k1]
        c2 = to4(e8_char(uchar(r2), 3), [2, 0, 1, 3])
        tot = Counter(); add(tot, c1); add(tot, c2)
        if tot == target: sols.append((r1, r2, i1, i2))
print('solutions:', len(sols))
for s in sols: print('  U1 =', s[0], 'idx', s[2], '| U2 =', s[1], 'idx', s[3])
