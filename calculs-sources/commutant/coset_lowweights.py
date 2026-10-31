"""Low-weight dimensions (weights <= 4) of coset vertex algebras, computed in the universal affine algebras (exact
below the first null vector: level k vacuum null vector at weight k+1).
 (a) parafermion K(sl2, 17) = Com(Heisenberg, V^17(sl2));
 (b) GKO coset Com(V^17_diag, V^4 (x) V^13);
 (c) product (a) (x) (b), to be compared with C = Com_S(A): dim C_n = 1, 0, 2, 2 (exact: lower bound from the two
     Virasoro algebras, upper bound modulo p), dim C_4 <= 7 (mod p) and >= 5."""
from collections import defaultdict
N = 4
def mul(A, B):   # sl2 characters as dict highest weight -> multiplicity
    out = defaultdict(int)
    for a, x in A.items():
        for b, y in B.items():
            for c in range(abs(a - b), a + b + 1, 2): out[c] += x * y
    return dict(out)
def add(A, B, s=1):
    out = defaultdict(int, A)
    for b, y in B.items(): out[b] += s * y
    return {k: v for k, v in out.items() if v}
def sym_power_adj(n):   # character of S^n(V(2)) via weights
    from itertools import combinations_with_replacement
    wts = defaultdict(int)
    for c in combinations_with_replacement([2, 0, -2], n): wts[sum(c)] += 1
    return from_weights(wts)
def from_weights(wts):
    wts = dict(wts); out = {}
    while any(wts.values()):
        m = max(w for w, v in wts.items() if v); mult = wts[m]; out[m] = mult
        for w in range(-m, m + 1, 2): wts[w] = wts.get(w, 0) - mult
    return out
def weights(A):
    w = defaultdict(int)
    for a, x in A.items():
        for j in range(-a, a + 1, 2): w[j] += x
    return w
# universal vacuum module U(g[t^-1]t^-1) as graded g-module: PBW over modes -1..-N of adj
def vacuum_chars(N):
    # weight-graded symmetric algebra: product over n of S(adj in degree n)
    res = {0: {0: 1}}
    for n in range(1, N + 1):
        new = defaultdict(dict)
        for w, A in res.items():
            for m in range(0, (N - w) // n + 1):
                S = sym_power_adj(m) if m else {0: 1}
                tgt = w + n * m
                new[tgt] = add(new.get(tgt, {}), mul(A, S))
        res = dict(new)
    return [res.get(w, {}) for w in range(N + 1)]
VAC = vacuum_chars(N)
# (b) GKO: b_lambda(n) from char(V(x)V)_n = sum_lambda sum_m VAC_{n-m} (x) V(lambda) * b_lambda(m)
TT = [{} for _ in range(N + 1)]
for i in range(N + 1):
    for j in range(N + 1 - i): TT[i + j] = add(TT[i + j], mul(VAC[i], VAC[j]))
b = [dict() for _ in range(N + 1)]
for n in range(N + 1):
    rest = dict(TT[n])
    for m in range(n):
        for lam, x in b[m].items(): rest = add(rest, mul(VAC[n - m], {lam: x}), -1)
    b[n] = rest
gko = [b[n].get(0, 0) for n in range(N + 1)]
# (a) parafermion: charge-0 weight spaces of vacuum module, divided by the Heisenberg Fock character
ch0 = [weights(VAC[n]).get(0, 0) for n in range(N + 1)]
fock = [1, 1, 2, 3, 5][:N + 1]
pf = []
for n in range(N + 1): pf.append(ch0[n] - sum(fock[n - m] * pf[m] for m in range(n)))
prod_ = [sum(pf[m] * gko[n - m] for m in range(n + 1)) for n in range(N + 1)]
vir2 = [1, 0, 2, 2, 5]
print('parafermion K(sl2,17)    :', pf)
print('GKO (sl2_4 x sl2_13)/sl2_17:', gko)
print('product                  :', prod_)
print('Vir x Vir (generic c)    :', vir2)
print('C = Com_S(A)             : [1, 0, 2, 2, <=7]')
