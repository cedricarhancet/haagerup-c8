"""Axis I2: Weyl group element of E6 (+ diagram automorphism) mapping Jiao's Cartan pair (P_H, B_H), transported by the
diagram isometry phi (Jiao 1,2,3,4,5,6,7,8 -> ours 1,3,4,5,6,2,7,8), onto our pair (H1, H2). Vectors in root coordinates."""
import sys, json
sys.path.insert(0, '../generation-forte'); import os; os.chdir('../generation-forte')
import lattice_voa as V
G = V.G
PERM = [0, 2, 3, 4, 5, 1, 6, 7]                       # Jiao index i (0-based) -> our index PERM[i]
def phi(v): out = [0]*8; [out.__setitem__(PERM[i], v[i]) for i in range(8)]; return tuple(out)
PH = (1, 2, 3, 2, 1, 0, 0, 0); BH = (3, 4, 3, 4, 3, 2, 2, 2)
H1 = (0, -1, 1, 0, 1, 0, 0, 0); H2 = (4, 5, 7, 10, 7, 4, 2, 2)
def ip(a, b): return sum(a[i]*G[i][j]*b[j] for i in range(8) for j in range(8))
# Jiao Gram (edges 12,23,34,45,36,78) must be the pull-back of ours
GJ = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
for a, b in [(0,1),(1,2),(2,3),(3,4),(2,5),(6,7)]: GJ[a][b] = GJ[b][a] = -1
assert all(G[PERM[i]][PERM[j]] == GJ[i][j] for i in range(8) for j in range(8)), 'phi is not an isometry'
P, B = phi(PH), phi(BH)
print('norms', ip(P, P), ip(H1, H1), ip(B, B), ip(H2, H2), ' cross', ip(P, B), ip(H1, H2))
def refl(v, i):
    c = sum(v[j]*G[j][i] for j in range(8))
    w = list(v); w[i] -= c; return tuple(w)
DIAG = [5, 1, 4, 3, 2, 0, 6, 7]                       # E6 diagram automorphism in our (Bourbaki) numbering: 1<->6, 3<->5
def diag(v): out = [0]*8; [out.__setitem__(DIAG[i], v[i]) for i in range(8)]; return tuple(out)
res = None
for d in (0, 1):
    start = (diag(P), diag(B)) if d else (P, B)
    seen = {start: []}; frontier = [start]
    while frontier and res is None:
        nxt = []
        for pair in frontier:
            for i in range(6):
                q = (refl(pair[0], i), refl(pair[1], i))
                if q not in seen:
                    seen[q] = seen[pair] + [i]; nxt.append(q)
                    if q == (H1, H2): res = (d, seen[q]); break
            if res: break
        frontier = nxt
    print('diagram', d, 'orbit explored', len(seen), 'found' if res else '')
    if res: break
print('result', res)
json.dump({'diagram': res[0] if res else None, 'word': res[1] if res else None}, open('../unicite/i2_weyl.json', 'w'))
