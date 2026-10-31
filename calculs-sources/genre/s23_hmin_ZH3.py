"""Axis G (v12), side result.  Let V be strongly rational with c(V) = 8 and Rep(V) braided equivalent to Z(H_3).
Since Z(H_3) has trivial universal grading group, its spherical structure is the pseudo-unitary one, and since
c(V) = 0 = c(Z(H_3)) mod 8 the character S-matrix is +s ([MoR24a, Lemma 3.6, Prop. 3.5]).  If h_min < 0, the set J of
simple objects of minimal conformal weight is a nonempty set of objects with equal twists, not containing 0 (all objects
are self-dual).  For every such J we look for a row i with s_ij <= 0 for all j in J and s_ij < 0 for some j in J; then
ch_i(iy) ~ (sum_j s_ij dim V^j_{h_min}) exp(2 pi (c/24 - h_min)/y) < 0 as y -> 0+, a contradiction.  Output: witnesses."""
import itertools, json, os
from mpmath import mp, mpf, sqrt, cos, pi
mp.dps = 50
HERE = os.path.abspath(os.path.dirname(__file__))
lab = ['0', 'b', 'a', 'c0', 'c1', 'c2'] + ['d%d' % l for l in range(1, 7)]
x = (13 - 3*sqrt(13))/26; y = 3/sqrt(13)
def cc(j): return -2*y*cos(2*pi*j/13)
M = [[mpf(0)]*12 for _ in range(12)]
M[0][:6] = [x, 1 - x, 1, 1, 1, 1]; M[1][:6] = [1 - x, x, 1, 1, 1, 1]
for l in range(6): M[0][6 + l] = y; M[1][6 + l] = -y
blk = [[2, -1, -1, -1], [-1, 2, -1, -1], [-1, -1, -1, 2], [-1, -1, 2, -1]]
for i in range(4):
    M[2 + i][0] = 1; M[2 + i][1] = 1
    for j in range(4): M[2 + i][2 + j] = blk[i][j]
for l in range(6):
    M[6 + l][0] = y; M[6 + l][1] = -y
    for m in range(6): M[6 + l][6 + m] = cc((l + 1)*(m + 1))
for i in range(12):
    for j in range(12):
        if i < 6 and j >= 6 and i >= 2: M[i][j] = 0
s = [[mpf(M[i][j])/3 for j in range(12)] for i in range(12)]
sym = max(abs(s[i][j] - s[j][i]) for i in range(12) for j in range(12))
uni = max(abs(sum(s[i][k]*s[k][j] for k in range(12)) - (1 if i == j else 0)) for i in range(12) for j in range(12))
classes = [[1, 2, 3], [4], [5]] + [[6 + l] for l in range(6)]     # twist classes without the unit (0 has twist 1 with b, a, c0)
res = []
ok = True
for cl in classes:
    for r in range(1, len(cl) + 1):
        for J in itertools.combinations(cl, r):
            wit = [lab[i] for i in range(12) if all(s[i][j] <= mpf(10)**-40 for j in J) and any(s[i][j] < -mpf(10)**-40 for j in J)]
            res.append({'J': [lab[j] for j in J], 'witness_rows': wit}); ok = ok and bool(wit)
for r in res: print(r)
print('S symmetric:', float(sym), ' S^2 = 1:', float(uni), ' every J has a witness:', ok)
json.dump({'checks': {'symmetry_defect': float(sym), 'S2_defect': float(uni)}, 'J_witnesses': res, 'all_J_excluded': ok},
          open(HERE + '/s23_hmin_ZH3.json', 'w'), indent=1)
