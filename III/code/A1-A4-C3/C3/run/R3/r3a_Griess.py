"""R3a step 2b: Griess algebra of G = Com_W(A x A)_2 (exact): structure constants u_(1) v, forms u_(3) v, and the
matrices of the zero modes o(g) on the highest weight spaces of all sl2^4-components of W_1.
Output: G_griess.pkl"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, time
from gmpy2 import mpq
import e8e8_voa as V
import common as CM
t0 = time.time()
D = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
G = D['G']; n = len(G)
prod = {}
for i in range(n):
    for j in range(i, n):
        p = V.clean(V.mode_vec(G[i], 1, G[j]))
        c = CM.span_coords(G, p)
        assert c is not None, 'G not closed under (1)'
        prod[(i, j)] = prod[(j, i)] = [CM.to_mq(x) for x in c]
    print('row', i, '(%.0fs)' % (time.time() - t0), flush=True)
form = [[V.mode_vec(G[i], 3, G[j]).get(V.VAC, mpq(0)) for j in range(n)] for i in range(n)]
comp = CM.w1_sectors()
Mk = {}
for kappa in sorted(comp):
    Mk[kappa] = [[[CM.to_mq(x) for x in row] for row in CM.zero_mode_matrix(g, kappa)] for g in G]
omCC = CM.span_coords(G, D['omCC'])
pickle.dump({'prod': prod, 'form': form, 'Mk': Mk, 'comp': dict(comp), 'omCC': [CM.to_mq(x) for x in omCC]},
            open(os.path.join(CM.HERE, 'G_griess.pkl'), 'wb'))
print('done (%.0fs)' % (time.time() - t0))
