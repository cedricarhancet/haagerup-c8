"""Write the exact relation lists of weights 2, 4, 6 (gamma = 0, 1) and the exact ranks of the coefficient maps."""
import json
import obstruction as O, vvmf_exact as V
names = V.names; d = V.d
labels = ['w0', 'w0p'] + ['w_' + names[z] for z in range(1, d)]
out = {}
for gamma in (0, 1):
    P = [0, names.index('d3')] if gamma == 1 else [0]
    A = [(0, -1), (0, 0)] + [(z, -1 if z in P else 0) for z in range(1, d)]
    for k in (2, 4, 6):
        R = O.relations(k, P, A, 1)
        F = O.exact_forms(k, P, A, 1)
        rk = V.rank(F, len(A)) if F else 0
        print('gamma %d weight %d: %d relations, exact forms give rank %d (sum %d)' % (gamma, k, len(R), rk, len(R) + rk))
        json.dump({'gamma': gamma, 'weight': k, 'labels': labels, 'relations': [[str(c) for c in r] for r in R],
                   'rank_of_exact_forms': rk}, open('weight%d_relations_exact_gamma%d.json' % (k, gamma), 'w'), indent=1)
