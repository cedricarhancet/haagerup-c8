"""Decomposition of the explicitly known lowest spaces R_b, R_a = R_c0, R_c1, R_c2 under s = g' + g'' in the two
cases of the classification (joint weights from e8_joint.py)."""
import json
from collections import Counter
import e8_joint as J
from u2_embedding_test import decompose2
sl3ad = [2, 2, -2, -2, 4, -4, 0, 0]; w3 = [2, 0, -2]
out = {}
for name, lab, i1p in (('A11+A139', (2, 1, 1, 0, 1, 2), 1), ('A13+A113', (1, 0, 0, 1, 0, 1), 3)):
    ad6, w27 = J.joint_weights(lab, i1p)
    Rb = Counter(decompose2(ad6 + [(0, s) for s in sl3ad]))
    Rb[(2, 0)] -= 1; Rb[(0, 2)] -= 1          # remove s = V(2,0) + V(0,2)
    Rb = +Rb
    Ra = decompose2([(x[0], x[1] + y) for x in w27 for y in w3])
    Rc1 = decompose2([(0, y) for y in w3]); Rc2 = decompose2(w27)
    out[name] = {z: {str(k): v for k, v in sorted(d.items())} for z, d in (('b', Rb), ('a', Ra), ('c1', Rc1), ('c2', Rc2))}
    print(name)
    for z, d in (('b', Rb), ('a=c0', Ra), ('c1', Rc1), ('c2', Rc2)):
        print('   R_%s =' % z, dict(sorted(d.items())), ' dim', sum((a + 1)*(b + 1)*m for (a, b), m in d.items()))
json.dump(out, open('known_sectors_table.json', 'w'), indent=1)
