"""Cache the Griess algebra B (case II) in the basis [e] + e^perp used by griess_subalgebras.json."""
import json, time, pickle
import griess_subalgebras as GS
t0 = time.time()
T, basis, prod, gram = GS.griess_data('II')
json.dump({'prod': {f'{i},{j}': [str(c) for c in v] for (i, j), v in prod.items()}, 'gram': [[str(c) for c in r] for r in gram]},
          open('griess_B_II.json', 'w'), indent=0)
pickle.dump(basis, open('griess_B_II_basis.pkl', 'wb'))
print('done %.0fs' % (time.time() - t0))
