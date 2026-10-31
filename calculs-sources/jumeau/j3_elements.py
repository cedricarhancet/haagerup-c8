"""Axis B1: the 8 elements of Z (case II) as matrices on V'_1 in the lattice basis; are they toral / monomial?"""
import sys, os, pickle, itertools, json
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath(HERE + '/../generation-forte'))
os.chdir(HERE + '/../generation-forte')
import sympy as sp
import lattice_voa as V
CASE = 'II'
d = pickle.load(open(HERE + '/j1_equations_%s.pkl' % CASE, 'rb'))
label = d['label']; P = sp.Matrix.hstack(*[sp.Matrix(c) for c in d['P']]); n = P.shape[0]
types = pickle.load(open(HERE + '/j1_types_%s.pkl' % CASE, 'rb'))['types']
B1 = V.basis(1)
def zmat(e0, e1, e2):
    z22 = e1*e2
    Zb = {(0, 2): sp.eye(3), (2, 0): sp.eye(1), (0, 4): sp.diag(e0, z22), (1, 1): sp.Matrix([[e1]]),
          (1, 3): sp.Matrix([[e2, 0], [sp.Rational(e2 - e1, 8), e1]]), (1, 5): sp.Matrix([[e2]]),
          (2, 2): sp.Matrix([[z22]]), (2, 4): sp.Matrix([[1]]), (3, 1): sp.Matrix([[e2]])}
    pos = {tuple(l) if not isinstance(l[0], tuple) else (tuple(l[0]), l[1], l[2], l[3]): i for i, l in enumerate(label)}
    Zf = sp.zeros(n, n)
    for (mu, k, i, j), p in pos.items():
        for l in range(Zb[mu].shape[0]):
            Zf[pos[(mu, l, i, j)], p] = Zb[mu][l, k]
    return P*Zf*P.inv()
out = {}
for e0, e1, e2 in itertools.product((1, -1), repeat=3):
    M = zmat(e0, e1, e2)
    diag = all(M[i, j] == 0 for i in range(n) for j in range(n) if i != j)
    # action on root vectors: B1 keys (beta, ()) for roots, Cartan keys (ZERO, ((1,i),))
    rootsigns = {}
    mono = True
    for j, k in enumerate(B1):
        col = [i for i in range(n) if M[i, j] != 0]
        if len(col) != 1: mono = False
    cart = M[[i for i, k in enumerate(B1) if k[0] == V.ZERO], [i for i, k in enumerate(B1) if k[0] == V.ZERO]]
    key = (e0, e1, e2)
    out[str(key)] = {'diagonal': diag, 'monomial': mono, 'cartan_block_identity': cart == sp.eye(8)}
    if diag:
        sg = {}
        for j, k in enumerate(B1):
            if k[0] != V.ZERO: sg[str(k[0])] = int(M[j, j])
        out[str(key)]['root_signs'] = sg
    print(key, out[str(key)]['diagonal'], out[str(key)]['monomial'], out[str(key)]['cartan_block_identity'])
    pickle.dump({str(k): v for k, v in out.items()}, open(HERE + '/j3_elements.pkl', 'wb'))
    sp_mat = [[str(x) for x in M.row(i)] for i in range(n)]
    pickle.dump(sp_mat, open(HERE + '/j3_z_%d%d%d.pkl' % ((e0 + 1)//2, (e1 + 1)//2, (e2 + 1)//2), 'wb'))
