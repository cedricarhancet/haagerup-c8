"""Realization of U_1 = g' + g'' inside V'_1 for the two cases of Theorem 8.5 (v03), with H', H'' in the lattice
Cartan subalgebra, so that every basis state of V' is a joint weight vector.
Case I : H''_e6 labels (2,1,1,0,1,2) (orbit A5, index 35), levels (1, 39).
Case II: H''_e6 labels (1,0,0,1,0,1) (orbit 2A2+A1, index 9), levels (3, 13).
H'_t: labels from cartan_conjugate.json (v03).  p2(g'') = principal so3 of sl3."""
import random, json
from gmpy2 import mpq
import sympy as sp
from flint import fmpq_mat, fmpq
import lattice_voa as V
CASES = {'I': {'lab2': (2, 1, 1, 0, 1, 2), 'levels': (1, 39)},
         'II': {'lab2': (1, 0, 0, 1, 0, 1), 'levels': (3, 13)}}
JS = json.load(open('cartan_conjugate.json'))
CASES['I']['lab1'] = tuple(JS['A11+A139'][0]); CASES['II']['lab1'] = tuple(JS['A13+A113'][0])
CE6inv = sp.Matrix(V.CE6).inv()
def cartan_from_labels(lab6, labA2=(0, 0)):
    c = list(sp.Matrix([lab6])*CE6inv) + list(sp.Matrix([labA2])*sp.Matrix([[2, -1], [-1, 2]]).inv())
    return [mpq(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in c]
def hstate(c):
    return {(V.ZERO, ((1, i),)): x for i, x in enumerate(c) if x}
def hval(c, beta):          # (H, beta) for H = sum c_i alpha_i
    return sum(c[i]*V.G[i][j]*beta[j] for i in range(8) for j in range(8) if c[i] and beta[j])
def br(x, y): return V.clean(V.mode_vec(x, 0, y))
def solve(cols, target, keys=None):
    """solve sum t_j cols[j] = target exactly; returns list or None"""
    if keys is None:
        keys = sorted(set().union(*[set(c) for c in cols]) | set(target))
    M = fmpq_mat(len(keys), len(cols), [fmpq(int(cols[j].get(k, 0).numerator), int(cols[j].get(k, 0).denominator)) for k in keys for j in range(len(cols))])
    b = fmpq_mat(len(keys), 1, [fmpq(int(target.get(k, 0).numerator), int(target.get(k, 0).denominator)) for k in keys])
    try:
        x = M.solve(b, algorithm=None)
    except Exception:
        # least: use rref on augmented
        A = fmpq_mat(len(keys), len(cols) + 1, [M[i, j] if j < len(cols) else b[i, 0] for i in range(len(keys)) for j in range(len(cols) + 1)])
        Rr, rank = A.rref()
        sol = [fmpq(0)]*len(cols); piv = []
        for i in range(rank):
            j = next(j for j in range(len(cols) + 1) if Rr[i, j] != 0)
            if j == len(cols): return None
            sol[j] = Rr[i, len(cols)]
        x = fmpq_mat(len(cols), 1, sol)
    sol = [mpq(int(x[j, 0].p), int(x[j, 0].q)) for j in range(len(cols))]
    chk = {}
    for j, c in enumerate(cols): V.vadd(chk, c, sol[j])
    if V.clean(chk) != V.clean(target): return None
    return sol
def lin(cols, coef):
    out = {}
    for c, t in zip(cols, coef):
        if t: V.vadd(out, c, t)
    return V.clean(out)

def build(case, seed0=0):
    D = CASES[case]
    c2 = cartan_from_labels(D['lab2']); cs = cartan_from_labels((0,)*6, (2, 2)); c1 = cartan_from_labels(D['lab1'])
    H2e6 = hstate(c2); Hsl3 = hstate(cs); H1 = hstate(c1)
    E6R = [a for a in V.ROOTS if a[6] == 0 and a[7] == 0]
    A2R = [a for a in V.ROOTS if all(x == 0 for x in a[:6])]
    pos2 = [a for a in E6R if hval(c2, a) == 2 and hval(c1, a) == 0]
    neg2 = [V.neg(a) for a in pos2]
    for seed in range(seed0, seed0 + 200):
        random.seed(seed)
        e = {(a, ()): mpq(random.randint(1, 5)) for a in pos2}
        cols = [br(e, {(b, ()): mpq(1)}) for b in neg2]
        sol = solve(cols, H2e6)
        if sol is None: continue
        f = V.clean({(b, ()): s for b, s in zip(neg2, sol)})
        if br(H2e6, e) == {k: 2*v for k, v in e.items()} and br(e, f) == V.clean(H2e6):
            break
    else:
        raise RuntimeError('no e6 triple')
    # sl3 principal triple
    es = {(a, ()): mpq(1) for a in A2R if hval(cs, a) == 2}
    negs = [V.neg(a) for a in A2R if hval(cs, a) == 2]
    sol = solve([br(es, {(b, ()): mpq(1)}) for b in negs], Hsl3)
    fs = V.clean({(b, ()): s for b, s in zip(negs, sol)})
    e2 = V.clean(V.vadd(dict(e), es)); f2 = V.clean(V.vadd(dict(f), fs)); H2 = V.clean(V.vadd(dict(H2e6), Hsl3))
    # g': root vectors of e6 with H'-eigenvalue 2, H''-eigenvalue 0, commuting with e, f
    p1 = [a for a in E6R if hval(c1, a) == 2 and hval(c2, a) == 0]
    # unknown e' = sum t_a e^a ; conditions [e', e] = 0 and [e', f] = 0
    imgs = [V.vadd(dict(br({(a, ()): mpq(1)}, e)), {}) for a in p1]
    imgs2 = [br({(a, ()): mpq(1)}, f) for a in p1]
    keys = sorted(set().union(*[set(x) for x in imgs]) )
    keys2 = sorted(set().union(*[set(x) for x in imgs2]))
    rows = [[imgs[j].get(k, 0) for j in range(len(p1))] for k in keys] + [[imgs2[j].get(k, 0) for j in range(len(p1))] for k in keys2]
    M = sp.Matrix(len(rows), len(p1), [sp.Rational(int(x.numerator), int(x.denominator)) for r in rows for x in r]) if rows else sp.zeros(1, len(p1))
    ns = M.nullspace()
    assert len(ns) == 1, ('dim of candidate e-prime space', len(ns))
    e1 = V.clean({(a, ()): mpq(int(sp.fraction(ns[0][j])[0]), int(sp.fraction(ns[0][j])[1])) for j, a in enumerate(p1) if ns[0][j] != 0})
    n1 = [V.neg(a) for a in p1]
    sol = solve([br(e1, {(b, ()): mpq(1)}) for b in n1], H1)
    assert sol is not None, 'H-prime not in [e-prime, g_-2]'
    f1 = V.clean({(b, ()): s for b, s in zip(n1, sol)})
    return {'H1': V.clean(H1), 'e1': e1, 'f1': f1, 'H2': H2, 'e2': e2, 'f2': f2, 'c1': c1, 'c2': [x + y for x, y in zip(c2, cs)]}

def check(T, levels):
    out = []
    for (H, e, f), k in (((T['H1'], T['e1'], T['f1']), levels[0]), ((T['H2'], T['e2'], T['f2']), levels[1])):
        out.append(br(H, e) == {kk: 2*v for kk, v in e.items()})
        out.append(br(H, f) == {kk: -2*v for kk, v in f.items()})
        out.append(br(e, f) == H)
        out.append(V.mode_vec(H, 1, H) == {V.VAC: mpq(2*k)})
        out.append(V.mode_vec(e, 1, f) == {V.VAC: mpq(k)})
    # the two ideals commute
    for x in ('H1', 'e1', 'f1'):
        for y in ('H2', 'e2', 'f2'):
            out.append(not br(T[x], T[y]))
            out.append(not V.clean(V.mode_vec(T[x], 1, T[y])))
    return all(out), out

def weight_of_key(T, key):
    b = key[0]
    return (int(hval(T['c1'], b)), int(hval(T['c2'], b)))

if __name__ == '__main__':
    for case in ('I', 'II'):
        T = build(case)
        ok, _ = check(T, CASES[case]['levels'])
        print('case', case, 'levels', CASES[case]['levels'], 'sl2 + sl2 relations, levels, commutation:', ok)
