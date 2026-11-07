# sanity tests of voa.py: Virasoro, commutator formula, skew-symmetry, on random states
import sys, random
sys.path.insert(0, '.')
from voa import *
from lat import ROOTS, NORM4
random.seed(1)
omega = {}
for i in range(16):
    omega[(((1, i), (1, i)), ZERO)] = Fr(1, 2)
def L(n, st):
    return mode_vec(omega, n+1, st)
# test: L(0) on weight-2 states, L(-1)=translation, c=16
c = mode_vec(omega, 3, omega)
print('omega_(3)omega =', c, ' (expect c/2 = 8)')
print('L0 omega = 2 omega:', L(0, omega) == {k: 2*v for k, v in omega.items()})
tests = [e_state(random.choice(ROOTS)), e_state(random.choice(NORM4)),
         {(((1, 3),), random.choice(ROOTS)): Fr(1)}, {(((2, 5),), ZERO): Fr(1)},
         {(((1, 2), (1, 7)), ZERO): Fr(1)}]
for t in tests:
    k = list(t)[0]
    print('L0 eigen', state_wt(k), L(0, t) == {kk: state_wt(k)*v for kk, v in t.items()})
# skew symmetry  a_(n) b = sum_j (-1)^{n+j+1} L(-1)^j/j! b_(n+j) a
from math import factorial
def Lm1pow(j, st):
    for _ in range(j):
        st = L(-1, st)
    return st
ok = True
for trial in range(30):
    a = random.choice(tests + [e_state(random.choice(ROOTS)) for _ in range(3)])
    b = random.choice(tests + [e_state(random.choice(ROOTS)) for _ in range(3)])
    # make lattice parts interact
    for n in range(0, 3):
        lhs = mode_vec(a, n, b)
        rhs = {}
        for j in range(0, 6):
            t = mode_vec(b, n+j, a)
            if t:
                sadd(rhs, Lm1pow(j, t), Fr((-1)**(n+j+1), factorial(j)))
        if clean(lhs) != clean(rhs):
            ok = False; print('skew fail', a, b, n)
print('skew symmetry tests ok:', ok)
# commutator formula [a_(m), b_(n)] v = sum_j binom(m,j) (a_(j) b)_(m+n-j) v, with a,b weight1, v weight 2
ok = True
for trial in range(40):
    a = e_state(random.choice(ROOTS)) if random.random() < 0.7 else h_state(random.randrange(16))
    # choose b with <a,b> = -1 often
    b = e_state(random.choice(ROOTS))
    v = random.choice(tests)
    for m in range(0, 2):
        for n in range(-1, 2):
            lhs = sadd(mode_vec(a, m, mode_vec(b, n, v)), mode_vec(b, n, mode_vec(a, m, v)), -1)
            rhs = {}
            for j in range(0, 3):
                ab = mode_vec(a, j, b)
                if ab:
                    sadd(rhs, mode_vec(ab, m+n-j, v), comb(m, j))
            if clean(lhs) != clean(rhs):
                ok = False; print('comm fail', a, b, v, m, n)
print('commutator tests ok:', ok)
# a test with a,b of weight 2 composite: [a_(m), b_(n)] with a = h h, b = e^beta
ok = True
for trial in range(15):
    a = random.choice(tests)
    b = random.choice(tests)
    v = random.choice(tests[:2] + [e_state(random.choice(ROOTS))])
    for m in range(0, 3):
        n = random.randrange(-1, 3)
        lhs = sadd(mode_vec(a, m, mode_vec(b, n, v)), mode_vec(b, n, mode_vec(a, m, v)), -1)
        rhs = {}
        for j in range(0, 6):
            ab = mode_vec(a, j, b)
            if ab:
                sadd(rhs, mode_vec(ab, m+n-j, v), comb(m, j))
        if clean(lhs) != clean(rhs):
            ok = False; print('comm2 fail', list(a), list(b), list(v), m, n)
print('commutator tests (weight 2) ok:', ok)
