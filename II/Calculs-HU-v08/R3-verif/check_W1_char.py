import sys, json, itertools, os
sys.path.insert(0, '.')
from lat import *
from collections import Counter
mult = Counter()
for i in range(16):
    mult[(0, 0, 0, 0)] += 1
for a in ROOTS:
    mult[tuple(ip(a, Hvec(X)) for X in range(4))] += 1
dec = {}
for t in list(mult):
    if all(x >= 0 for x in t):
        m = 0
        for e in itertools.product((0, 1), repeat=4):
            m += (-1)**sum(e) * mult.get(tuple(t[i]+2*e[i] for i in range(4)), 0)
        if m:
            dec[t] = m
exp = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 't12_W1.json')))['W1']
exp = {tuple(int(x) for x in k.split(',')): v for k, v in exp.items()}
print(sorted(dec.items()))
print('n types', len(dec), 'dim', sum(m*__import__('math').prod(x+1 for x in t) for t, m in dec.items()))
print('MATCH expected W1 decomposition:', dec == exp)
# per copy
for cp in range(2):
    mc = Counter()
    for i in range(8): mc[(0,0,0,0)] += 1
    for a in (ROOTS[:240] if cp == 0 else ROOTS[240:]):
        mc[tuple(ip(a, Hvec(X)) for X in range(4))] += 1
    d = {}
    for t in list(mc):
        if all(x >= 0 for x in t):
            m = sum((-1)**sum(e)*mc.get(tuple(t[i]+2*e[i] for i in range(4)), 0) for e in itertools.product((0,1), repeat=4))
            if m: d[t] = m
    print('copy', cp+1, sorted(d.items()))
