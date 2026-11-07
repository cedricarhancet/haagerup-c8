from collections import defaultdict
from fractions import Fraction as Fr
R = {'c1': {(0,2):1}, 'c2': {(0,0):1,(0,4):1,(1,1):1,(1,3):1,(2,2):1}, 'd1': {(0,2):1,(1,1):1},
     'd5': {(0,4):1,(1,3):1}, 'd4': {(0,0):1,(1,1):1}, 'd6': {(0,2):1,(1,3):1,(2,0):1}}
PI = {'c1':'c2','c2':'c1','d1':'d5','d5':'d1','d4':'d6','d6':'d4'}
HZ = {'c1': Fr(1,3), 'c2': Fr(2,3), 'd1': Fr(6,13), 'd5': Fr(7,13), 'd4': Fr(5,13), 'd6': Fr(8,13)}
hA = lambda a,b: Fr(a*(a+2),20) + Fr(b*(b+2),60)
comp = defaultdict(list)
comp[(2,0,0,0)].append('S1x1'); comp[(0,2,0,0)].append('S1x1'); comp[(0,0,2,0)].append('1xS1'); comp[(0,0,0,2)].append('1xS1')
for z, pz in PI.items():
    for (a,b) in R[z]:
        for (c,d) in R[pz]:
            comp[(a,b,c,d)].append(z)
print('number of components:', len(comp))
for k in sorted(comp):
    secs = comp[k]
    lam = [0 if s in ('S1x1','1xS1') else HZ[s]-hA(k[0],k[1]) for s in secs]
    print(k, secs, [str(x) for x in lam])
print('max non-S sectors:', max(sum(1 for s in v if s not in ('S1x1','1xS1')) for v in comp.values()))
