# Comparison with the first computation (R3) OUTPUTS only: restricted unpickler (gmpy2.from_binary only),
# no code from R3 is imported.
import pickle, json, sys, io, os
import gmpy2
from fractions import Fraction as Fr
class RU(pickle.Unpickler):
    def find_class(self, mod, name):
        if (mod, name) == ('gmpy2.gmpy2', 'from_binary'):
            return gmpy2.from_binary
        raise pickle.UnpicklingError('forbidden global %s.%s' % (mod, name))
R3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'R3') + os.sep
Gg = RU(open(R3 + 'G_griess.pkl', 'rb')).load()
print('G_griess.pkl content type:', type(Gg), list(Gg.keys()) if isinstance(Gg, dict) else len(Gg))
for k, v in (Gg.items() if isinstance(Gg, dict) else []):
    print('  ', k, type(v), (len(v) if hasattr(v, '__len__') else v))
form = [[Fr(int(x.numerator), int(x.denominator)) for x in row] for row in Gg['form']]
omCC = [Fr(int(x.numerator), int(x.denominator)) for x in Gg['omCC']]
sols = json.load(open(R3 + 'omegaC1.json'))
# elements of Q(sqrt(-6)) as pairs (a, b) = a + b s, s^2 = -6
def mul(p, q): return (p[0]*q[0] - 6*p[1]*q[1], p[0]*q[1] + p[1]*q[0])
def add(p, q): return (p[0]+q[0], p[1]+q[1])
X = [[(Fr(a), Fr(b)) for a, b in s['coords']] for s in sols]
def bil(x, y):
    acc = (Fr(0), Fr(0))
    for i in range(16):
        for j in range(16):
            if form[i][j]:
                acc = add(acc, mul(mul(x[i], (form[i][j], Fr(0))), y[j]))
    return acc
print('R3 omega_CC . omega_CC =', bil([(c, Fr(0)) for c in omCC], [(c, Fr(0)) for c in omCC]))
print('R3 Gram matrix (u_k, u_l) of the 4 solutions (pairs a + b sqrt(-6)):')
for k in range(4):
    print('   ', [tuple(map(str, bil(X[k], X[l]))) for l in range(4)])
print('R3 (u_k, omega_CC):', [tuple(map(str, bil(X[k], [(c, Fr(0)) for c in omCC]))) for k in range(4)])
