"""Exact weight-two relations by the obstruction principle (Axis A).

W in M_2(rho, mu): holomorphic weight-two form for rho = rho_8 with W_z in q^{mu_z} C[[q]],
mu_z = Lam_z - 1 for z in P (P = {0} if gamma = 0, {0, d3} if gamma = 1), mu_z = Lam_z otherwise.
G weakly holomorphic weight-zero form for rho^* (built exactly in vvmf_exact).  Then sum_z W_z G_z
is a weakly holomorphic scalar form of weight two for SL_2(Z), hence has constant term 0.
If G_z vanishes to the right order, this constant term only involves the 13 coefficients
   w0 = [q^{-1/3}] W_0, w0p = [q^{2/3}] W_0, w_z = [q^{mu_z}] W_z (z != 0),
and gives an EXACT linear relation.  We produce the space of all relations obtained from
G = pi^* (eta^{-8} Y), Y in the span of
   E4 Xi_j J^m,  serre_2(D(Xi_j J^m)),  DJ * D(Xi_j J^m),  DJ^2 Xi_j J^m,  (j < 12, m <= MJ),
and compare it with the relations of v02 (weight2_relations_gamma*.json, obtained numerically).
Second part: exact elements of M_2(rho, mu), from D(Xi_j J^m) and DJ * Xi_j J^m, and the rank of
the 13 functionals on them.  Equality of the two bounds gives the exact dimension of the image.
Usage: python3 weight2_obstruction_exact.py gamma [MJ]
"""
import sys, json, time
from fractions import Fraction as Fr
import vvmf_exact as V

t0 = time.time()
gamma = int(sys.argv[1]) if len(sys.argv) > 1 else 1
MJ = int(sys.argv[2]) if len(sys.argv) > 2 else 2
d, Lam, names = V.d, V.Lam, V.names
D3 = names.index('d3')
P = [0, D3] if gamma == 1 else [0]
labels = ['w0', 'w0p'] + ['w_' + names[z] for z in range(1, d)]
assert V.check_pi_exact()
print('gamma =', gamma, ' MJ =', MJ, ' pi =', [names[V.PI[z]] for z in range(d)])

# ---------- Part 1: relations ----------
gens = []
for j in range(d):
    for m in range(MJ + 1):
        B = V.scal(V.Jpow(m), V.col(j))
        DB = V.D(B)
        gens.append(('E4*Xi_%s*J^%d' % (names[j], m), V.scal(V.E4, B)))
        gens.append(('serre(D(Xi_%s*J^%d))' % (names[j], m), V.serre(DB, 2)))
        gens.append(('DJ*D(Xi_%s*J^%d)' % (names[j], m), V.scal(V.DJ, DB)))
        gens.append(('DJ^2*Xi_%s*J^%d' % (names[j], m), V.scal(V.smul(V.DJ, V.DJ), B)))
Hs = [V.H_from_Y(Y) for _, Y in gens]
# conditions on G_z = H_{pi z}: exponent of G_z is -Lam_z + k_z + n (n = index in H_{pi z})
# pairing coefficient of W_z at q^{Lam_z + m} is G_z at n = -k_z - m
#   z not in P: only m = 0 allowed -> kill n < -k_z ; relation coefficient at n = -k_z
#   z = 0: m in {-1 (w0), 0 (w0p)} -> kill n < -k_0 ; coefficients at n = -k_0 + 1 (w0), -k_0 (w0p)
#   z = d3 in P (gamma=1): only m = -1 -> kill n <= -k_z ; coefficient at n = -k_z + 1
kill = []          # list of (z, n)
coef = {}          # label -> (z, n)
for z in range(d):
    k = -V.kz(z)
    if z == 0:
        coef['w0'] = (0, k + 1); coef['w0p'] = (0, k)
        lim = k - 1
    elif z in P:
        coef['w_' + names[z]] = (z, k + 1)
        lim = k
    else:
        coef['w_' + names[z]] = (z, k)
        lim = k - 1
    # all n <= lim must vanish; the H's have finitely many such n
    nmin = min([min(H[V.PI[z]].keys()) for H in Hs if H[V.PI[z]]] + [lim + 1])
    for n in range(nmin, lim + 1):
        kill.append((z, n))
ng = len(gens)
rows = [[H[V.PI[z]].get(n, Fr(0)) for H in Hs] for (z, n) in kill]
null = V.nullspace(rows, ng)
print('generators:', ng, ' vanishing conditions:', len(kill), ' admissible G (dim):', len(null), ' (%.1fs)' % (time.time() - t0))
rels = []
for vec in null:
    rel = []
    for lab in labels:
        z, n = coef[lab]
        rel.append(sum(c*H[V.PI[z]].get(n, Fr(0)) for c, H in zip(vec, Hs) if c != 0))
    rels.append(rel)
R, piv = V.rref(rels, len(labels))
print('independent exact relations obtained:', len(R))
for row in R:
    print('  0 = ' + ' + '.join('(%s)*%s' % (c, labels[j]) for j, c in enumerate(row) if c != 0))
# comparison with v02
old = json.load(open('../gamma1/weight2_relations_gamma%d.json' % gamma))
assert old['labels'] == labels
O = [[Fr(c) for c in row] for row in old['relations']]
rk_old = V.rank(O, len(labels)); rk_new = len(R); rk_join = V.rank(O + [list(r) for r in R], len(labels))
print('rank v02 relations =', rk_old, ' rank exact relations =', rk_new, ' rank of union =', rk_join)
same = (rk_old == rk_new == rk_join)
print('=> the exact relations', 'COINCIDE with' if same else 'DIFFER from', 'the relations of v02')
json.dump({'gamma': gamma, 'labels': labels, 'relations': [[str(c) for c in row] for row in R],
           'method': 'obstruction principle, exact', 'MJ': MJ},
          open('weight2_relations_exact_gamma%d.json' % gamma, 'w'), indent=1)

# certificate: for each relation, one explicit G (coefficients on the generators)
cert = []
for vec in null:
    cert.append({gens[i][0]: str(c) for i, c in enumerate(vec) if c != 0})
json.dump({'gamma': gamma, 'MJ': MJ, 'G_basis_on_generators': cert}, open('weight2_obstruction_certificate_gamma%d.json' % gamma, 'w'), indent=0)

# ---------- Part 2: exact elements of M_2(rho, mu) ----------
fgens = []
for j in range(d):
    for m in range(MJ + 1):
        B = V.scal(V.Jpow(m), V.col(j))
        fgens.append(('D(Xi_%s*J^%d)' % (names[j], m), V.D(B)))
        fgens.append(('DJ*Xi_%s*J^%d' % (names[j], m), V.scal(V.DJ, B)))
fk = []
for z in range(d):
    lo = -1 if z in P else 0
    nmin = min([min(Wf[z].keys()) for _, Wf in fgens if Wf[z]] + [lo])
    for n in range(nmin, lo):
        fk.append((z, n))
rows2 = [[Wf[z].get(n, Fr(0)) for _, Wf in fgens] for (z, n) in fk]
null2 = V.nullspace(rows2, len(fgens))
fun = {'w0': (0, -1), 'w0p': (0, 0)}
for z in range(1, d):
    fun['w_' + names[z]] = (z, -1 if z in P else 0)
img = []
forms = []
for vec in null2:
    Wv = V.zero()
    for c, (_, Wf) in zip(vec, fgens):
        if c != 0:
            Wv = V.vadd(Wv, Wf, 1, c)
    forms.append(Wv)
    img.append([Wv[fun[l][0]].get(fun[l][1], Fr(0)) for l in labels])
rk_img = V.rank(img, len(labels))
print('exact elements of M_2(rho,mu) found (dim of span):', len(null2), ' rank of the 13 functionals on them:', rk_img)
# all exact relations must vanish on these forms
bad = sum(1 for row in R for v in img if sum(a*b for a, b in zip(row, v)) != 0)
print('relations evaluated on the exact forms: failures =', bad)
print('=> exact: rank of the functionals on M_2(rho,mu) is', rk_img if rk_img + len(R) == len(labels) else 'NOT determined',
      '(lower bound from forms %d, upper bound from relations %d)' % (rk_img, len(labels) - len(R)))
print('total time %.1fs' % (time.time() - t0))

# ---------- Part 3: the 13 functionals are injective on M_2(rho, mu) ----------
# Let W be in M_2(rho,mu) with all 13 functionals zero.  Then F = W E10/Delta = Xi(a + J b) (Theorem EG)
# has b = 0 and a_z = 0 for z != d3, hence W = a_{d3} Xi_{d3} Delta/E10 (gamma = 1), and W = 0 (gamma = 0).
# For gamma = 1, W_{d3} = a_{d3} q^{Lam_{d3}} + ..., W_0 = O(q^{5/3}), W_z = O(q^{Lam_z + 1}) otherwise.
# A form G for rho^* killing all these coefficients except [q^{Lam_{d3}}] W_{d3} and pairing nontrivially
# with it gives a_{d3} = 0.
if gamma == 1:
    kill3 = []
    for z in range(d):
        k = -V.kz(z)
        lim = (k - 1) if z == D3 else (k - 1 if z != 0 else k - 1)
        # z != d3: W_z[m] may be nonzero for m >= 1 only -> kill n <= -k_z - 1 ; z = d3: m >= 0 -> kill n < -k_z
        lim = k - 1
        nmin = min([min(H[V.PI[z]].keys()) for H in Hs if H[V.PI[z]]] + [lim + 1])
        for n in range(nmin, lim + 1):
            kill3.append((z, n))
    rows3 = [[H[V.PI[z]].get(n, Fr(0)) for H in Hs] for (z, n) in kill3]
    null3 = V.nullspace(rows3, ng)
    kd3 = -V.kz(D3)
    vals = [sum(c*H[V.PI[D3]].get(kd3, Fr(0)) for c, H in zip(vec, Hs) if c != 0) for vec in null3]
    ok = any(v != 0 for v in vals)
    print('Part 3 (gamma = 1): admissible G:', len(null3), '; some G pairs nontrivially with [q^Lam_d3] W_d3:', ok)
    print('=> dim M_2(rho, mu) =', rk_img if ok and rk_img + len(R) == len(labels) else 'undetermined')
else:
    print('Part 3 (gamma = 0): the functionals are injective by the structure argument; dim M_2(rho, mu) =',
          rk_img if rk_img + len(R) == len(labels) else 'undetermined')
