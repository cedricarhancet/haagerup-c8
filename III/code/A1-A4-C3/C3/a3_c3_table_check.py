#!/usr/bin/env python3
"""Check the exact input tables for Block C3 of Article III.

Compare the twelve lowest weight spaces and conformal weights with Table 2
and the lowest-spaces equation in Article I, with the block table in the Zhu
appendix of Article III, and with the corresponding data in the extracted
Certificats-Zhu-III-v05 archive. Recompute the W_1 character from the tables.

Only the Python standard library is needed. An optional --out writes a JSON
report. Input hashes in the report are keyed by their roles in the check.
"""
import argparse
from collections import Counter
from fractions import Fraction as Fr
import hashlib
import json
import math
from pathlib import Path
import pickle
import re
import sys

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--calc", type=Path, default=Path(__file__).resolve().parent / "run",
                help="C3 calculation directory containing data/t12_table.json and data/t12_W1.json")
ap.add_argument("--artI", type=Path, required=True, help="LaTeX source of Article I")
ap.add_argument("--artIII", type=Path, required=True, help="LaTeX source of Article III")
ap.add_argument("--zhu", type=Path, required=True, help="root of the extracted Certificats-Zhu-III-v05 archive")
ap.add_argument("--out", type=Path, help="optional path for the JSON report")
args = ap.parse_args()

LABELS = ['0', 'b', 'a', 'c0', 'c1', 'c2', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6']
report = {'inputs': {}, 'checks': [], 'discrepancies': []}
def record_input(role, path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    report["inputs"][role] = digest
    print(f"input {role}  sha256 {digest}")
def note(name, ok, detail=''):
    report['checks'].append({'check': name, 'ok': bool(ok), 'detail': detail})
    print(('OK   ' if ok else 'ECART') + ' ' + name + ((' : ' + detail) if detail else ''))
    if not ok:
        report['discrepancies'].append({'check': name, 'detail': detail})
def dimV(mod):
    return sum((a + 1) * (b + 1) * m for (a, b), m in mod.items())
def fmt(mod):
    return ' + '.join(('%d' % m if m > 1 else '') + '(%d,%d)' % k for k, m in sorted(mod.items()))

# ------------------------------------------------------------------ (T12) data/t12_table.json, data/t12_W1.json
p_t12 = args.calc / "data" / "t12_table.json"
p_w1 = args.calc / "data" / "t12_W1.json"
record_input("data/t12_table.json", p_t12)
record_input("data/t12_W1.json", p_w1)
record_input("Article I", args.artI)
record_input("Article III", args.artIII)
T12 = json.loads(p_t12.read_text(encoding="utf-8"))
W1J = json.loads(p_w1.read_text(encoding="utf-8"))
def parse_mod_json(d):
    return {tuple(map(int, k.split(','))): int(v) for k, v in d.items()}
t12_R = {z: parse_mod_json(T12[z]['R']) for z in LABELS}
t12_h = {z: Fr(T12[z]['h']) for z in LABELS}
t12_dim = {z: int(T12[z]['dimR']) for z in LABELS}
t12_U2 = parse_mod_json(T12['U2'])
note('T12: dimR coherent avec la decomposition R_z (12 labels)', all(dimV(t12_R[z]) == t12_dim[z] for z in LABELS),
     str({z: (dimV(t12_R[z]), t12_dim[z]) for z in LABELS if dimV(t12_R[z]) != t12_dim[z]}) or 'toutes')
note('T12: dim U2 = 120', dimV(t12_U2) == 120, str(dimV(t12_U2)))
note('T12: bornes des niveaux (a <= 3, b <= 13) sur les R_z', all(a <= 3 and b <= 13 for z in LABELS for (a, b) in t12_R[z]))

# ------------------------------------------------------------------ (A-I) Article I : Table 2 (case II) et eq:lowest-spaces
src = args.artI.read_text(encoding="utf-8")
i_lab = src.index(r'\label{tab:lowest-spaces}')
i_beg = src.rindex(r'\begin{tabular}', 0, i_lab)
i_end = src.index(r'\end{tabular}', i_beg)
tab = src[i_beg:i_end]
# sanity: header names the two cases, column 3 is case (II), levels (3,13)
hdr = [l for l in tab.split('\n') if 'Case (I)' in l][0]
cols = [c.strip() for c in hdr.split('&')]
note('A-I Table 2: en-tete « Case (I), levels (1,39) | Case (II), levels (3,13) » (colonne 3 = cas II)',
     len(cols) == 3 and '(1,39)' in cols[1] and '(3,13)' in cols[2], hdr.strip())
rows = {}; cur = None
for line in tab.split('\n'):
    if '&' not in line or 'Case (I)' in line:
        continue
    parts = line.split('&')
    if len(parts) != 3:
        continue
    head, cI, cII = (p.strip() for p in parts)
    if head:
        cur = head
        rows[cur] = [cII]
    elif cur is not None:
        rows[cur].append(cII)
def parse_mod_tex(cells):
    txt = ' '.join(cells).replace('\\\\', ' ').replace('$', ' ').replace('\\ ', ' ').replace('\\hline', ' ')
    txt = txt.replace(' ', '')
    # the cell must be a sum of terms m(a,b) only
    assert re.fullmatch(r'\+?(\d*\(\d+,\d+\)\+?)+', txt), 'cellule inattendue: ' + txt
    d = Counter()
    for mult, a, b in re.findall(r'(\d*)\((\d+),(\d+)\)', txt):
        d[(int(a), int(b))] += int(mult) if mult else 1
    return dict(d)
artI_R = {}; artI_U2 = None
for head, cells in rows.items():
    mod = parse_mod_tex(cells)
    if head == '$U_2$':
        artI_U2 = mod; continue
    for nm in re.findall(r'R_\{?([a-z0-9_]+)\}?', head):
        artI_R[nm.replace('_', '')] = mod
note('A-I Table 2: onze lignes R_z lues (b, a, c0, c1, c2, d1..d6) et la ligne U_2', sorted(artI_R) == sorted(LABELS[1:]) and artI_U2 is not None,
     'lues: ' + ', '.join(sorted(artI_R)))
artI_R['0'] = {(0, 0): 1}     # R_0 = C 1 (not a row of Table 2); dim R_0 = 1 by eq:lowest-spaces below
# eq:lowest-spaces with gamma = 0 (epsilon = 1)
blk = src[src.index(r'\label{eq:lowest-spaces}'):]
blk = blk[:blk.index(r'\end{equation}')]
def tex_frac_eval(s, gamma=Fr(0), epsi=Fr(1)):
    s = s.strip().replace('$', '')
    s = re.sub(r'\\frac\{(\d+)\}\{(\d+)\}', r'(\1/\2)', s)
    s = re.sub(r'\\frac(\d)(\d)', r'(\1/\2)', s)
    s = s.replace(r'\gamma', 'g').replace(r'\epsi', 'e').replace('(1-g)', '*(1-g)')
    s = re.sub(r'(\d)\s*g', r'\1*g', s).replace('g+', 'g+')
    # safe evaluation: only digits, + - * / ( ) g e
    assert re.fullmatch(r'[0-9+\-*/() ge]+', s), 'expression inattendue: ' + s
    import ast
    def ev(n):
        if isinstance(n, ast.Expression): return ev(n.body)
        if isinstance(n, ast.Constant): return Fr(n.value)
        if isinstance(n, ast.Name): return {'g': gamma, 'e': epsi}[n.id]
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            return {ast.Add: a + b, ast.Sub: a - b, ast.Mult: a * b, ast.Div: a / b if b else None}[type(n.op)]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return -ev(n.operand)
        raise ValueError(ast.dump(n))
    return ev(ast.parse(s, mode='eval'))
eq_h = {}; eq_dim = {}
lines = [l.strip() for l in blk.split('\n') if '&' in l]
zrows = [l for l in lines if l.startswith('z &')]
hrows = [l for l in lines if l.startswith('h_z &')]
drows = [l for l in lines if l.startswith(r'\dim R_z &')]
assert len(zrows) == len(hrows) == len(drows) == 2, (zrows, hrows, drows)
def tex_label(s):
    s = s.strip().rstrip('\\').strip().replace('$', '')
    return s.replace('_', '').replace('{', '').replace('}', '')
for zr, hr, dr in zip(zrows, hrows, drows):
    zs = [tex_label(x) for x in zr.split('&')[1:]]
    hs = [tex_frac_eval(x.rstrip('\\')) for x in hr.split('&')[1:]]
    ds = [tex_frac_eval(x.rstrip('\\')) for x in dr.split('&')[1:]]
    for z, h, d in zip(zs, hs, ds):
        eq_h[z] = h; eq_dim[z] = d
note('A-I eq:lowest-spaces: douze labels lus', sorted(eq_h) == sorted(LABELS), ', '.join(sorted(eq_h)))
note('A-I eq:lowest-spaces (gamma = 0, epsilon = 1): dim R_z = dim de la ligne de Table 2 (12 labels)',
     all(eq_dim[z] == dimV(artI_R[z]) for z in LABELS), str({z: (str(eq_dim[z]), dimV(artI_R[z])) for z in LABELS}))
note('A-I: dim U_2 (Table 2, cas II) = 120 = eq:lowest-spaces (gamma = 0)', dimV(artI_U2) == 120, str(dimV(artI_U2)))

# ------------------------------------------------------------------ T12 vs A-I
bad = [z for z in LABELS if t12_R[z] != artI_R[z]]
note('T12 = A-I Table 2 (cas II): decompositions R_z (types et multiplicites) pour les 12 labels', not bad,
     'ecarts: ' + str({z: (fmt(t12_R[z]), fmt(artI_R[z])) for z in bad}) if bad else 'les 12 egalites exactes')
bad = [z for z in LABELS if t12_h[z] != eq_h[z]]
note('T12 = A-I eq:lowest-spaces: poids h_z (gamma = 0, epsilon = 1) pour les 12 labels', not bad,
     'ecarts: ' + str({z: (str(t12_h[z]), str(eq_h[z])) for z in bad}) if bad else str({z: str(t12_h[z]) for z in LABELS}))
note('T12 = A-I: U_2 (cas II)', t12_U2 == artI_U2, fmt(t12_U2) if t12_U2 == artI_U2 else '%s  vs  %s' % (fmt(t12_U2), fmt(artI_U2)))

# ------------------------------------------------------------------ (Article III) table zhu:table:blocks du source LaTeX
article_iii = args.artIII.read_text(encoding="utf-8")
j_lab = article_iii.index(r'\label{zhu:table:blocks}')
j_beg = article_iii.rindex(r'\begin{tabular}', 0, j_lab)
j_end = article_iii.index(r'\end{tabular}', j_beg)
ptab = article_iii[j_beg:j_end]
article_iii_rows = {}    # label -> (n_i, module, dim, h)
for line in ptab.split('\n'):
    if '&' not in line or 'label' in line:
        continue
    parts = [p.strip() for p in line.split('&')]
    assert len(parts) == 5, line
    labs = [tex_label(x) for x in parts[0].split(',')]
    n_i = int(parts[1].replace('$', ''))
    mod = parse_mod_tex([parts[2]])
    dim = int(parts[3].replace('$', ''))
    h = tex_frac_eval(parts[4].rstrip('\\'))
    for lb in labs:
        article_iii_rows[lb] = (n_i, mod, dim, h)
note('Article III zhu:table:blocks: douze labels lus', sorted(article_iii_rows) == sorted(LABELS), ', '.join(sorted(article_iii_rows)))
note('Article III zhu:table:blocks: colonne dim = dim de la decomposition', all(article_iii_rows[z][2] == dimV(article_iii_rows[z][1]) for z in LABELS))
note('Article III zhu:table:blocks: n_i = somme des multiplicites m_{i,v} (dimension du module simple du coin D)',
     all(article_iii_rows[z][0] == sum(article_iii_rows[z][1].values()) for z in LABELS),
     str({z: (article_iii_rows[z][0], sum(article_iii_rows[z][1].values())) for z in LABELS}))
bad = [z for z in LABELS if article_iii_rows[z][1] != t12_R[z]]
note('T12 = Article III zhu:table:blocks: decompositions R_z pour les 12 labels', not bad,
     'ecarts: ' + str({z: (fmt(t12_R[z]), fmt(article_iii_rows[z][1])) for z in bad}) if bad else 'les 12 egalites exactes')
bad = [z for z in LABELS if article_iii_rows[z][3] != t12_h[z]]
note('T12 = Article III zhu:table:blocks: poids h_i pour les 12 labels', not bad,
     'ecarts: ' + str({z: (str(t12_h[z]), str(article_iii_rows[z][3])) for z in bad}) if bad else 'les 12 egalites exactes')
# the vertex list lambda_0..lambda_11 of Article III
m = re.search(r'\\lambda_0, \\ldots, \\lambda_\{11\}\$ are the twelve pairs ([^.]*)\.', article_iii)
lam_text = [tuple(map(int, x)) for x in re.findall(r'\((\d+),(\d+)\)', m.group(1))] if m else None
note('Article III texte: liste lambda_0..lambda_11 lue (12 paires)', lam_text is not None and len(lam_text) == 12, str(lam_text))

# ------------------------------------------------------------------ (ZHU) certificats de l'Article III
zd = args.zhu / "Calculs" / "zhu-exact"
p_blocks = zd / "specfib" / "specfib_blocks.json"
p_upper = zd / "specfib" / "z4_upper_algebra_modp.json"
p_corner = zd / "rows_exact_reduced_modp.pkl"
p_exactZ = args.zhu / "reference" / "exactZ.json"
for p in (p_blocks, p_upper, p_corner, p_exactZ):
    record_input(p.relative_to(args.zhu).as_posix(), p)
class RU(pickle.Unpickler):      # restricted unpickler: builtins only (the pickle holds ints, tuples, lists, dicts)
    def find_class(self, mod, name):
        raise pickle.UnpicklingError('forbidden global %s.%s' % (mod, name))
with p_corner.open("rb") as stream:
    CR = RU(stream).load()
CONS, VERT = [tuple(c) for c in CR['CONS']], list(CR['VERT'])
# vertex -> list of constituents (types); vertices 0..11 have one constituent, 12 and 13 two (killed in Lambda)
v2types = {}
for c, v in zip(CONS, VERT):
    v2types.setdefault(v, []).append(c)
note('ZHU rows_exact_reduced_modp.pkl: 16 constituants, sommets 0..11 simples, 12 et 13 doubles',
     len(CONS) == 16 and all(len(v2types[v]) == 1 for v in range(12)) and len(v2types[12]) == 2 and len(v2types[13]) == 2,
     'v -> types: ' + str({v: v2types[v] for v in sorted(v2types)}))
note('ZHU: types des sommets 0..11 = liste lambda_0..lambda_11 du texte de l’Article III',
     lam_text is not None and [v2types[v][0] for v in range(12)] == lam_text)
def ranks_to_mod(ranks):
    """vertex ranks {v: m} -> g-module {(r,l): m} (vertices 0..11 only; 12, 13 must be absent)"""
    mod = Counter()
    for v, mlt in ranks.items():
        v = int(v)
        assert v <= 11, 'sommet tue %d avec rang non nul' % v
        mod[v2types[v][0]] += int(mlt)
    return dict(mod)
SB = json.loads(p_blocks.read_text(encoding="utf-8"))
h2z = {}
for z in LABELS:
    h2z.setdefault(t12_h[z], []).append(z)
roots_seen = {}
all_ok = True; details = []
for ent in SB:
    c = Fr(ent['rational'])
    roots_seen[c] = ent
    zs = h2z.get(c, [])
    if not zs:
        all_ok = False; details.append('racine %s sans label' % c); continue
    if len(zs) == 1:
        z = zs[0]
        # single block: corner dims are the squares of the vertex ranks
        sq = {int(v): d for v, d in ent['corner_dims'].items()}
        rk = {}
        for v, d in sq.items():
            r = math.isqrt(d)
            if r * r != d:
                all_ok = False; details.append('c=%s: dimension de coin %d non carree' % (c, d))
            rk[v] = r
        mod = ranks_to_mod(rk)
        n_i = sum(rk.values())
        ok = (mod == t12_R[z]) and (ent['kernel_dim'] == n_i ** 2) and (ent['multiplicity_charpoly'] == n_i ** 2) and n_i == article_iii_rows[z][0]
        details.append('c=%s -> %s : module %s, n_i=%d, dim ker=%d : %s' % (c, z, fmt(mod), n_i, ent['kernel_dim'], 'egal' if ok else 'ECART (T12: %s)' % fmt(t12_R[z])))
        all_ok &= ok
    else:
        # c = 1: three blocks b, a, c0 ; the corner dimensions are sums of squares of the vertex ranks
        exp = Counter()
        for z in zs:
            for (r, l), mlt in t12_R[z].items():
                v = [vv for vv in range(12) if v2types[vv][0] == (r, l)][0]
                exp[v] += mlt ** 2
        got = {int(v): d for v, d in ent['corner_dims'].items()}
        n2 = sum(article_iii_rows[z][0] ** 2 for z in zs)
        ok = dict(exp) == got and ent['kernel_dim'] == n2 and ent['multiplicity_charpoly'] == n2
        details.append('c=%s -> %s : sommes de carres des rangs %s, dim ker=%d=%s : %s' % (c, '+'.join(zs), got, ent['kernel_dim'], '+'.join('%d^2' % article_iii_rows[z][0] for z in zs), 'egal' if ok else 'ECART (attendu %s)' % dict(exp)))
        all_ok &= ok
note('ZHU specfib_blocks.json: les 10 racines c de L_Z et les coins e_v ker(L_Z - c) e_v redonnent les R_z de T12 (9 blocs simples, et le triple b+a+c0 pour c = 1)',
     all_ok and set(roots_seen) == set(t12_h.values()), '\n      ' + '\n      '.join(details))
# z4_upper_algebra_modp.json : per-block vertex ranks (12 blocks)
UP = json.loads(p_upper.read_text(encoding="utf-8"))
blocks_cert = sorted((b['size'], tuple(sorted(ranks_to_mod(b['vertex_ranks']).items()))) for b in UP['blocks'])
blocks_t12 = sorted((article_iii_rows[z][0], tuple(sorted(t12_R[z].items()))) for z in LABELS)
note('ZHU z4_upper_algebra_modp.json: multi-ensemble des 12 blocs (taille n_i, rangs de sommets -> types) = les 12 R_z de T12 (a et c0 indistinguables, tous deux de taille 12)',
     blocks_cert == blocks_t12 and UP['dimension'] == 762 and sum(b['size'] ** 2 for b in UP['blocks']) == 762,
     'tailles: ' + str(sorted(b['size'] for b in UP['blocks'])) + ' ; dimension ' + str(UP['dimension']) + ' ; statut du fichier: ' + str(UP.get('status')))
note('ZHU: taille n_i = somme des rangs de sommets dans chaque bloc', all(b['size'] == sum(b['vertex_ranks'].values()) for b in UP['blocks']))
EZ = json.loads(p_exactZ.read_text(encoding="utf-8"))
note('ZHU exactZ.json: ensemble H des racines du polynome annulateur exact de Z = {h_z}, product_vanishes',
     set(Fr(x) for x in EZ['H']) == set(t12_h.values()) and EZ.get('product_vanishes') is True and not EZ.get('absent'),
     'H = ' + str(sorted(EZ['H'], key=Fr)))

# ------------------------------------------------------------------ t12_W1.json vs t12_table.json
pi = W1J['pi']
note('t12_W1.json: pi est une involution de l ensemble des 12 labels', sorted(pi) == sorted(LABELS) and all(pi[pi[z]] == z for z in LABELS), str({z: pi[z] for z in LABELS if pi[z] != z}))
note('t12_W1.json: h_z + h_pi(z) entier pour tout z', all((t12_h[z] + t12_h[pi[z]]).denominator == 1 for z in LABELS),
     str({z: str(t12_h[z] + t12_h[pi[z]]) for z in LABELS}))
w1pairs = [z for z in LABELS if t12_h[z] + t12_h[pi[z]] == 1]
W1 = Counter()
for t in [(2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2)]:
    W1[t] += 1
for z in w1pairs:
    for (a, b), m1 in t12_R[z].items():
        for (c_, d_), m2 in t12_R[pi[z]].items():
            W1[(a, b, c_, d_)] += m1 * m2
W1_file = {tuple(map(int, k.split(','))): v for k, v in W1J['W1'].items()}
dimW1 = sum(m * math.prod(x + 1 for x in t) for t, m in W1.items())
note('t12_W1.json: W_1 = (S_1 x 1) + (1 x S_1) + sum_{h_z + h_pi(z) = 1} R_z x R_pi(z) recalcule a partir de T12 = champ W1 du fichier (20 types), dim 496',
     dict(W1) == W1_file and dimW1 == 496 == W1J['dimW1'], 'paires de poids 1: %s ; dim %d ; types %d' % ([(z, pi[z]) for z in w1pairs], dimW1, len(W1)))
# Dynkin indices (3,13,3,13) x 60 as a control
def ind1(n): return n * (n + 1) * (n + 2) // 6
idx = [0, 0, 0, 0]
for t, m in W1.items():
    for k in range(4):
        idx[k] += m * (math.prod(x + 1 for x in t) // (t[k] + 1)) * ind1(t[k])
note('t12_W1.json: indices de Dynkin de W_1 par facteur = 60 x (3, 13, 3, 13)', idx == [180, 780, 180, 780] and W1J['index_W1'] == idx, str(idx))

# ------------------------------------------------------------------ omega_certificate.json : secteurs et poids utilises par les equations
p_cert = args.calc / "omega_certificate.json"
if p_cert.is_file():
    record_input("omega_certificate.json", p_cert)
    C = json.loads(p_cert.read_text(encoding="utf-8"))
    lw = {z: Fr(v) for z, v in C['lowest_weights'].items()}
    note('certificat: lowest_weights = h_z de T12 pour les six secteurs de poids un (c1, c2, d1, d5, d4, d6)',
         sorted(lw) == sorted(w1pairs) and all(lw[z] == t12_h[z] for z in lw), str({z: str(v) for z, v in lw.items()}))
    # multiset of sectors of each sl2^4-type of W_1, derived from T12 + pi
    sect = {}
    for t in [(2, 0, 0, 0), (0, 2, 0, 0)]: sect.setdefault(t, []).append('S1x1')
    for t in [(0, 0, 2, 0), (0, 0, 0, 2)]: sect.setdefault(t, []).append('1xS1')
    for z in w1pairs:
        for (a, b), m1 in t12_R[z].items():
            for (c_, d_), m2 in t12_R[pi[z]].items():
                sect.setdefault((a, b, c_, d_), []).extend([z] * (m1 * m2))
    cert_sect = {tuple(e['highest_weight']): sorted(e['sectors']) for e in C['multiplicity_spaces']}
    ok = len(cert_sect) == 20 and set(cert_sect) == set(sect) and all(sorted(sect[t]) == cert_sect[t] for t in sect)
    note('certificat: les 20 espaces de multiplicite (type sl2^4, multi-ensemble des secteurs S1x1 / 1xS1 / R_z x R_pi(z)) = ceux deduits de T12 et de pi',
         ok, '' if ok else str({t: (sorted(sect.get(t, [])), cert_sect.get(t)) for t in set(sect) | set(cert_sect) if sorted(sect.get(t, [])) != cert_sect.get(t)}))
    # eigenvalue of o(omega_{C x 1}) on sector z of type (a,b,c,d): h_z - a(a+2)/20 - b(b+2)/60 (0 on S1x1 and 1xS1)
    def eig(z, t):
        if z in ('S1x1', '1xS1'): return Fr(0)
        return t12_h[z] - Fr(t[0] * (t[0] + 2), 20) - Fr(t[1] * (t[1] + 2), 60)
    report['sector_eigenvalues'] = {','.join(map(str, t)): [(z, str(eig(z, t))) for z in sorted(sect[t])] for t in sorted(sect)}
    print('      valeurs propres attendues h_z - h_A(a,b) de o(omega_{C x 1}) par type (a,b,c,d) et secteur :')
    for t in sorted(sect):
        print('        %-14s %s' % (','.join(map(str, t)), ' ; '.join('%s: %s' % (z, eig(z, t)) for z in sorted(sect[t]))))

# ------------------------------------------------------------------ bilan
ok = not report['discrepancies']
report['result'] = 'EGALITE EXACTE des douze decompositions (types, multiplicites) et des douze poids h_z entre T12, Article I (Table 2 cas II + eq:lowest-spaces gamma = 0), Article III zhu:table:blocks et certificats Zhu' if ok else 'ECARTS: ' + str(report['discrepancies'])
report['labels'] = {z: {'h': str(t12_h[z]), 'dim': t12_dim[z], 'R': fmt(t12_R[z]), 'n_i': article_iii_rows[z][0]} for z in LABELS}
if args.out is not None:
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
print('\nRESULTAT :', report['result'])
sys.exit(0 if ok else 1)
