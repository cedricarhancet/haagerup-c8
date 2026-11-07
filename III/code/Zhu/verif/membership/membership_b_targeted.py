"""Verify J_K <= I_K using only the overlaps needed to span the target rules.

The input G0 is the set already certified by membership_a.py. Reducing a
target by G0 preserves its class modulo the ideal of G0. Every overlap
remainder used below belongs to that ideal. Echelon steps and the final
span tests use exact arithmetic in K.
"""
from collections import defaultdict
from functools import lru_cache
import pickle
import sys
import time

sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu-exact')
import kfast as field

ROOT = '/home/claude/hj/work/Calculs/zhu-exact'
with open(ROOT + '/upper_exact.pkl', 'rb') as f:
    rules = pickle.load(f)['rules']
with open(ROOT + '/rows_exact.pkl', 'rb') as f:
    arrows = pickle.load(f)['ARROWS']
with open('/home/claude/hj/work/verif/membership/G0.pkl', 'rb') as f:
    good = {q: rules[q] for q in pickle.load(f)['G0']}
assert len(good) == 5105 and all(len(q) == 2 for q in good)


def source(a):
    return arrows[a]['source']


def target(a):
    return arrows[a]['target']


def endpoints(word):
    return target(word[0]), source(word[-1])


def order(word):
    return -len(word), word


def add(out, row, coeff):
    for word, value in row.items():
        term = value if coeff == field.ONE else field.mul(coeff, value)
        if word in out:
            term = field.add(out[word], term)
        if field.iszero(term):
            out.pop(word, None)
        else:
            out[word] = term


@lru_cache(maxsize=30000)
def normal(word):
    for i in range(len(word) - 2, -1, -1):
        q = word[i:i + 2]
        if q in good:
            out = {}
            for tail, coeff in good[q].items():
                add(out, normal(word[:i] + tail + word[i + 2:]), coeff)
            return out
    return {word: field.ONE}


def reduce_by_table(row, table):
    while row:
        lead = min(row, key=order)
        if lead not in table:
            break
        add(row, table[lead], field.neg(row[lead]))
    return row


def overlaps(start, end, by_end, by_first):
    for q in by_end[start]:
        for r in by_first[q[-1]]:
            if source(r[-1]) != end:
                continue
            delta = {}
            for w, c in good[q].items():
                add(delta, normal(w + (r[-1],)), c)
            for w, c in good[r].items():
                add(delta, normal((q[0],) + w), field.neg(c))
            yield delta


t0 = time.monotonic()
targets_by_endpoint = defaultdict(list)
by_first = defaultdict(list)
by_end = defaultdict(list)
for q in good:
    by_first[q[0]].append(q)
    by_end[target(q[0])].append(q)
for q, rhs in rules.items():
    if q in good:
        continue
    row = dict(normal(q))
    for w, c in rhs.items():
        add(row, normal(w), field.neg(c))
    targets_by_endpoint[endpoints(q)].append(row)

tested = 0
overlap_count = 0
for ep, goals in sorted(targets_by_endpoint.items()):
    table = {}
    remaining = [g for g in goals if g]
    for row in overlaps(*ep, by_end, by_first):
        if not remaining:
            break
        overlap_count += 1
        if not reduce_by_table(row, table):
            continue
        lead = min(row, key=order)
        inv = field.inv(row[lead])
        table[lead] = {w: field.mul(c, inv) for w, c in row.items()}
        remaining = [g for g in remaining if reduce_by_table(g, table)]
    if remaining:
        raise RuntimeError(f"Could not certify {len(remaining)} rules in corner {ep}")
    tested += len(goals)
    print(f"Corner {ep}, {len(goals)} certified, cumulative {tested}, {time.monotonic()-t0:.0f}s", flush=True)

assert tested == 554
print(f"Completion rules certified: {tested} of 554 failures 0; overlaps used {overlap_count}", flush=True)
