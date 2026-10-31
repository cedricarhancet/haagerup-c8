#!/usr/bin/env python3
"""Exact supplementary certificates for the revision of Sections 5 and 8.
Run from any directory, with --sources pointing to the supplied classification
folder. The imported code computes the Evans--Gannon expansion in rational
arithmetic. No source files or stored computational outputs are overwritten.
"""
import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--sources', type=Path, default=Path(__file__).resolve().parents[1]/'calculs-sources'/'classification')
parser.add_argument('--output', type=Path, default=Path(__file__).with_name('modular-certificates-v02.json'))
args = parser.parse_args()
sys.path.insert(0, str(args.sources.resolve()))
import vvmf_exact as V
import obstruction as O


def enc(matrix):
    return [[str(c) for c in row] for row in matrix]


def describe(k, P, A, mj):
    gens = O.Y_generators(6-k, mj)
    Hs = [V.H_from_Y(y) for _, y in gens]
    positions = []
    for z in range(V.d):
        lower = -1 if z in P else 0
        nmin = min([min(h[V.PI[z]]) for h in Hs if h[V.PI[z]]] + [10**6])
        for m in range(lower, -V.kz(z)-nmin+1):
            if (z,m) not in set(A):
                positions.append((z,-V.kz(z)-m))
    matrix = [[h[V.PI[z]].get(n,Fraction(0)) for h in Hs] for z,n in positions]
    null = V.nullspace(matrix,len(gens))
    evaluation = [[sum(c*h[V.PI[z]].get(-V.kz(z)-m,Fraction(0)) for c,h in zip(vec,Hs)) for z,m in A] for vec in null]
    relations,_ = V.rref(evaluation,len(A))
    form_values = O.exact_forms(k,P,A,mj)
    rank_forms = V.rank(form_values,len(A))
    assert all(sum(c*v for c,v in zip(row,value)) == 0 for row in relations for value in form_values)
    assert len(relations)+rank_forms == len(A)
    return {
        'weight': k, 'max_power_J': mj, 'number_generators': len(gens),
        'vanishing_matrix_shape': [len(matrix),len(gens)],
        'vanishing_matrix_rank': len(gens)-len(null),
        'admissible_coefficient_nullity': len(null),
        'coefficient_positions': [[V.names[z],m] for z,m in A],
        'relations_rank': len(relations), 'evaluation_image_rank': rank_forms,
        'relations_rref': enc(relations),
        'generator_order': [name for name,_ in gens],
        'G_basis_coefficients': [{name:str(c) for c,(name,_) in zip(vec,gens) if c} for vec in null],
    }

out = {'arithmetic': 'fractions.Fraction', 'source': 'Evans--Gannon Xi recurrence via vvmf_exact/obstruction'}
d3=V.names.index('d3')
for gamma in [0,1]:
    P=[0,d3] if gamma else [0]
    A=[(0,-1),(0,0)]+[(z,-1 if z in P else 0) for z in range(1,V.d)]
    data=describe(2,P,A,1)
    out[f'weight2_gamma{gamma}']=data
    print(f'gamma={gamma}: constraints={data["vanishing_matrix_shape"]}, rank={data["vanishing_matrix_rank"]}; relations={data["relations_rank"]}, image={data["evaluation_image_rank"]}',flush=True)

# An explicit dual form pairing with the sole remaining possible coefficient.
gens=O.Y_generators(4,1)
Hs=[V.H_from_Y(y) for _,y in gens]
positions=[]
for z in range(V.d):
    lim=-V.kz(z)-1
    nmin=min([min(h[V.PI[z]]) for h in Hs if h[V.PI[z]]] + [lim+1])
    positions.extend((z,n) for n in range(nmin,lim+1))
matrix=[[h[V.PI[z]].get(n,Fraction(0)) for h in Hs] for z,n in positions]
null=V.nullspace(matrix,len(gens))
candidates=[]
for vec in null:
    pairing=sum(c*h[V.PI[d3]].get(-V.kz(d3),Fraction(0)) for c,h in zip(vec,Hs))
    if pairing:
        vec=[c/pairing for c in vec]
        candidates.append(vec)
assert candidates
vec=min(candidates,key=lambda v:sum(bool(c) for c in v))
assert all(sum(c*x for c,x in zip(vec,row)) == 0 for row in matrix)
pairing=sum(c*h[V.PI[d3]].get(-V.kz(d3),Fraction(0)) for c,h in zip(vec,Hs))
assert pairing == 1
out['injectivity_gamma1']={
    'vanishing_matrix_shape':[len(matrix),len(gens)],
    'vanishing_matrix_rank':len(gens)-len(null),
    'pairing':str(pairing),
    'G_coefficients':{name:str(c) for c,(name,_) in zip(vec,gens) if c},
    'killed_positions':[[V.names[z],n] for z,n in positions],
}
print('injectivity witness:',out['injectivity_gamma1']['G_coefficients'],flush=True)

A=[(0,-1),(0,0)]+[(z,0) for z in range(1,V.d)]+[(0,1)]
for k in [2,4,6]:
    data=describe(k,[0],A,2)
    assert any(row[-1] for row in data['relations_rref'] if Fraction(row[-1]))
    out[f'augmented_weight{k}']=data
    print(f'augmented k={k}: constraints={data["vanishing_matrix_shape"]}, rank={data["vanishing_matrix_rank"]}; relations={data["relations_rank"]}, image={data["evaluation_image_rank"]}',flush=True)
args.output.write_text(json.dumps(out,indent=2)+'\n')
print('All exact checks passed; wrote',args.output,flush=True)
