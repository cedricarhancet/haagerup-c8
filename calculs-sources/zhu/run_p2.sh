#!/bin/bash
# second prime: whole modular chain (products, relations, corner, upper algebra)
export HAAG_PSTART=4194304 HAAG_PDIR=/home/claude/work/Calculs/zhu/p4194277
cd /home/claude/work/Calculs/zhu
python3 z1_products.py > p4194277/z1.out 2>&1 && python3 z2_relations_modp.py > p4194277/z2.out 2>&1 && \
python3 z3_corner_modp.py > p4194277/z3.out 2>&1 && python3 z4_upper_algebra_modp.py > p4194277/z4.out 2>&1
echo finished > p4194277/FINISHED
