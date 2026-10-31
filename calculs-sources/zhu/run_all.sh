#!/bin/bash
# Zhu computation for case II (v07). Requires Calculs/generation-forte (exact engine, exact_S3_basis.pkl, piv_w4/).
# Memory < 1 GB except z1 (~0.6 GB). Times on 2 cores: z1_J 12 s, z0 11 s, z1 25 min, z2 1 min, z3 1 min, z4 4 min.
cd "$(dirname "$0")"
(cd ../generation-forte && python3 ../zhu/z1_J_in_S.py) && python3 z0_highest.py && python3 z1_products.py && python3 z2_relations_modp.py && \
python3 z3_corner_modp.py && python3 z4_upper_algebra_modp.py
# second prime: ./run_p2.sh
