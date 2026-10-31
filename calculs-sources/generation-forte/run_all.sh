#!/bin/bash
# Ordre d'exécution (voir README-generation-forte.txt). Deux calculs au plus en parallèle.
set -e
python3 certify_w3.py > certify_w3.out
python3 count_generators.py > count_generators.out
python3 certify_w4b.py > certify_w4.out            # pivots (piv_w4/) puis tests (tests_w4/) ; peut être interrompu après les pivots
python3 certify_w4_part.py '(3, 3);(2, 6);(1, 7);(1, 5);(2, 2);(1, 3)' tests_w4_part > certify_w4_part.out
python3 solve_w4.py > solve_w4.out
python3 phase_d_w4.py > phase_d_w4.out
python3 commutant_dims.py > commutant_dims.out
