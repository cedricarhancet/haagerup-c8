#!/bin/bash
# Ordre d'exécution des calculs de la version 03 (dossier Calculs/classification).
set -e
python3 weight2_obstruction_exact.py 1 1 > weight2_obstruction_exact_gamma1_output.txt
python3 weight2_obstruction_exact.py 0 1 > weight2_obstruction_exact_gamma0_output.txt
python3 dump_relations.py > dump_relations_output.txt
python3 gamma1_exclusion_exact.py > gamma1_exclusion_exact_output.txt
python3 weight4_formula_check.py > weight4_formula_check_output.txt
python3 weight6_formula_check.py > weight6_formula_check_output.txt
python3 e6_tools.py > e6_tools_output.txt
python3 weight4_solve.py > weight4_solve_output.txt
python3 riesen_quartic_check.py > riesen_quartic_check_output.txt
python3 centre_gamma0.py > centre_gamma0_output.txt
python3 sl2_e6_superset.py > sl2_e6_superset_output.txt
python3 e6_centralizers.py > e6_centralizers_output.txt
python3 e6_orbit_labels.py > e6_orbit_labels_output.txt
python3 sl2_ideal_feasibility.py > sl2_ideal_feasibility_output.txt
python3 sl2_ideal_feasibility.py --level > sl2_ideal_feasibility_level_output.txt
python3 sl2_pair_feasibility.py > sl2_pair_feasibility_output.txt
python3 sl2_pair_exact.py > sl2_pair_exact_output.txt
python3 weight6_single.py > weight6_single_output.txt
python3 weight6_single.py --pass2 > weight6_single_pass2_output.txt
python3 weight6_mixed.py > weight6_mixed_output.txt
python3 sextic_joint.py > sextic_joint_output.txt
python3 weightn_single_check.py 8 10 > weightn_single_check_output.txt
python3 degree8_feasibility.py > degree8_feasibility_output.txt
python3 joint_uniqueness.py --d3 > joint_uniqueness_output.txt
python3 known_sectors_table.py > known_sectors_table_output.txt
python3 u2_moments.py > u2_moments_output.txt
python3 u2_embedding_test.py > u2_embedding_test_output.txt
python3 u2_parameter_count.py > u2_parameter_count_output.txt
python3 cartan_conjugate.py > cartan_conjugate_output.txt
python3 lattice_branching.py > lattice_branching_output.txt
python3 u3_moments.py > u3_moments_output.txt
python3 u3_embedding_test.py > u3_embedding_test_output.txt
echo "tous les calculs ont abouti"
