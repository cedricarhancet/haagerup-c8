#!/bin/bash
# Ordre d'exécution des calculs du cas (II) (version 04). Environ 30 minutes au total sur deux coeurs.
set -e
python3 test_lattice_voa.py > test_lattice_voa_output.txt
python3 cases.py > cases_output.txt
python3 check_invariants.py > check_invariants_output.txt
python3 commutant_griess.py > commutant_griess_output.txt
python3 griess_subalgebras.py > griess_subalgebras_output.txt
python3 primaries.py > primaries_output.txt
python3 isotropy_test.py > isotropy_test_output.txt
./run_exact.sh II 0 1 2 3 6 7 8 11 12 13 14
./run_exact.sh I 1 2 4
for a in "II 9" "II 10" "I 5" "I 6"; do python3 exact_isotropy_s13.py $a > exact_isotropy_s13_${a// /_}.out; done
python3 isotropic_lines.py II 9 > isotropic_lines_II_9.out
python3 isotropic_lines.py I 5 > isotropic_lines_I_5.out
python3 closure_explore.py II 9 > closure_explore_II_9.out
python3 closure_certify.py II 9 > closure_certify_II_9.out
python3 closure_full.py II 9 > closure_full_II_9.out
python3 closure_modp.py II 9 0 > closure_modp_II_9_0.out
python3 closure_modp.py II 9 -1 > closure_modp_II_9_rejete.out
C3U1=1 C3NAMES='e;x;(0, 4);(1, 1);(1, 3);(1, 5);(1, 7);(2, 2);(2, 6);(3, 3)' python3 closure_modp.py II 9 0 w3 > closure_modp_II_9_w3_C3all.out
PSTART=$((2**27-1000000)) NOC3=1 python3 closure_modp.py II 9 0 w3 > closure_modp_II_9_w3_p2.out
python3 weight3_explore.py II 9 > weight3_explore_II_9.out
python3 closure_explore2.py I 5 > closure_explore2_I_5.out
echo "tous les calculs ont abouti"
