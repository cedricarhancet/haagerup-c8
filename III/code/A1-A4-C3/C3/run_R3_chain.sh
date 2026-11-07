#!/usr/bin/env bash
set -eu
ROOT=$(cd -- "$(dirname -- "$0")" && pwd)
# Chaine R3 (premiere implementation), ordre du README-R3.txt. Les fichiers pkl/json regenerables ont ete supprimes
# de la copie avant le lancement ; chaque script les regenere.
R="$ROOT/run/R3"
RUN="$ROOT/run_one.sh"
bash "$RUN" "$R" test_e8e8        test_e8e8.py
bash "$RUN" "$R" sl2x4            sl2x4.py
bash "$RUN" "$R" weights_W        weights_W.py
bash "$RUN" "$R" t12b_optionB     t12b_optionB.py
bash "$RUN" "$R" r3a_G            r3a_G.py
bash "$RUN" "$R" r3a_Griess       r3a_Griess.py
bash "$RUN" "$R" r3a_groebner     r3a_groebner.py
bash "$RUN" "$R" r3a_omega_exact  r3a_omega_exact.py
bash "$RUN" "$R" r3b_chi          r3b_chi.py
for CI in 0 1 2 3; do bash "$RUN" "$R" r3a_S2_$CI r3a_S2.py $CI; done
bash "$RUN" "$R" r3a_S2_full_0    r3a_S2_full.py 0 full
bash "$RUN" "$R" r3b_verify_0     r3b_verify.py 0
bash "$RUN" "$R" r2_signature_0   r2_signature.py 0
bash "$RUN" "$R" w_positivity     w_positivity.py
bash "$RUN" "$R" r3a_invariants_0 r3a_invariants.py 0
bash "$RUN" "$R" r3a_1S2_0        r3a_1S2.py 0
for CI in 1 2 3; do bash "$RUN" "$R" r3a_S2_full_$CI r3a_S2_full.py $CI; done
for CI in 1 2 3; do bash "$RUN" "$R" r3b_verify_$CI r3b_verify.py $CI; done
for CI in 1 2 3; do bash "$RUN" "$R" r2_signature_$CI r2_signature.py $CI; done
echo "R3 chain done $(date -u +%FT%TZ)" >> "$ROOT/replay/TIMINGS.tsv"
