#!/usr/bin/env bash
set -eu
ROOT=$(cd -- "$(dirname -- "$0")" && pwd)
# Controles supplementaires (non utilises par les preuves du par. 7.5, mais compares aux sorties d'origine quand elles existent)
R="$ROOT/run/R3"
RUN="$ROOT/run_one.sh"
bash "$RUN" "$R" r2_reduced_0_d2_d3_b r2_reduced.py 0 d2,d3,b
bash "$RUN" "$R" r2_reduced_0_d4_d1_d5_d6_c1_c2 r2_reduced.py 0 d4,d1,d5,d6,c1,c2
for CI in 1 2 3; do bash "$RUN" "$R" r2_reduced_${CI}_all r2_reduced.py $CI d1,d2,d3,d4,d5,d6,c1,c2,b; done
for CI in 1 2 3; do bash "$RUN" "$R" r3a_invariants_$CI r3a_invariants.py $CI; done
for CI in 1 2 3; do bash "$RUN" "$R" r3a_1S2_$CI r3a_1S2.py $CI; done
echo "extra chain done $(date -u +%FT%TZ)" >> "$ROOT/replay/TIMINGS.tsv"
