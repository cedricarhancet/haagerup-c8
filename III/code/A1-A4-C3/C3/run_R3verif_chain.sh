#!/usr/bin/env bash
set -eu
ROOT=$(cd -- "$(dirname -- "$0")" && pwd)
# Chaine R3-verif (seconde implementation, ecrite de zero), ordre du README.txt de R3-verif.
# Les scripts ecrivent dans le repertoire courant : on les lance depuis run/R3-verif.
R="$ROOT/run/R3-verif"
RUN="$ROOT/run_one.sh"
bash "$RUN" "$R" verif_test_voa          test_voa.py
bash "$RUN" "$R" verif_check_W1_char     check_W1_char.py
bash "$RUN" "$R" verif_solve_signs       solve_signs.py
bash "$RUN" "$R" verif_embed             embed.py
bash "$RUN" "$R" verif_g_space           g_space.py
bash "$RUN" "$R" verif_griess            griess.py
bash "$RUN" "$R" verif_setup_system      setup_system.py
bash "$RUN" "$R" verif_solve_groebner    solve_groebner.py
bash "$RUN" "$R" verif_verify_solutions  verify_solutions.py
bash "$RUN" "$R" verif_phi_chi           phi_chi.py
bash "$RUN" "$R" verif_s2_space_0_all    s2_space.py 0 all
bash "$RUN" "$R" verif_s2_check_0        s2_check.py 0
bash "$RUN" "$R" verif_twist_check       twist_check.py
echo "R3-verif chain done $(date -u +%FT%TZ)" >> "$ROOT/replay/TIMINGS.tsv"
