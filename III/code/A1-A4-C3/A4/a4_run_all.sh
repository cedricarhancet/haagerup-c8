#!/bin/sh
set -eu
cd "$(dirname "$0")"
mkdir -p replay
PYTHON=${PYTHON:-python3}
"$PYTHON" -I -X utf8 a4_kfield.py > replay/a4_kfield.out
for name in a4_smoke a4_item1_lattice a4_item2_section5 a4_item3_section6 a4_item4_section7 a4_item5_prop81 a4_item6_riesen; do
    "$PYTHON" -I -X utf8 "$name.py" > "replay/$name.out"
    if grep -q '^FAUX' "replay/$name.out"; then
        echo "$name.py reported a failed check" >&2
        exit 1
    fi
done
