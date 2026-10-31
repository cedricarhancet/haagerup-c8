#!/bin/bash
# exact isotropy analysis, one process per candidate (Axis B3)
case=$1; shift; real=${case%b}
for c in "$@"; do timeout 1500 python3 exact_isotropy.py $real $c > exact_isotropy_${real}_$c.out 2>&1; done
echo DONE > exact_${case}.done
