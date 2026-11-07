#!/usr/bin/env bash
set -eu
ROOT=$(cd -- "$(dirname -- "$0")" && pwd)
DIR=$1
NAME=$2
shift 2
OD="$ROOT/replay"
mkdir -p "$OD"
cd -- "$DIR"
PYTHON=${PYTHON:-python3}
if "$PYTHON" -I -X utf8 -B -u -c 'import os,runpy,sys; sys.argv=sys.argv[1:]; sys.path.insert(0,os.getcwd()); runpy.run_path(sys.argv[0],run_name="__main__")' "$@" > "$OD/$NAME.out" 2> "$OD/$NAME.err"; then
    printf '%s\t%s\n' "$NAME" ok >> "$OD/TIMINGS.tsv"
else
    RC=$?
    printf '%s\texit=%s\n' "$NAME" "$RC" >> "$OD/TIMINGS.tsv"
    exit "$RC"
fi
