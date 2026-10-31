#!/bin/bash
# Identification of Jiao's algebra with S (v08). Run from this folder; requires ../generation-forte and ../zhu.
# Times: i2 < 1 min, i3 < 1 min, i4 < 1 min, i5 ~ 2 min, i6 ~ 3 min.
cd "$(dirname "$0")"
python3 i2_weyl.py && python3 i3_transport.py && python3 i4_torus.py && python3 i5_membership.py && python3 i6_generation.py
