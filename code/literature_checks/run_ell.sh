#!/bin/bash
for u in 1000 100; do
 for pair in "24 200" "32 150" "48 80" "64 40" "96 16"; do
  set -- $pair; python3 ell_fss.py $u $1 $2 $((8000+$1+${u%.*})); done; done
