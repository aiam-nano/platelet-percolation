#!/bin/bash
# usage: run_shapes.sh shape  (seeds 100*index+L)
s=$1; k=0
for pair in "16 600" "24 400" "32 300" "48 150" "64 100" "96 40" "128 16"; do
  set -- $pair; python3 shape_fss.py $s $1 $2 $((5000+$1)) ; done
