#!/bin/bash
s=$1
for pair in "24 400" "32 300" "48 150" "64 100" "96 40" ; do
  set -- $pair; python3 span_open.py $s $1 $2 $((7000+$1)) ; done
