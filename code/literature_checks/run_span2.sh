#!/bin/bash
s=$1
for pair in "24 300" "32 200" "48 100" "64 60" ; do
  set -- $pair; python3 span_open.py $s $1 $2 $((7500+$1)) ; done
