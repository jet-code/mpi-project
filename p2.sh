#!/bin/bash

# MATH 6378 Project #2
# Alexander Hebert
# 

script=$( basename $0 )
if [ $# -ne 6 ]; then
  echo "Usage: $script xl xr ul ur u(x) f(x) nits" >&2
  exit 1
fi

HOMEDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
SRCFILE=$HOMEDIR/p2_v2-time.c
# EXEFILE=/tmp/project2-$USER-$$
EXEFILE=$HOMEDIR/p2

gcc -DXL=$1 -DXR=$2 -DUL=$3 -DUR=$4 -DUOFX="$5" -DFOFX="$6" $SRCFILE -o $EXEFILE -lm6378 -lm


# rm -f $EXEFILE












