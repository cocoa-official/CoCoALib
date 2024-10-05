#! /bin/bash

SCRIPT_NAME=[[`basename "$0"`]]

# Script to test whether libnormaliz wants OpenMP.
# If OpenMP is needed the script prints "OpenMP", otherwise nothing.

if [ $# != 1 ]
then
  echo "ERROR: expected 1 arg (full path of libnormaliz.a)   $SCRIPT_NAME" > /dev/stderr
  exit 1
fi

# Simply check for the presence of symbols beginning with "omp"
# Note redirection of error output from nm to /dev/null
nm "$1" 2>/dev/null | grep -q -F " omp" > /dev/null
if [ $? = 0 ]
then
  echo "OpenMP"
fi
