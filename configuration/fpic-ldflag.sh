#!/bin/bash

SCRIPT_NAME=[[$(basename "$0")]]
## UNUSED? SCRIPT_DIR=$(dirname "$0")

# Auxiliary script for CoCoALib configuration process.
# Expects env variable CXX to be set (to compiler's name).

# Script to see whether platform is MacOS X; if so, print out extra compiler flags.

# Exit code = 0 means all OK
# Exit code != 0 means an error occurred.


if [ $# -ne 0 ]
then
  echo "ERROR: expected no args   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

# Check environment variable CXX --- [2022-11-16] currently not used!
if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX not set.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


OS=$(uname -s)
PLATFORM=$(uname -m)
if [ "$OS" = "Darwin" -a "$PLATFORM" = "x86_64" ]
then
    echo "-Wl,-no_pie"  # note the underscore!
fi
