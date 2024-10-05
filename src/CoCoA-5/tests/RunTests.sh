#! /bin/bash

# Script for running the CoCoA-5 tests (copied from src/tests/RunTests.sh)
# This script is normally called by make as part of the target "check";
# it is not really intended to be called manually.

# For each file called test-XYZ (where XYZ can be changed) there may
# be some other related files:
#   test-XYZ.in   -- sample input for the test program
#   test-XYZ.out  -- expected output (on cout) from the test program
#   test-XYZ.err  -- expected output (on cerr) from the test program

# When test-XYZ is run, its output is placed into two files:
#   test-XYZ.found -- actual output (on cout) from the test program
#   test-XYZ.cerr  -- actual output (on cerr) from the test program
# A correct test should always exit with code 0; if the executable exits
# with a non-zero code, this is regarded as a failure.  If the exit code
# is 0 then the actual output is compared with the expected output.
# If they match then test-XYZ.found and test-XYZ.cerr are deleted,
# otherwise a failure message is printed, and the files are not deleted.

# The output files are compared using "diff -w" to work around a
# gratuitous incompatibility introduced by Microsoft.

if [ $# -eq 0 ]
then
  echo "$0: ERROR: no tests specified."
  echo "$0: usage: $0 <list of CoCoA test scripts>"
  exit 1
fi

if [ \! -d ../../../src/CoCoA-5/tests ]
then
    echo "$0: curr dir must be src/CoCoA-5/tests/"
    exit 1
fi

source ../../../configuration/shell-fns.sh
#source ../../../configuration/autoconf.mk

unset progs
for arg in "$@"; do
  prog=$(basename "$arg" .cocoa5)
  HAVE_EXTLIB=HAVE_"${prog:6:100}"
  if [ "ExtLib" != "${prog:0:6}" ] || [ "${!HAVE_EXTLIB}" = yes ]; then
    progs=("${progs[@]}" "$arg")
  fi
done
NUM_TESTS=${#progs[@]}-CoCoA-5

# Before running the tests check that empty input does not trigger an error (e.g. package error)
../CoCoAInterpreter --no-readline --no-preamble --packageDir ../packages /dev/null > NULL.found 2> NULL.cerr
if [ $? -ne 0 ] || [ -s NULL.found ] || [ -s NULL.cerr ]
then
    if [ -s NULL.found ] || [ -s NULL.cerr ]
    then
	echobox "EMPTY TEST FAILED (unexpected output)"
    else
	echobox "EMPTY TEST FAILED (non-zero exit status)"
    fi
    cat NULL.found  NULL.cerr
    /bin/rm  NULL.found  NULL.cerr
    echo "----------------------------------------------------"
    echo
    sleep 2
    exit 1
fi
/bin/rm  NULL.found  NULL.cerr

echo
echounderline "Running the CoCoA-5 tests (${#progs[@]} tests altogether)"

# Keep track of which tests failed, to print a summary at the end.
failures=""

# This loop iterates through the names of tests to run.
COUNTER=0
for p in "${progs[@]}"; do
  COUNTER=$(( 1 + $COUNTER ))
  prog=$(basename "$p" .cocoa5)
  if [ $? -ne 0 ] || [ ! -f "$p" ]
  then
      echo "!!!!! Bad test source file \`$p' !!!!!"
      exit 1
  fi
  /bin/rm -f  "$prog.found"  "$prog.cerr"
  ../CoCoAInterpreter --no-readline --no-preamble --packageDir ../packages "$prog.cocoa5"  > "$prog.found"  2> "$prog.cerr"
  if [ $? -ne 0 ]
  then
    echo "[$COUNTER/$NUM_TESTS] *****  $prog FAILED  ***** (non-zero exit status)"
    failures="$failures  $prog"
  else
    if [ -f "$prog.out" ]
    then
      diff -w  "$prog.found"  "$prog.out" > /dev/null
    else
      diff -w  "$prog.found"  /dev/null > /dev/null
    fi
    if [ $? -ne 0 ]
    then
      echo "[$COUNTER/$NUM_TESTS] *****  $prog FAILED  ***** (wrong output)"
      failures="$failures  $prog"
    else
      if [ -f "$prog.err" ]
      then
        diff -w  "$prog.cerr"  "$prog.err" > /dev/null
      else
        diff -w  "$prog.cerr"  /dev/null > /dev/null
      fi
      if [ $? -ne 0 ]
      then
        echo "[$COUNTER/$NUM_TESTS] *****  $prog FAILED  ***** (wrong output on cerr/clog)"
        failures="$failures  $prog"
      else
        /bin/rm  "$prog.found"  "$prog.cerr"
        echo "[$COUNTER/$NUM_TESTS] $prog.cocoa5 ..... OK"
      fi
    fi
  fi
done
if [ -z "$failures" ]
then
  echo "==================================="
  echo "Good news: all CoCoA-5 tests passed"
  echo "==================================="
  echo
  exit 0
fi


echo "**********************"
echo "*****  Bad news  *****"
echo "**********************"
echo "*****  The following CoCoA-5 tests failed, please tell us about it."
echo "*****  $failures"
exit 1
