
# Script for running the CoCoALib tests (after they've been compiled).
# This script is normally called by make as part of the target "check";
# it is not really intended to be called manually.

# For each executable called test-XYZ (where XYZ can be changed) there may
# be some other related files:
#   test-XYZ.in   -- sample input for the test program
#   test-XYZ.out  -- expected output (on cout) from the test program
#   test-XYZ.err  -- expected output (on cerr) from the test program

# When test-XYZ is run, its output is placed into two files:
#   test-XYZ.cout-found  -- actual output (on cout) from the test program
#   test-XYZ.cerr-found  -- actual output (on cerr) from the test program
# A correct test should always exit with code 0; if the executable exits
# with a non-zero code, this is regarded as a failure.  If the exit code
# is 0 then the actual output is compared with the expected output.
# If they match then test-XYZ.cout-found and test-XYZ.cerr-found are deleted,
# otherwise a failure message is printed, and the files are not deleted.

# The output files are compared using "diff -w" to work around a
# gratuitous incompatibility introduced by Microsoft.

if [ $# -eq 0 ]
then
  echo "$0: ERROR: no tests specified."
  echo "$0: usage: $0 <list of executables>"
  exit 1
fi

source ../../configuration/shell-fns.sh

NUM_TESTS=$#-CoCoALib

echo
echounderline "Running the CoCoALib tests ($NUM_TESTS tests altogether)"

# Keep track of which tests failed, to print a summary at the end.
failures=""

# This loop iterates through the names supplied as arguments to the script.
COUNTER=0
while [ $# -ne 0 ]
do
  COUNTER=`expr 1 + $COUNTER`
  prog="$1"; shift
  /bin/rm -f  $prog.cout-found  $prog.cerr-found
  if [ \! -x $prog ]
  then
    echo "[$COUNTER/$NUM_TESTS] *****  $prog MISSING or NOT EXECUTABLE  *****"
  fi
  if [ -f $prog.in ]
  then
    ./$prog < $prog.in > $prog.cout-found 2> $prog.cerr-found
  else
    ./$prog > $prog.cout-found 2> $prog.cerr-found
  fi
  if [ $? -ne 0 ]
  then
    echo "[$COUNTER/$NUM_TESTS] *****  $prog FAILED  ***** (non-zero exit status)"
    failures="$failures $prog"
  else
    if [ -f $prog.out ]
    then
      diff -w  $prog.cout-found  $prog.out  > /dev/null
    else
      diff -w  $prog.cout-found  /dev/null  > /dev/null
    fi
    if [ $? -ne 0 ]
    then
      echo "[$COUNTER/$NUM_TESTS] *****  $prog FAILED  ***** (wrong output)"
      failures="$failures $prog"
    else
      if [ -f $prog.err ]
      then
        diff -w  $prog.cerr-found  $prog.err > /dev/null
      else
        diff -w  $prog.cerr-found  /dev/null > /dev/null
      fi
    if [ $? -ne 0 ]
      then
        echo "[$COUNTER/$NUM_TESTS] *****  $prog FAILED  ***** (wrong output on cerr/clog)"
        failures="$failures $prog"
      else
        /bin/rm  $prog.cout-found  $prog.cerr-found
        echo "[$COUNTER/$NUM_TESTS] $prog..... OK"
      fi
    fi
  fi
done
if [ -z "$failures" ]
then
  echo "===================================="
  echo "Good news: all CoCoALib tests passed"
  echo "===================================="
  echo
  exit 0
fi
#if [ "$failures" != " test-RingTwinFloat3 test-RingTwinFloat6" ]
#then
  echo "**********************"
  echo "*****  Bad news  *****"
  echo "**********************"
  echo "*****  The following CoCoALib tests failed, please tell us about it."
  echo "*****  $failures"
  exit 1
#fi
