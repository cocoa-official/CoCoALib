
# Script for running the CoCoALib tests (after they've been compiled).
# This script is normally called by make as part of the target "check";
# it is not really intended to be called manually.

# For each executable called test-XYZ (where XYZ can be changed) there may
# be some other related files:
#   test-XYZ.in   -- sample input for the test program
#   test-XYZ.out  -- expected output (on cout) from the test program
#   test-XYZ.err  -- expected output (on cerr) from the test program

# When test-XYZ is run, its output is placed into two files:
#   test-XYZ.cout -- actual output (on cout) from the test program
#   test-XYZ.cerr -- actual output (on cerr) from the test program
# A correct test should always exit with code 0; if the executable exits
# with a non-zero code, this is regarded as a failure.  If the exit code
# is 0 then the actual output is compared with the expected output.
# If they match then test-XYZ.cout and test-XYZ.cerr are deleted,
# otherwise a failure message is printed, and the files are not deleted.

# The output files are compared using "diff -wq" to work around a
# gratuitous incompatibility introduced by Microsoft.

if [ $# -eq 0 ]
then
  echo "$0: ERROR: no tests specified."
  echo "$0: usage: $0 <list of executables>"
  exit 1
fi

source ../configuration/shell-fns.sh

echo
echounderline "Running valgrind on the CoCoALib examples ($# programs altogether)"

# Keep track of which tests failed, to print a summary at the end.
failures=""

# This loop iterates through the names supplied as arguments to the script.
while [ $# -ne 0 ]
do
  prog="$1"; shift
  echo " --$prog "
  if [ -f "$prog.in" ]
  then
    valgrind "./$prog" < "$prog.in"  2>&1 | grep "definitely lost\|no leaks are possible\|Invalid write"
  else
    valgrind "./$prog"               2>&1 | grep "definitely lost\|no leaks are possible\|Invalid write"
  fi
  /bin/rm -f "$prog.cout"
done
