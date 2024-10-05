#! /bin/bash

# Script to choose the best ULL2LL defn (in ULongLong2LongLong.H).
# Expects env variables CXX and CXXFLAGS to be set.

SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

# Try the three possible settings for the flag CoCoA_ULONGLONG2LONGLONG.
# Print out the first which works (and exit with code 0): prints 1, 2 or 3
# (this is the value to assign to preprocessor symbol CoCoA_ULONGLONG2LONGLONG).
# If none works, print an error message on stderr, and exit with non-zero code.
# Exit code 1 means some input was bad.
# Exit code 2 means no working defn of ULong2Long was found.
# Exit code 3 means a problem in compilation.

# taken from StackExchange 256434
is_absolute()
{
    case "$1" in
	///* | //) true;;
	//*) false;;
	/*) true;;
	*) false;;
    esac
}


if [ $# -ne 1 ]
then
  echo "ERROR: expected 1 arg (full path of CoCoALib include directory)   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

# The following is a cryptic if...then block
is_absolute "$1" ||
(
  echo "ERROR: arg is not absolute \"$1\"   $SCRIPT_NAME"   > /dev/stderr
  exit 1
)

if [ \! -d "$1" ]
then
    echo "ERROR: arg is not a directory \"$1\"   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi
    
COCOA_INC_DIR=`dirname "$1"`
BASE=`basename "$1"`
if [ "$BASE" \!= "CoCoA" -o \! -f "$COCOA_INC_DIR/CoCoA/ULongLong2LongLong.H" ]
then
  echo "ERROR: cannot find CoCoA header file \"$COCOA_INC_DIR/CoCoA/ULongLong2LongLong.H\"   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

if [ -z "$CXX" ]
then
  echo "ERROR: CXX environment variable is not defined   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir cpp-flags-ulonglong2longlong`

pushd "$TMP_DIR"  > /dev/null
/bin/cat > CheckULongLong2LongLong.C  <<EOF
#include "CoCoA/ULongLong2LongLong.H"
using CoCoA::ULongLong2LongLong;

#include <limits>
using std::numeric_limits;

// Check whether ULongLong2LongLong works at four negative values:
//   MinLongLong, MinLongLong+1, (MinLongLong+1)/2, and -1
// Exits with 0 iff the fn works correctly at all four points.
// A non-zero exit code indicates which check failed, but this
// information is not currently used anywhere.

// This program is compiled and run by cpp-flags-ulonglong2longlong.sh
// trying the various possible values for CoCoA_ULONGLONG2LONGLONG.

int main()
{
  volatile unsigned long long ull;
  long long ll = numeric_limits<long long>::min();
  ull = ll;
  if (ULongLong2LongLong(ull) != ll) return 1;

  ll = ll+1;
  ull = ll;
  if (ULongLong2LongLong(ull) != ll) return 2;

  ll = ll/2;
  ull = ll;
  if (ULongLong2LongLong(ull) != ll) return 3;

  ll = -1;
  ull = ll;
  if (ULongLong2LongLong(ull) != ll) return 4;

  return 0;
}
EOF


# Now try the various possible settings for CoCoA_ULONGLONG2LONGLONG
ULL2LL_DEFN=
for defn in 1 2 3
do
  CPPFLAG="-DCoCoA_ULONGLONG2LONGLONG=$defn"
  $CXX $CXXFLAGS $CPPFLAG -I"$COCOA_INC_DIR"  CheckULongLong2LongLong.C -o CheckULongLong2LongLong  > LogFile  2>&1
  if [ $? -ne 0 ]
  then
    echo "ERROR: Compilation of test program failed   $SCRIPT_NAME"  > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"  > /dev/stderr
    exit 3  # do not clean TMP_DIR, for possible debugging
  fi
  ./CheckULongLong2LongLong  >> LogFile  2>&1
  if [ $? = 0 ]
  then
    ULL2LL_DEFN="$defn"
    break
  fi
done

if [ -z "$ULL2LL_DEFN" ]
then
  echo "ERROR: Failed to determine a working defn of ULongLong2LongLong   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"  > /dev/stderr
  exit 2 # do not clean TMP_DIR, for possible debugging
fi

# Clean up TMP_DIR
popd  > /dev/null
/bin/rm -rf "$TMP_DIR"
echo $ULL2LL_DEFN
exit 0
