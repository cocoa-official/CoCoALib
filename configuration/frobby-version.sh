#! /bin/bash

# Script to print out libfrobby version number.
# Expects env variables CXX and CXXFLAGS to be set.


SCRIPT_NAME=[[$(basename "$0")]]
SCRIPT_DIR=$(dirname "$0")

if [ -z "$CXX" ]
then
  echo "ERROR: CXX environment variable is not defined   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

if [ $# -ne 4 ]
then
    echo "ERROR: expected 4 args: FROBBY_INC_DIR & FROBBY_LIB & GMP_LIB & GMPXX_LIB (full paths)    $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

FROBBY_INC_DIR="$1"
if ! { [ -d "$FROBBY_INC_DIR" ] && [ -f "$FROBBY_INC_DIR/frobby.h" ]; }
then
    echo "ERROR: unable to find/read Frobby header in $FROBBY_INC_DIR   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

FROBBY_LIB="$2"
if ! { [ -f "$FROBBY_LIB" ] && [ -r "$FROBBY_LIB" ]; }
then
    echo "ERROR: specified Frobby lib ($FROBBY_LIB) not a readable file   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

GMP_LIB="$3"
if ! { [ -f "$GMP_LIB" ] && [ -r "$GMP_LIB" ]; }
then
    echo "ERROR: specified GMP lib ($GMP_LIB) not a readable file   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

GMPXX_LIB="$4"
if ! { [ -f "$GMPXX_LIB" ] && [ -r "$GMPXX_LIB" ]; }
then
    echo "ERROR: specified GMPXX lib ($GMPXX_LIB) not a readable file   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=$(mktempdir frobby-version)

pushd "$TMP_DIR"  > /dev/null   || ( echo "ERROR: push failed   $SCRIPT_NAME" > /dev/stderr; exit 2 )
/bin/cat > FrobbyVersion.C  <<EOF
#include "frobby.h"

#include <iostream>

int main()
{
  std::cout << frobby_version << std::endl;
  return 0;
}
EOF


echo "$CXX $CXXFLAGS  -I \"$FROBBY_INC_DIR\"  FrobbyVersion.C  \"$FROBBY_LIB\"  \"$GMP_LIB\"  \"$GMPXX_LIB\"   -o FrobbyVersion"  > LogFile
$CXX $CXXFLAGS  -I "$FROBBY_INC_DIR"   FrobbyVersion.C  "$FROBBY_LIB"  "$GMP_LIB"  "$GMPXX_LIB"  -o FrobbyVersion         >> LogFile  2>&1
if [ $? -ne 0 ]
then
  echo "ERROR: Compilation of test program failed   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"                  > /dev/stderr
  exit 3  # do not clean TMP_DIR, for possible debugging
fi
if ! FROBBY_VER=$(./FrobbyVersion 2> LogFile)
then
  echo "ERROR: test program crashed   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"    > /dev/stderr
  exit 1
fi

# Clean up TMP_DIR
popd  > /dev/null  || ( echo "ERROR: popd failed   $SCRIPT_NAME" > /dev/stderr; exit 2 )
/bin/rm -rf "$TMP_DIR"
echo "$FROBBY_VER"
exit 0
