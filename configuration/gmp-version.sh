#!/bin/bash

SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

# This script gets the version number of GMP.
# Assumes that env variables CXX and EXTLIB_DIR_FULL are set, and
# that a link to GMP header file is in $EXTLIB_DIR_FULL/include if
# there is no system default GMP.

# If an error occurs, exit code is non-zero (& mesg is printed on stderr).
# Otherwise exit code is zero, and output is version number (e.g. 6.0.1)

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

if [ $# -ne 0 ]
then
  echo "ERROR: expected no args   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX must be set to a C++ compiler compatible with GMP   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

if [ -z "$EXTLIB_DIR_FULL" ]
then
    echo "ERROR: environment variable EXTLIB_DIR_FULL not set.   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi


# The following is a cryptic if...then block
is_absolute "$EXTLIB_DIR_FULL" ||
(
  echo "ERROR: environment variable EXTLIB_DIR_FULL is not absolute: \"$EXTLIB_DIR_FULL\"   $SCRIPT_NAME"  > /dev/stderr
  exit 1
)

if [ \! -d "$EXTLIB_DIR_FULL" -o \! -d "$EXTLIB_DIR_FULL/include" -o \! -d "$EXTLIB_DIR_FULL/lib" ]
then
  echo "ERROR: environment variable EXTLIB_DIR_FULL is implausible: \"$EXTLIB_DIR_FULL\"   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


# Get version number from the header file; we (ab)use the compiler.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir gmp-version`

pushd "$TMP_DIR"  > /dev/null

/bin/cat > test-gmp-version.C <<EOF
#include "gmp.h"
#include <iostream>

int main()
{
  std::cout << __GNU_MP_VERSION << "." << __GNU_MP_VERSION_MINOR << "." << __GNU_MP_VERSION_PATCHLEVEL << std::endl;
}
EOF

# Use c++ compiler specified in CXX; no need to specify libgmp as all info is in header file!!
echo "$CXX -I"$EXTLIB_DIR_FULL/include"  test-gmp-version.C  -o test-gmp-version"  > LogFile
$CXX -I"$EXTLIB_DIR_FULL/include"  test-gmp-version.C  -o test-gmp-version  >> LogFile 2>&1

# Check whether compilation failed; if so, complain.
if [ $? -ne 0 ]
then
  echo "ERROR: unable to determine version of GMP library   $SCRIPT_NAME"   > /dev/stderr
  echo "ERROR: (compilation failed in gmp-version.sh)       $SCRIPT_NAME"   > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile     $SCRIPT_NAME"                         > /dev/stderr
  exit 1
fi

# Compilation succeeded, so run $PROG which will print out the version.
GMP_LIB_VERSION=`./test-gmp-version`

# Clean up TMP_DIR
popd  > /dev/null
/bin/rm -rf "$TMP_DIR"

# If we get here, all tests have passed, so print version number and exit with code 0.
echo $GMP_LIB_VERSION
