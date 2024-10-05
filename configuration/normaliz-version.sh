#! /bin/bash

# Script to print out libnormaliz version number.
# Expects env variables CXX and CXXFLAGS to be set.


SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

if [ -z "$CXX" ]
then
  echo "ERROR: CXX environment variable is not defined   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

if [ $# -ne 2 ]
then
    echo "ERROR: expected 2 args: NORMALIZ_INC_DIR & NORMALIZ_LIB (full paths)    $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

NMZ_INC_DIR="$1"
if [ \! -d "$NMZ_INC_DIR" -o \! -f "$NMZ_INC_DIR/libnormaliz/libnormaliz.h" ]
then
    echo "ERROR: unable to find/read Normaliz headers in $NMZ_INC_DIR   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

NMZ_LIB="$2"
if [ \! -f "$NMZ_LIB" -o \! -r "$NMZ_LIB" ]
then
    echo "ERROR: specified Normaliz lib ($NMZ_LIB) not a readable file   $SCRIPT_NAME"  > /dev/stderr
    exit 1
fi

NMZ_LIB_DIR=`dirname "$NMZ_LIB"`
NMZ_LIB_BASE=`basename "$NMZ_LIB"`

# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir normaliz-version`

pushd "$TMP_DIR"  > /dev/null
/bin/cat > NormalizVersion.C  <<EOF
#include "libnormaliz/libnormaliz.h"

#include <iostream>

int main()
{
  std::cout << libnormaliz::getVersion() << std::endl;
  return 0;
}
EOF


echo "$CXX $CXXFLAGS  -I \"$NMZ_INC_DIR\"  NormalizVersion.C  \"$NMZ_LIB\"   -o NormalizVersion"  > LogFile 
$CXX $CXXFLAGS  -I "$NMZ_INC_DIR"   NormalizVersion.C  "$NMZ_LIB"  -o NormalizVersion         >> LogFile  2>&1
if [ $? -ne 0 ]
then
  echo "ERROR: Compilation of test program failed   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"                  > /dev/stderr
  exit 3  # do not clean TMP_DIR, for possible debugging
fi
NMZ_VER=`./NormalizVersion 2> LogFile`
if [ $? -ne 0 ]
then
  echo "ERROR: test program crashed   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"    > /dev/stderr
  exit 1
fi

# Clean up TMP_DIR
popd  > /dev/null
/bin/rm -rf "$TMP_DIR"
echo $NMZ_VER
exit 0
