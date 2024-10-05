#!/bin/bash

SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

# Auxiliary script for CoCoALib configuration process.
# Script expects the env variables CXX and CXXFAGS to be set.

# Script to see whether the -fPIC flag produces annoying compiler warnings.
# If no warning is produced, the script prints "-fPIC"; otherwise it prints nothing.

if [ $# -ne 0 ]
then
  echo "ERROR: expected no args   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

# Check environment variable CXX
if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX not set.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


FPIC_FLAG=-fPIC

# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir fpic-flag`

pushd "$TMP_DIR"  > /dev/null

/bin/cat > test-fpic-flag.C <<EOF
int f(int x)
{
  return (x+1)*x+41;
}
EOF


COMPILER_MESG=`"$CXX" $FPIC_FLAG -c  -o test-fpic-flag.o  test-fpic-flag.C  2>& 1`
if [ $? -ne 0 ]
then
  echo "ERROR: test compilation failed   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

# Clean up TMP_DIR
popd  > /dev/null
/bin/rm -rf "$TMP_DIR"
if [ -z "$COMPILER_MESG" ]
then
  echo $FPIC_FLAG
fi
