#!/bin/bash

SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

# Auxiliary script for CoCoALib configuration process.
# Script expects the env variables CXX and CXXFLAGS to be set.

# Script to see whether the -std=c++14 compiler flag is needed/recognised.
# Exit with non-0 code if we did not find a way to compile C++14 code.
# Exit with code 0 if we found a way to compile C++14 code; printed
# value is flag to give compiler to get C++14 compilation
# (printed value may be empty string or "-std=c++14")


if [ $# -ne 0 ]
then
  echo "ERROR: expected no args.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

# Check environment variable CXX
if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX not set.   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir cxx14`

pushd "$TMP_DIR"  > /dev/null


/bin/cat > language-version.C <<EOF
// Inspired by https://stackoverflow.com/questions/2324658/how-to-determine-the-version-of-the-c-standard-used-by-the-compiler

#include<iostream>

int main() {
    if (__cplusplus == 202002L) std::cout << "C++20\n";
    else if (__cplusplus == 201703L) std::cout << "C++17\n";
    else if (__cplusplus == 201402L) std::cout << "C++14\n";
    else if (__cplusplus == 201103L) std::cout << "C++11\n";
    else if (__cplusplus == 199711L) std::cout << "C++98\n";
    else std::cout << "non-standard C++\n";
}
EOF

# First try with no compiler flag...
"$CXX"  language-version.C  -o language-version  >> LogFile  2>& 1 
if [ $? -ne 0 ]
then
    echo "ERROR: compilation unexpectedly failed; is $CXX a c++ compiler?   $SCRIPT_NAME" > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"                                       > /dev/stderr
    exit 1
fi
./language-version | tee language-version.out  >> LogFile
if [ $? -ne 0 ]
then
    echo "ERROR: language-version program crashed unexpectedly    $SCRIPT_NAME"           > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile    $SCRIPT_NAME"                                      > /dev/stderr
    exit 1
fi
CXXVER=`/bin/cat language-version.out`
if [ "$CXXVER" = "CXX14" ]
then
    popd  > /dev/null
    /bin/rm -rf "$TMP_DIR"
    exit 0; # exit without printing (no flag needed for C++14)
fi

# Compilation without flag is not C++14 standard; try with -std=c++14

CXX14="-std=c++14"
"$CXX"  $CXX14  language-version.C  -o language-version  >> LogFile  2>& 1 
if [ $? -ne 0 ]
then
    echo "ERROR: compilation with flag $CXX14 failed   $SCRIPT_NAME"                      > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"                                       > /dev/stderr
    exit 1
fi
./language-version | tee language-version.out  >> LogFile
if [ $? -ne 0 ]
then
    echo "ERROR: language-version program crashed unexpectedly   $SCRIPT_NAME"            > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile    $SCRIPT_NAME"                                      > /dev/stderr
    exit 1
fi
CXXVER=`/bin/cat language-version.out`
if [ "$CXXVER" = "C++14" ]
then
    popd  > /dev/null
    /bin/rm -rf "$TMP_DIR"
    echo "$CXX14"
    exit 0; # Success (flag for C++14 sent via echo to stdout)
fi

# With luck we never get here
echo "ERROR: failed to find flag for C++14 compilation   $SCRIPT_NAME"                    > /dev/stderr
echo "LOGFILE: $TMP_DIR/LogFile    $SCRIPT_NAME"                                          > /dev/stderr
exit 2
