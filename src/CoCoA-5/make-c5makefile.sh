#! /bin/bash

# This script builds the CoCoA-5 GUI; it assumes that CoCoALib
# has already been built, and that the Qt libraries have been installed
# (and that the $QMAKE command is available).

# if [ $# -ne 5 ]
# then
#   echo "ERROR: $0 requires 5 args (full paths of BOOST_INC_DIR,  BOOST_LIB_DIR, libgmpxx, libgmp, gmp.h)"  > /dev/stderr
#   exit 1
# fi

##################################################################
# Check $QMAKE

if [ -z "$QMAKE" ]; then
  echo "ERROR: QMAKE not set -- unable to build Qt GUI" > /dev/stderr
  echo "ERROR: Please re-run configure script *without* the option \"--no-qt-gui\"" > /dev/stderr
  exit 2
fi

if [ ! -r "$QMAKE" -o ! -x "$QMAKE" ]; then
  echo "ERROR: $QMAKE not executable" > /dev/stderr
  exit 3
fi

##################################################################
# Set the variable ARCH if we are on a MacOS X system...

SYS_TYPE=`uname`

# On MacOSX we must explicitly state the "architecture".
# The block below should pick the right one (as dictated by libgmp).
if [ "$SYS_TYPE" = "Darwin" ]
then
  ARCH=`arch`  # either "i386" or "ppc"
  if [ "$ARCH" = "i386" ]
  then
    (cd ../../examples; make ex-empty) > /dev/null 2>&1
    arch -x86_64 ../../examples/ex-empty > /dev/null 2>&1
    if [ $? = 0 ]
    then
      ARCH=x86_64
    fi
  fi
  if [ "$ARCH" = "ppc" ]
  then
    (cd ../../examples; make ex-empty) > /dev/null 2>&1
    arch -ppc64 ../../examples/ex-empty > /dev/null 2>&1
    if [ $? = 0 ]
    then
      ARCH=ppc64
    fi
  fi
fi


##################################################################
# Generate C5.pro file

/bin/rm -f C5.pro

echo "## This file was generated automatically by $0"  > C5.pro
echo                                                  >> C5.pro
echo "INCLUDEPATH += ../../include"                   >> C5.pro
echo "INCLUDEPATH += ../../configuration/ExternalLibs/include"    >> C5.pro
echo "LIBS += ../../lib/libcocoa.a"                   >> C5.pro

# Process params: ASSUME all compilation flags start with minus, or
#                        if it's a dir then its an include dir
#                        o/w it's a library.
# WARNING: need to check how to guarantee paths with spaces:
# \"$1\" does not work because $1 might be two libraries
while [ $# -gt 0 ]
do
  if [ "XX-l" = "XX${1:0:2}" ] ## first 2 characters in $1 are "-l"
  then
    echo "LIBS += $1"                                     >> C5.pro
  elif [ "XX-L" = "XX${1:0:2}" ] ## first 2 characters in $1 is "-L"
  then
    echo "LIBS += $1"                                     >> C5.pro
  elif [ "XX-D" = "XX${1:0:2}" ] ## first 2 characters in $1 are "-D"
  then
    echo "DEFINES += ${1:2}  ##### ExternalLib #####"     >> C5.pro
  elif [ "XX-" = "XX${1:0:1}" ] ## first character in $1 is "-"
  then
    echo "QMAKE_CXXFLAGS += $1"  >> C5.pro
    if [ "XX$1" = "XX-fopenmp" ]
    then
      echo "QMAKE_LFLAGS += -fopenmp  ### for Normaliz"   >> C5.pro
    fi
  elif [ -d "$1" ] ## $1 is a directory
  then
    echo "INCLUDEPATH += $1"                              >> C5.pro
  else
    echo "LIBS += $1"                                     >> C5.pro
  fi
  shift
done


if [ "$SYS_TYPE" = "Darwin" ]
then
  echo "CONFIG += $ARCH    #MacOSX architecture"      >> C5.pro
fi
cat C5.pro.in                                         >> C5.pro


if [ "$SYS_TYPE" = "Darwin" ]
then
    # 2022-11-18 Don't really know what I'm doing here... I just hope it works
    if [ "$CXX" = "g++" ]
    then
	DARWIN_OPTS="-spec macx-g++"
    else
	DARWIN_OPTS="-spec macx-clang"
    fi
fi
"$QMAKE"  $DARWIN_OPTS  C5.pro  -o C5Makefile
if [ "$?" -ne 0 ]
then
  echo "ERROR: $0 failed because $QMAKE failed for C5.pro" > /dev/stderr
  exit 6
fi
