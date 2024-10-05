#!/bin/bash

# This script looks for the necessary BOOST sub-libraries.
# If found, it prints out two SHELL assignments (to be processed
# using the eval commmand) to BOOST_LIB_DIR & BOOST_LDLIBS.
# It also creates symbolic links in $EXTLIB_5_DIR_FULL [var must be set]
# If not, it prints out a warning and returns with exit code 1.

SCRIPT_NAME=[[$(basename "$0")]]

# 2022-09-14: previously needed "system thread filesystem"
### 2022-11-12: my SmallBlue computer seems to need "system" to build the GUI
# BOOST 1.67 and earlier need both filesystem and system
SUBLIBS="filesystem  system"   # ORDER IS IMPORTANT!  system needed for older BOOST (before 1.74)

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

##################################################################
# Check we have 1 arg, and that it is sane.

if [ "$#" \!= "1" ]
then
    echo "ERROR: expected 1 arg (full path of BOOST_HDR_DIR)   $SCRIPT_NAME"   > /dev/stderr
    exit 1
fi

if [ \! -d "$1" -o \! -d "$1/boost" ]
then
  echo "ERROR: arg must be a directory containing subdir \"boost\"   $SCRIPT_NAME"   > /dev/stderr
  exit 2
fi

BOOST_HDR_DIR="$1"


# Check environment variable EXTLIB_5_DIR_FULL: expect full path to dir with subdirs .../lib & .../include

if [ -z "$EXTLIB_5_DIR_FULL" ]
then
    echo "ERROR: environment variable EXTLIB_5_DIR_FULL not set.   $SCRIPT_NAME"   > /dev/stderr
    exit 1
fi

# The following is a cryptic if...then block
is_absolute "$EXTLIB_5_DIR_FULL" ||
(
  echo "ERROR: environment variable EXTLIB_5_DIR_FULL is not absolute: \"$EXTLIB_5_DIR_FULL\"   $SCRIPT_NAME"   > /dev/stderr
  exit 1
)

if [ \! -d "$EXTLIB_5_DIR_FULL" -o \! -d "$EXTLIB_5_DIR_FULL/include" -o \! -d "$EXTLIB_5_DIR_FULL/lib" ]
then
    echo "ERROR: environment variable EXTLIB_5_DIR_FULL is implausible: \"$EXTLIB_5_DIR_FULL\"   $SCRIPT_NAME"   > /dev/stderr
    echo "ERROR: (expected subdirs .../include/ and .../lib/ to exist)"  > /dev/stderr
  exit 1
fi


# 2022-11-16: This block below was needed to be able to call the script extn.sh -- maybe no longer needed?
# #############################################################################
# # Find directory of this script (taken from link below)
# # http://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself

# pushd . > /dev/null
# SCRIPT_PATH="${BASH_SOURCE[0]}"
# while [ -L "$SCRIPT_PATH" ]
# do
#   cd "$(dirname "$SCRIPT_PATH")"
#   SCRIPT_PATH="$(readlink "$(basename "$SCRIPT_PATH")")"
# done
# cd "$(dirname "$SCRIPT_PATH")" > /dev/null
# SCRIPT_DIR="$(pwd)"
# popd  > /dev/null


# 2024-03-20 MinGW still uses the -mt suffix
has_libboost_mt ()
{
    MISSING_SUBLIBS=
    for sublib in $SUBLIBS
    do
	/bin/ls  "${1}"/libboost_$sublib-mt.*  >/dev/null  2>&1
	if [ $? -ne 0 ]
	then
	    MISSING_SUBLIBS="libboost_$sublib-mt  $MISSING_SUBLIBS"
	fi
    done
    test -z "$MISSING_SUBLIBS"
}

has_libboost ()
{
    MISSING_SUBLIBS=
    for sublib in $SUBLIBS
    do
	/bin/ls  "${1}"/libboost_$sublib.*  >/dev/null  2>&1
	if [ $? -ne 0 ]
	then
	    MISSING_SUBLIBS="libboost_$sublib  $MISSING_SUBLIBS"
	fi
    done
    test -z "$MISSING_SUBLIBS"
}


ARCH=$(uname -m)      #so far seen:  x86_64, aarch64
# CHEAP HACK: prefer lib/$ARCH-linux-gnu over lib64 over lib
# Remember: GMP probably chose 64-bits
BOOST_LIB_DIR_DIR="$(dirname "$BOOST_HDR_DIR")"
for subdir in  "lib/$ARCH-linux-gnu"  lib64  lib  lib/i386-linux-gnu
do
    DIR="$BOOST_LIB_DIR_DIR/$subdir"
    if [ -d "$DIR" ]
    then
	## [2022-11-16] BOOST 1.80 seems no longer to use -mt suffix (but their doc is out of date)
	# Next if stmt needed for MinGW
        if has_libboost_mt "$DIR"
        then
	    BOOST_LIB_SUFFIX="-mt"
            BOOST_LIB_DIR="$DIR"
            break
	fi
        if has_libboost "$DIR"
        then
	    BOOST_LIB_SUFFIX=""
            BOOST_LIB_DIR="$DIR"
            break
	fi
    fi
done

if [ -z "$BOOST_LIB_DIR" ]
then
    echo "ERROR: BOOST headers found, but not the required BOOST libs ($SUBLIBS)   $SCRIPT_NAME"   > /dev/stderr
    echo > /dev/stderr
  exit 1
fi


# Create symlinks to the libraries found, and put the names in BOOST_LDLIBS
BOOST_LDLIBS=""
for sublib in $SUBLIBS
do
    BOOST_LIB_ORIG=$(/bin/ls "$BOOST_LIB_DIR/libboost_$sublib$BOOST_LIB_SUFFIX".* | head -1)  # prefers blah.a over blah.so ;-)
    #    BOOST_LIB_EXTN=$("$SCRIPT_DIR/extn.sh" "$BOOST_LIB_ORIG")
    BOOST_LIB_NAME=$(basename "$BOOST_LIB_ORIG")
    BOOST_LIB_EXTN=$(echo "$BOOST_LIB_NAME" | sed -e "s/libboost_$sublib$BOOST_LIB_SUFFIX.//" )
    case "$BOOST_LIB_EXTN" in
	so.*) BOOST_LIB_EXTN=so ;;
    esac
    # case "$BOOST_LIB_ORIG" in
    # 	*.a) BOOST_LIB_EXTN=a ;;
    # 	*.so) BOOST_LIB_EXTN=so ;;
    # 	*) echo "ERROR: unknown library suffix in $BOOST_LIB_ORIG" > /dev/stderr; exit 2 ;;
    # 	esac
    /bin/ln -s "$BOOST_LIB_ORIG"  "$EXTLIB_5_DIR_FULL/lib/libboost_$sublib-symlink.$BOOST_LIB_EXTN"
    BOOST_LDLIBS="$BOOST_LDLIBS  -lboost_$sublib-symlink"
done


# Hack for MinGW as it needs libws2_32 for Boost ASIO
MINGW=`uname -s | grep -F MINGW`
if [ -n "$MINGW" ]
then
  BOOST_LDLIBS="$BOOST_LDLIBS  -lws2_32"
fi



# Print two assignments -- strictly do not need BOOST_LIB_DIR, but info is useful (see configure script)
echo "BOOST_LIB_DIR=\"$BOOST_LIB_DIR\""
echo "BOOST_LDLIBS=\"$BOOST_LDLIBS\""
#echo "$BOOST_LDLIBS"
exit 0
