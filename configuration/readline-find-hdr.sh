#! /bin/bash

SCRIPT_NAME=[[`basename "$0"`]]

# This script expects the (full) path of libreadline.a (or the string "-lreadline").
# If the arg is "-lreadline", the script simply exits with code 0.
# It prints out the (full) path of (hopefully) the corresponding readline.h.

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
  echo "ERROR: expected 1 arg (full path of libreadline)   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi

READLINE_LIB="$1"

# If READLINE is in std directory, assume header is too.  Output nothing.
if [ "X$READLINE_LIB" = "X-lreadline" ]
then
  exit 0
fi

# The following is a cryptic if...then block
is_absolute "$READLINE_LIB" ||
(
  echo "ERROR: arg is not an absolute path: \"$READLINE_LIB\"   $SCRIPT_NAME"  > /dev/stderr
  exit 1
)

if [ \! -f "$READLINE_LIB" -o \! -r "$READLINE_LIB" ]
then
  echo "ERROR: specified READLINE library is unreadable $READLINE_LIB   $SCRIPT_NAME"  > /dev/stderr
  exit 1
fi


# Usual paths are  /usr/lib/libreadline.a   or  /usr/lib/platform/libreadline.a
# We try both cases, hoping to find /usr/include/readline/readline.h
# [BUG: really need a better algorithm]

READLINE_LIB_DIR=`dirname "$READLINE_LIB"`
READLINE_LIB_DIR_DIR=`dirname "$READLINE_LIB_DIR"`
READLINE_LIB_DIR_DIR_DIR=`dirname "$READLINE_LIB_DIR_DIR"`
READLINE_HDR_DIR1="$READLINE_LIB_DIR_DIR/include/readline"
### On some systems there is an "intermediate" directory:
INTERMEDIATE_DIR="$(basename "$READLINE_LIB")"
INTERMEDIATE_DIR="${INTERMEDIATE_DIR%.*}"
INTERMEDIATE_DIR="${INTERMEDIATE_DIR#lib}"
READLINE_HDR_DIR2="$READLINE_LIB_DIR_DIR/include/$INTERMEDIATE_DIR/readline"
###
READLINE_HDR_DIR3="$READLINE_LIB_DIR_DIR_DIR/include/readline"
READLINE_HDR1="$READLINE_HDR_DIR1/readline.h"
READLINE_HDR2="$READLINE_HDR_DIR2/readline.h"
READLINE_HDR3="$READLINE_HDR_DIR3/readline.h"

if [ -f "$READLINE_HDR1" -a -r "$READLINE_HDR1" ]
then
  # We've found a plausible readline.h.
  echo "$READLINE_HDR1"
  exit 0
fi

if [ -f "$READLINE_HDR2" -a -r "$READLINE_HDR2" ]
then
  # We've found a plausible readline.h.
  echo "$READLINE_HDR2"
  exit 0
fi

if [ -f "$READLINE_HDR3" -a -r "$READLINE_HDR3" ]
then
  # We've found a plausible readline.h.
  echo "$READLINE_HDR3"
  exit 0
fi

echo "ERROR: Trouble finding READLINE header file:    $SCRIPT_NAME"  > /dev/stderr
echo "ERROR: + tried  $READLINE_HDR1   $SCRIPT_NAME"                 > /dev/stderr
echo "ERROR: + tried  $READLINE_HDR2   $SCRIPT_NAME"                 > /dev/stderr
echo "ERROR: + tried  $READLINE_HDR3   $SCRIPT_NAME"                 > /dev/stderr

exit 2
