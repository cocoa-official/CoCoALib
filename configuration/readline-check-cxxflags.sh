#! /bin/bash

SCRIPT_NAME=[[`basename "$0"`]]
SCRIPT_DIR=`dirname "$0"`

# Expects args: include paths and libs for READLINE
# expects env variables CXX and CXXFLAGS inherited from parent shell.

# Check READLINE lib is compatible with CXXFLAGS (from GMP)
# Exit code is 0 if compatible, and output is one of -ltermcap, -lncurses, -lcurses.
# Exit code is 1 is not compatible (no mesg is output)
# Exit code is 2 is input was bad (a diagnostic is output on /dev/stderr).

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


if [ $# -ne 2 ]
then
    echo "ERROR: expected 2 args (abs paths of readline header and libreadline)   $SCRIPT_NAME" > /dev/stderr
    exit 2
fi

READLINE_HDR="$1"
READLINE_LIB="$2"
# The following is a cryptic if...then block
is_absolute "$READLINE_HDR" || is_absolute "$READLINE_LIB" ||
(
  echo "ERROR: args must be absolute paths (readline header and libreadline)   $SCRIPT_NAME"   > /dev/stderr
  exit 1
)


READLINE_HDR_DIR=`dirname "$1"`
READLINE_HDR_DIR_DIR=`dirname "$READLINE_HDR_DIR"`

if [ -z "$CXX" ]
then
  echo "ERROR: environment variable CXX must be set.   $SCRIPT_NAME"   > /dev/stderr
  exit 1
fi

# Create tmp directory, put test prog in it, compile and run.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=`mktempdir readline-check-cxxflags`

pushd "$TMP_DIR"  > /dev/null
/bin/cat > test-readline.c <<EOF
#include "stdlib.h"
#include "stdio.h"
#include "unistd.h"
#include "readline/readline.h"
#include "readline/history.h"

#include <string>

int main()
{
  char *line;
  line = readline("");
  std::string str(line);
  free(line);
  // We shall test on an input of 3 chars...
  if (str.size() != 3) return 1;
  return 0;
}
EOF

for TERMCAP in termcap ncurses curses
do
  echo "Trying TERMCAP=$TERMCAP"  >> LogFile
  $CXX $CXXFLAGS -I"$READLINE_HDR_DIR_DIR" test-readline.c -o test-readline "$READLINE_LIB" -l$TERMCAP  >> LogFile 2>&1
  if [ $? = 0 ]
  then
    LIBTERMCAP=-l$TERMCAP
    break
  fi
done

if [ -z "$LIBTERMCAP" ]
then
  echo "ERROR: did not find suitable termcap/curses library   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE:  $TMP_DIR/LogFile   $SCRIPT_NAME"                           > /dev/stderr
  exit 1
fi

/bin/cat > test-readline.in <<EOF
abc
EOF

./test-readline < test-readline.in  >> LogFile  2>&1
if [ $? -ne 0 ]
then
  echo "ERROR: test program gave run-time error   $SCRIPT_NAME"  > /dev/stderr
  echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"                > /dev/stderr
  exit 1
fi

# Clean up TMP_DIR
popd  > /dev/null
/bin/rm -rf "$TMP_DIR"
echo "$LIBTERMCAP"
exit 0
