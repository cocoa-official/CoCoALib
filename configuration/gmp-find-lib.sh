#!/bin/bash

SCRIPT_NAME=[[$(basename "$0")]]

# This script looks for a GMP library in a standard location.
# If a single GMP library is found, it prints out the full path of
# the (static) library and returns with exit code 0.
# If  several are found, it prints out an "info" message (on stderr),
# and prints the full path one of them, and returns with exit code 1.
# If none is found, it prints out an error message (on stderr)
# and returns with exit code 2.

# EXIT CODES: 0 = "found 1 GMP installation";  prints full path
#             1 = "found several installations"; prints full path of one
#             2 = "found no GMP installation"


SCRIPT_NAME=[[$(basename "$0")]]
#??? SCRIPT_DIR=$(dirname "$0")

# Check there were no args
if [ $# -ne 0 ]
then
  echo "ERROR: expecting no args   $SCRIPT_NAME"  > /dev/stderr
  exit 2
fi

##################################################################
# Use find to search through various standard directories.
# NB look through all directories, even if a GMP has already been found.

# List of directories under which libgmp.a and/or libgmp.so is normally found.
ARCH=$(uname -m)      #so far seen:  x86_64, aarch64, i686, i386
STD_GMP_LIBDIRS="/usr/local/lib  /usr/lib  /usr/lib/$ARCH-linux-gnu   /usr/lib64  /usr/lib32  /usr/lib/i386-linux-gnu  /opt/homebrew/lib  /opt/local/lib  /sw/lib  /usr/sfw/lib"


LIBGMPPATHS=libgmp-paths
/bin/rm -rf "$LIBGMPPATHS"
COUNTER=0
for dir in $STD_GMP_LIBDIRS
do
  if [ -d "$dir" ]
  then
      # BUG? check also that the files are readable???
    if [ -f "$dir/libgmp.a" ];      then echo "$dir/libgmp.a"  >> "$LIBGMPPATHS";     COUNTER=$((1 + COUNTER)); continue; fi
    if [ -f "$dir/libgmp.so" ];     then echo "$dir/libgmp.so" >> "$LIBGMPPATHS";     COUNTER=$((1 + COUNTER)); continue; fi
    if [ -f "$dir/libgmp.dll.a" ];  then echo "$dir/libgmp.dll.a"  >> "$LIBGMPPATHS"; COUNTER=$((1 + COUNTER)); continue; fi
  fi
done

if [ "$COUNTER" -eq 0 ]
then
  # Did not find any plausible GMP installation, so return empty handed.
    echo "ERROR: No GMP installation found; looked inside $STD_GMP_LIBDIRS   $SCRIPT_NAME"   > /dev/stderr
    echo ">>>>> HINT:  try installing the linux package libgmp-dev  <<<<<"  >> /dev/stderr
    echo ">>>>> HINT:  or get sources from https://www.gmplib.org/  <<<<<"  >> /dev/stderr
#    /bin/rm -f "$LIBGMPPATHS"  # never created the file -- no need to remove it!
  exit 2
fi


if [ "$COUNTER" -eq 1 ]
then
    EXIT_CODE=0
else
    EXIT_CODE=1
  echo "INFO:  Found $COUNTER GMP libraries   $SCRIPT_NAME"   > /dev/stderr
  sed -e "s/^/INFO: /"  "$LIBGMPPATHS"                        > /dev/stderr
fi

# We have found >= 1 files called libgmp.a or libgmp.so;
# pick the first one, and print it (it is a full path).
GMP_LIB=$(head -1 "$LIBGMPPATHS")
/bin/rm  "$LIBGMPPATHS"
echo "$GMP_LIB"
exit "$EXIT_CODE"
