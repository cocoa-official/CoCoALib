#!/bin/bash

# This script looks for a BOOST installation in a standard location.
# If a single BOOST installation is found, it prints out the full path of
# the dir containing subdir "boost" and returns with exit code 0.
# If none is found, it prints out an appropriate message (on /dev/stderr)
# and exits with code 2 (used by "configure" in CoCoA root dir).
# If several are found, it prints out the first one then exits with code 1
# (used by "configure" in CoCoA root dir).

# EXIT CODES: 0 = "found 1 boost insallation";  prints full path
#             1 = "found several installations"; prints full path of one
#             2 = "found no boost installation"

SCRIPT_NAME=[[$(basename "$0")]]
SCRIPT_DIR=$(dirname "$0")

# Check there were no args
if [ $# -ne 0 ]
then
  echo "ERROR: expecting no args   $SCRIPT_NAME"  > /dev/stderr
  exit 2
fi

##################################################################
# List of most common directories in which boost header dir is normally found.
STD_BOOST_HDR_DIRS="/usr/local/include  /usr/include  /opt/local/include  /sw/include  /usr/sfw/include"



# We create a temp dir and work in there.
umask 22
source "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR=$(mktempdir boost-find-hdrs)
/bin/rm -rf "$TMP_DIR"  &&  /bin/mkdir -p "$TMP_DIR"
if [ $? -ne 0 ]
then
  echo "ERROR: failed to create temporary directory \"$TMP_DIR\"   $SCRIPT_NAME"   > /dev/stderr
  exit 2
fi

/bin/rm -rf "$TMP_DIR/BOOST-HDR-DIR-LIST"
COUNTER=0
for dir in $STD_BOOST_HDR_DIRS
do
  if [ -d "$dir/boost" ] && [ -r "$dir/boost" ]
  then
    echo "$dir" >> "$TMP_DIR/BOOST-HDR-DIR-LIST"
    COUNTER=$(( 1 + COUNTER ))
  fi
done

if [ "$COUNTER" -eq 0 ]
then
  # Did not find any plausible BOOST installation, so return empty handed.
  echo "ERROR: No BOOST installation found; looked inside $STD_BOOST_HDR_DIRS   $SCRIPT_NAME"   > /dev/stderr
  exit 2
fi

if [ "$COUNTER" = 1 ]
then
    EXIT_CODE=0
else
    EXIT_CODE=1
    echo "INFO: Found $COUNTER BOOST installations:   $SCRIPT_NAME"   > /dev/stderr
    sed -e "s/^/INFO: /" "$TMP_DIR/BOOST-HDR-DIR-LIST"                > /dev/stderr
##??    /bin/cat  "$TMP_DIR/BOOST-HDR-DIR-LIST"                              > /dev/stderr
#    echo "INFO: SELECTING THE FIRST ONE"                            > /dev/stderr
fi

# We have found at least one suitable directory (it exists & is readable,
# and contains readable subdir called "boost"), so pick first one.
BOOST_HDR_DIR=$(head -1 "$TMP_DIR/BOOST-HDR-DIR-LIST")
/bin/rm -rf "$TMP_DIR"
echo "$BOOST_HDR_DIR"
exit $EXIT_CODE
