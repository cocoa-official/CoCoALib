#!/bin/bash

# Various shell functions used in some of the Makefiles and other
# scripts included with CoCoALib.

# Copyright 2006,2023 John Abbott.
# You are free to use any part of this code in your own programs.

# tempdir name must depend on user in case two users simultaneously
# want to run the configuration scripts.
# BUG (not serious): suboptimal if configure is run just before midnight!
mktempdir()
{
    TODAY=`date "+%Y%m%d"`
    TIME=`date "+%H%M%S"`
    TMP_DIR="/tmp/CoCoALib-config-$USER/$TODAY/$1-$TIME-$$"
    /bin/rm -rf "$TMP_DIR"  &&  /bin/mkdir -p "$TMP_DIR"
    if [ $? -ne 0 ]
    then
	echo "ERROR: failed to create temporary directory \"$TMP_DIR\"   $SCRIPT_NAME"   > /dev/stderr
	exit 1
    fi
    echo "$TMP_DIR"
}


echounderline()
{
  echo "$*"
  echo "$*" | tr "\040-\377" "[-*]"
}

echobox()
{
  mesg=">>>>  $*  <<<<"
  dashes=`echo "$mesg" | tr "\040-\377" "[-*]"`
  echo "$dashes"
  echo "$mesg"
  echo "$dashes"
}

echoerror()
{
  mesg=">>>>>  $*  <<<<<"
  equals=`echo "$mesg" | tr "\040-\377" "[=*]"`
  echo "$equals"
  echo "$mesg"
  echo "$equals"
}

# # https://stackoverflow.com/questions/3915040/how-to-obtain-the-absolute-path-of-a-file-via-shell-bash-zsh-sh
# function abspath() {
#     # generate absolute path from relative path
#     # $1     : relative filename
#     # return : absolute path
#     if [ -d "$1" ]; then
#         # dir
#         (cd "$1"; pwd)
#     elif [ -f "$1" ]; then
#         # file
#         if [[ $1 = /* ]]; then
#             echo "$1"
#         elif [[ $1 == */* ]]; then
#             echo "$(cd "${1%/*}"; pwd)/${1##*/}"
#         else
#             echo "$(pwd)/$1"
#         fi
#     fi
# }
