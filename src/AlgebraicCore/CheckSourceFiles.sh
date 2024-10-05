#! /bin/bash

if [ $# = 0 ]
then
  echo "ERROR: $0 expects args (the names of the CoCoALib C++ source files)"
  exit 1
fi

OLD_SRCS="$*"
IGNORE_SRCS="leak_checker.C debug_new.C"

TMPFILE1=.old_srcs
TMPFILE2=.new_srcs

/bin/rm -f "$TMPFILE1"
/bin/rm -f "$TMPFILE2"
for file in $OLD_SRCS
do
  echo "$file" >> "$TMPFILE2"
done
for file in $IGNORE_SRCS
do
  echo "$file" >> "$TMPFILE2"
done
sort "$TMPFILE2" > "$TMPFILE1"

/bin/ls -d ./*.C | cut -b 3- > "$TMPFILE2"

cmp "$TMPFILE1" "$TMPFILE2" > /dev/null
if [ $? = 0 ]
then
  /bin/rm "$TMPFILE1" "$TMPFILE2"
  exit 0
fi


NEW_SRCS=`diff "$TMPFILE2" "$TMPFILE1" | grep "^<" | tr -d "<"`
LOST_SRCS=`diff "$TMPFILE2" "$TMPFILE1" | grep "^>" | tr -d ">"`
/bin/rm -f $TMPFILE1 $TMPFILE2

echo "ERROR:"                                                          >/dev/stderr
if [ -n "$LOST_SRCS" ]
then
  echo "ERROR:  ***LOST SOURCE FILES*** "                              >/dev/stderr
  for lostfile in $LOST_SRCS
  do
    echo "ERROR: --> $lostfile"                                        >/dev/stderr
  done
  echo "ERROR:"                                                        >/dev/stderr
  echo "ERROR:  Please recover the lost files!"                        >/dev/stderr
  echo "ERROR:  (or edit SRCS in CoCoALib/src/AlgebraicCore/Makefile)" >/dev/stderr
  echo "ERROR:"                                                        >/dev/stderr
fi

if [ -n "$NEW_SRCS" ]
then
  echo "ERROR:  ***NEW SOURCE FILES*** "                                      >/dev/stderr
  for newfile in $NEW_SRCS
  do
    echo "ERROR: --> $newfile"                                                >/dev/stderr
  done
  echo "ERROR:"                                                               >/dev/stderr
  echo "ADVICE: If you want to add these files to the CoCoALib sources"       >/dev/stderr
  echo "ADVICE: you must edit CoCoALib/src/AlgebraicCore/Makefile, and"       >/dev/stderr
  echo "ADVICE: add the new name(s) to the variable SRCS; you will"           >/dev/stderr
  echo "ADVICE: find it on line 15."                                          >/dev/stderr
  echo "ADVICE: Otherwise move the files elsewhere (to examples directory?)"  >/dev/stderr
fi

exit 2
