#! /bin/bash

# This script should be called from the Makefile in CoCoALib/src/CoCoA-5/ directory.

# Do some finishing if we are on MacOSX; nothing to do on other platforms.

SYS_TYPE=`uname`

if [ "$SYS_TYPE" = "Darwin" ]
then
  mv C5.app/Contents/MacOS/C5 C5.app/Contents/MacOS/C5.bin
  cp C5.mac.script C5.app/Contents/MacOS/C5
  PKG_DIR=C5.app/Contents/MacOS/packages
  MAN_DIR=C5.app/Contents/MacOS/CoCoAManual
  /bin/rm -rf $PKG_DIR $MAN_DIR
  mkdir $PKG_DIR $MAN_DIR
  cp -r packages/*.cpkg5  $PKG_DIR
  cp packages/init.cocoa5  $PKG_DIR
  cp -r CoCoAManual/CoCoAHelp.xml  $MAN_DIR
  rm -rf C5.app/Contents/Frameworks/
  mkdir  C5.app/Contents/Frameworks
  cp -R /Library/Frameworks/QtXml.framework C5.app/Contents/Frameworks
  cp -R /Library/Frameworks/QtGui.framework C5.app/Contents/Frameworks
  cp -R /Library/Frameworks/QtCore.framework C5.app/Contents/Frameworks
fi
