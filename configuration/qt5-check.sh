#!/bin/bash

# This script checks sufficiency of installed Qt5 system for the CoCoA5 GUI.
# Exit code is 0 if sufficient, o/w non-zero.

# This script expects no arguments.

SCRIPT_NAME="[[$(basename "$0")]]"
SCRIPT_DIR="$(dirname "$0")"

if [ $# -ne 0 ]; then
    echo "ERROR: expected no args   $SCRIPT_NAME" > /dev/stderr
    exit 1
fi

# Create tmp directory, put test prog in it, compile and run.
umask 22
. "$SCRIPT_DIR/shell-fns.sh"
TMP_DIR="$(mktempdir qt5-check)"

pushd "$TMP_DIR" > /dev/null || exit 255

# Create project file for qmake
/bin/cat > qt5-check.pro << 'EOF'
QMAKE_MAKEFILE = qt5-check-makefile

TEMPLATE = app
TARGET = qt5-check
QT += xml webkitwidgets printsupport # crucial line with dependencies
SOURCES += qt5-check.cpp
EOF

# Create C++ file
/bin/cat > qt5-check.cpp << 'EOF'
#include <QtWidgets/QApplication> // a header file necessary for QCodeEdit
#include <iostream>

int main(int,char**) {
    std::cout << "Qt5 version " << qVersion() << " found" << std::endl;
    return 0;
}
EOF

# Set qmake flag for macOS (see src/CoCoA-5/make-c5makefile.sh
#                           and src/CoCoA-5/make-qcodeeditmakefile.sh)
if [ "$(uname)" = "Darwin" ]
then
    DARWIN_OPTS="-spec macx-g++"
fi

# Check for qmake
echo "# Checking for qmake-qt5 or qmake..."  > LogFile
if ! QMAKE=$(   ( command -v qmake-qt5 2>> LogFile ) \
		    || ( command -v qmake     2>> LogFile ) )
then
    # Deliberately leave $TMP_DIR to assist debugging.
    echo "ERROR: qmake (and qmake-qt5) not found  $SCRIPT_NAME" > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"             > /dev/stderr
    exit 3
fi

# Basic check of $QMAKE
echo "$QMAKE -v" >> LogFile
if ! "$QMAKE" -v >> LogFile 2>&1
then
    # Deliberately leave $TMP_DIR to assist debugging.
    echo "ERROR: $QMAKE not working properly: does not report version number    $SCRIPT_NAME" > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile    $SCRIPT_NAME"                                          > /dev/stderr
    exit 4
fi

# Use $QMAKE to create makefile
echo "$QMAKE $DARWIN_OPTS qt5-check.pro" >> LogFile
if ! "$QMAKE" $DARWIN_OPTS qt5-check.pro >> LogFile 2>&1
then
    # Deliberately leave $TMP_DIR to assist debugging.
    echo "ERROR: $QMAKE exited abnormally   $SCRIPT_NAME"     > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"           > /dev/stderr
    exit 5
fi

# If MAKE is unset, try to set it reasonably
if [ -z "$MAKE" ]
then
    echo "Warning: MAKE not set, trying to find it   $SCRIPT_NAME" >> LogFile
    if ! MAKE="$(command -v make 2>> LogFile)"
    then
        # Deliberately leave $TMP_DIR to assist debugging.
        echo "ERROR: make not found   $SCRIPT_NAME"     > /dev/stderr
        echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME" > /dev/stderr
        exit 6
    fi
fi

# Use $MAKE to compile cpp file
echo "Using MAKE = $MAKE   $SCRIPT_NAME" >> LogFile
echo "$MAKE -f qt5-check-makefile" >> LogFile
if ! "$MAKE" -f qt5-check-makefile >> LogFile 2>&1
then
    # Deliberately leave $TMP_DIR to assist debugging.
    echo "ERROR: $MAKE exited abnormally   $SCRIPT_NAME"     > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"          > /dev/stderr
    exit 7
fi

# Find Qt5 program - on some platforms it is located in release/
if [ -e "release/qt5-check" ]
then
    QT_EXECUTABLE="release/qt5-check"
else
    # Otherwise, it should be in the current directory
    QT_EXECUTABLE="./qt5-check"
fi


# Run Qt5 program (since we succeeded in building it)
echo "$QT_EXECUTABLE" >> LogFile
if ! $QT_EXECUTABLE >> LogFile 2>&1
then
    # Deliberately leave $TMP_DIR to assist debugging.
    echo "ERROR: qt5-check exited abnormally   $SCRIPT_NAME"      > /dev/stderr
    echo "LOGFILE: $TMP_DIR/LogFile   $SCRIPT_NAME"               > /dev/stderr
    exit 8
fi

# Clean up $TMP_DIR
popd > /dev/null  ||  exit 255
/bin/rm -rf "$TMP_DIR"

echo "$QMAKE"
