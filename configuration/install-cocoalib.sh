#! /bin/bash

# EXPECT INPUTS
# COCOALIB_INSTALL_DIR  (will be put into .../include and .../lib)
# INSTALL_CMD           (cp or install, probably)


# This script assumes it is run from the CoCoALib ROOT directory, and that libcocoa.a has been built.

if [ $? -ne 2 ]
then
    echo "$0: expects 2 args: COCOALIB_INSTALL_DIR & INSTALL_CMD"  > /dev/stderr
    exit 1
fi

COCOALIB_INSTALL_DIR=$1
INSTALL_CMD=$2

which "$INSTALL_CMD" > /dev/null
if [ $? -ne 0 ]
then
    echo "$0: ERROR: \`$INSTALL_CMD\' is not a recognized command"  > /dev/stderr
    exit 1
fi

echo "======================================================================"
echo ">>> WARNING  CoCoALib installation procedure is still PRELIMINARY! <<<"
echo "======================================================================"
echo
echo "Continuing with installation after 5 secs..."
echo
sleep 5

# Ad hoc check that we are in CoCoA root directory, and that libcocoa.a has been built.
if [ \! -f lib/libcocoa.a ];
then
    echo "***** INSTALLATION ERROR: CoCoALib is not built!           *****";
    echo "***** Please run \"make library doc\" before installation.   *****";
    echo;
    exit 1;
fi

if [ \! -f doc/CoCoALib.pdf -o \! -f examples/index.html ];
then
    echo "***** INSTALLATION ERROR: CoCoALib documentation is missing! *****";
    echo "***** Please run \"make library doc\" before installation.     *****";
    echo;
    exit 1;
fi

if [ \! -d "$(COCOALIB_INSTALL_DIR)/include" -o \! -w "$(COCOALIB_INSTALL_DIR)/include" ];
then
    echo;
    echo "***** ERROR: Installation directory \"$(COCOALIB_INSTALL_DIR)/include\" is not writable *****";
    echo "***** >>>>  Consider using \"sudo\" command  <<<<";
    echo;
    exit 1;
fi

if [ \! -d "$(COCOALIB_INSTALL_DIR)/lib" -o \! -w "$(COCOALIB_INSTALL_DIR)/lib" ];
then
    echo;
    echo "***** ERROR: Installation directory \"$(COCOALIB_INSTALL_DIR)/lib\" is not writable *****";
    echo "***** >>>>   Consider using \`sudo' command   <<<<";
    exit 1;
fi

if [ -e "$(COCOALIB_INSTALL_DIR)/include/CoCoA" -a \! -L "$(COCOALIB_INSTALL_DIR)/include/CoCoA" ];
then
    echo;
    echo "***** ERROR: $(COCOALIB_INSTALL_DIR)/include/CoCoA exists but is not a symlink *****";
    echo "***** >>>> Please remove it or rename it before installing CoCoALib <<<<";
    exit 2;
fi

if [ -e "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a" -a \! -f "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a" ];
then
    echo;
    echo "***** ERROR: $(COCOALIB_INSTALL_DIR)/lib/libcocoa.a exists but is not a file *****";
    echo "***** >>>> Please remove it or rename it before installing CoCoALib <<<<";
    exit 2;
fi

if [ -e "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)" ];
then
    echo ">>>>  ?? CoCoALib ALREADY INSTALLED ??  <<<<";
    /bin/ls -ld "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)";
    echo;
    read -p "Really overwrite existing installation? " yn;
    if [ "X$$yn" \!= "Xy" -a "X$$yn" \!= "Xyes" ]; then exit 3; fi;
fi

/bin/rm -rf "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"
/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"
$(INSTALL_CMD) -m 644 include/CoCoA/*.H "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"
echo "COCOA_ROOT=$(COCOALIB_INSTALL_DIR)" > "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/Makefile"
/bin/cat configuration/autoconf.mk >> "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/Makefile"
/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/examples"
$(INSTALL_CMD) -m 644 examples/ex-*.C "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/examples"
$(INSTALL_CMD) -m 644 examples/index.html "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/examples"
echo "Installed CoCoA examples in \"$(COCOALIB_INSTALL_DIR)/include/CoCoA/examples/\""
/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc"
$(INSTALL_CMD) -m 644 doc/CoCoALib.pdf "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc"
/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc/html"
$(INSTALL_CMD) -m 644 doc/html/*.html "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc/html"
$(INSTALL_CMD) -m 644 doc/html/cocoalib-doc.css "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc/html"
echo "Installed CoCoA documentation in \"$(COCOALIB_INSTALL_DIR)/include/CoCoA/doc/\""
/bin/rm -f "$(COCOALIB_INSTALL_DIR)/include/CoCoA"
/bin/ln -s "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)" "$(COCOALIB_INSTALL_DIR)/include/CoCoA"
echo "Installed CoCoA headers in \"$(COCOALIB_INSTALL_DIR)/include/CoCoA/\""
/bin/rm -f "$(COCOALIB_INSTALL_DIR)/lib/libcocoa-$(COCOALIB_VERSION).a"
$(INSTALL_CMD) -m 644 lib/libcocoa.a         "$(COCOALIB_INSTALL_DIR)/lib/libcocoa-$(COCOALIB_VERSION).a"
/bin/rm -f "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a"
/bin/ln -s "$(COCOALIB_INSTALL_DIR)/lib/libcocoa-$(COCOALIB_VERSION).a" "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a"
echo "Installed CoCoA library in \"$(COCOALIB_INSTALL_DIR)/lib/\""
echo
if [ "$(HAVE_BOOST)" = "yes" ]; then echo "**IMPORTANT** To install also CoCoA-5:  cd src/CoCoA-5; make install"; fi
