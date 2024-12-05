# Makefile for CoCoALib root directory.

# This line will cause make to complain if you run make
# before running configure.
include configuration/autoconf.mk

.PHONY: default
default: check-platform
	@echo "======================================================="
	@echo "Compiling CoCoALib-$(COCOALIB_VERSION)"
	@echo "  CXX = $(CXX)"
	@echo "  CXXFLAGS = $(CXXFLAGS)"
	@echo "  CXXFLAGS_DEFINES = $(CXXFLAGS_DEFINES)"
	@echo "======================================================="
	@$(MAKE) all


.PHONY: check-platform
check-platform:
	@PLATFORM="`uname -s -r -m`"; if [ "$$PLATFORM" \!= "$(PLATFORM)" ]; then . configuration/shell-fns.sh; echobox "ERROR: new platform!  Please run \"configure\" script."; exit 1; fi


.PHONY: install
install:
	@echo "======================================================================"
	@echo ">>> WARNING  CoCoALib installation procedure is still PRELIMINARY! <<<"
	@echo "======================================================================"
	@echo
	@echo "Continuing with installation after 5 secs..."
	@echo
	@sleep 5
	@if [ \! -f lib/libcocoa.a ]; \
	 then \
	   echo "***** INSTALLATION ERROR: CoCoALib is not built!           *****"; \
	   echo "***** Please run \"make library doc\" before installation.   *****"; \
	   echo; \
	   exit 1; \
	  fi
	@if [ \! -f doc/CoCoALib.pdf -o \! -f examples/index.html ]; \
	 then \
	   echo "***** INSTALLATION ERROR: CoCoALib documentation is missing! *****"; \
	   echo "***** Please run \"make library doc\" before installation.     *****"; \
	   echo; \
	   exit 1; \
	  fi
	@if [ \! -d "$(COCOALIB_INSTALL_DIR)/include" -o \! -w "$(COCOALIB_INSTALL_DIR)/include" ]; \
	 then \
	   echo; \
	   echo "***** ERROR: Installation directory \"$(COCOALIB_INSTALL_DIR)/include\" is not writable *****"; \
	   echo "***** >>>>  Consider using \"sudo\" command  <<<<"; \
	   echo; \
	   exit 1; \
	 fi
	@if [ \! -d "$(COCOALIB_INSTALL_DIR)/lib" -o \! -w "$(COCOALIB_INSTALL_DIR)/lib" ]; \
	 then \
	   echo; \
	   echo "***** ERROR: Installation directory \"$(COCOALIB_INSTALL_DIR)/lib\" is not writable *****"; \
	   echo "***** >>>>   Consider using \`sudo' command   <<<<"; \
	   exit 1; \
	 fi
	@if [ -e "$(COCOALIB_INSTALL_DIR)/include/CoCoA" -a \! -L "$(COCOALIB_INSTALL_DIR)/include/CoCoA" ]; \
	then \
	  echo; \
	  echo "***** ERROR: $(COCOALIB_INSTALL_DIR)/include/CoCoA exists but is not a symlink *****"; \
	  echo "***** >>>> Please remove it or rename it before installing CoCoALib <<<<"; \
	  exit 2; \
	fi
	@if [ -e "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a" -a \! -f "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a" ]; \
	then \
	  echo; \
	  echo "***** ERROR: $(COCOALIB_INSTALL_DIR)/lib/libcocoa.a exists but is not a file *****"; \
	  echo "***** >>>> Please remove it or rename it before installing CoCoALib <<<<"; \
	  exit 2; \
	fi
	@if [ -e "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)" ]; \
	 then \
	   echo ">>>>  ?? CoCoALib ALREADY INSTALLED ??  <<<<"; \
	   /bin/ls -ld "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"; \
	   echo; \
	   read -p "Really overwrite existing installation? " yn; \
	   if [ "X$$yn" \!= "Xy" -a "X$$yn" \!= "Xyes" ]; then exit 3; fi; \
	 fi
	@/bin/rm -rf "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"
	@/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"
	@$(INSTALL_CMD) -m 644 include/CoCoA/*.H "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)"
	@echo "COCOA_ROOT=$(COCOALIB_INSTALL_DIR)" > "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/Makefile"
	@/bin/cat configuration/autoconf.mk >> "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/Makefile"
	@/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/examples"
	@$(INSTALL_CMD) -m 644 examples/ex-*.C "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/examples"
	@$(INSTALL_CMD) -m 644 examples/index.html "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/examples"
	@echo "Installed CoCoA examples in \"$(COCOALIB_INSTALL_DIR)/include/CoCoA/examples/\""
	@/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc"
	@$(INSTALL_CMD) -m 644 doc/CoCoALib.pdf "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc"
	@/bin/mkdir -p "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc/html"
	@$(INSTALL_CMD) -m 644 doc/html/*.html "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc/html"
	@$(INSTALL_CMD) -m 644 doc/html/cocoalib-doc.css "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)/doc/html"
	@echo "Installed CoCoA documentation in \"$(COCOALIB_INSTALL_DIR)/include/CoCoA/doc/\""
	@/bin/rm -f "$(COCOALIB_INSTALL_DIR)/include/CoCoA"
	@/bin/ln -s "$(COCOALIB_INSTALL_DIR)/include/CoCoA-$(COCOALIB_VERSION)" "$(COCOALIB_INSTALL_DIR)/include/CoCoA"
	@echo "Installed CoCoA headers in \"$(COCOALIB_INSTALL_DIR)/include/CoCoA/\""
	@/bin/rm -f "$(COCOALIB_INSTALL_DIR)/lib/libcocoa-$(COCOALIB_VERSION).a"
	@$(INSTALL_CMD) -m 644 lib/libcocoa.a         "$(COCOALIB_INSTALL_DIR)/lib/libcocoa-$(COCOALIB_VERSION).a"
	@/bin/rm -f "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a"
	@/bin/ln -s "$(COCOALIB_INSTALL_DIR)/lib/libcocoa-$(COCOALIB_VERSION).a" "$(COCOALIB_INSTALL_DIR)/lib/libcocoa.a"
	@echo "Installed CoCoA library in \"$(COCOALIB_INSTALL_DIR)/lib/\""
	@echo
	@if [ "$(HAVE_BOOST)" = "yes" ]; then echo "**IMPORTANT** To install also CoCoA-5:  cd src/CoCoA-5; make install"; fi


.PHONY: all
all: check doc cocoa5-check examples server


# This target will be "built" only if you try to run make before
# having run configure (which creates autoconf.mk).
configuration/autoconf.mk: configuration/version
	@if [ -f configuration/autoconf.mk ]; \
	 then \
	   echo; \
	   echo "=========================================================================="; \
	   echo ">>>>> ERROR: Version has changed: please run   ./configure --again   <<<<<"; \
	   echo "=========================================================================="; \
	   echo; \
	   exit 1; \
	 fi
	@. configuration/shell-fns.sh; \
	 echoerror "ERROR: Cannot build CoCoALib: not yet configured."
	@echo
	@echo "====================================================="
	@echo ">>>            R U D E   M E S S A G E            <<<"
	@echo "====================================================="
	@echo ">>>  You must run configure before running make.  <<<"
	@echo ">>>  Please read the file README carefully.       <<<"
	@echo "====================================================="
	@exit 1


.PHONY: unified-header
unified-header:
	@cd include/CoCoA;  $(MAKE) -s

# This target will print out some harmless warning messages if
# you manage to delete the dependencies file in some subdirectory.
.PHONY: dependencies
dependencies:  unified-header
	@cd src;  $(MAKE) -s dependencies

# Note:  we are not not using "make -C" for compatibility with Solaris make
.PHONY: library
library:  dependencies
	@cd src;  $(MAKE) -s library

# Just an alias for library
.PHONY: lib
lib:  library

# Just an alias for library
.PHONY: cocoalib
cocoalib: library


.PHONY: check
check:  library
	@cd src;  $(MAKE) -s check


.PHONY: cocoa5
cocoa5:  library
	@if [ \! -f src/CoCoA-5/Main.C ]; \
	 then \
	   echo "--------------------------------------------"; \
	   echo ">>>>  CoCoA-5 not in this distribution  <<<<"; \
	   echo "--------------------------------------------"; \
	   exit; \
	 fi; \
	 if [ $(HAVE_BOOST) = "yes" ]; \
	 then \
	   cd src/CoCoA-5; \
	   $(MAKE) -s cocoa5; \
	 else \
	   echo "*********************************************************"; \
	   echo "*** Skipping CoCoA-5: configured with --only-cocoalib ***"; \
	   echo "*********************************************************"; \
	 fi

.PHONY: cocoa5-check
cocoa5-check:  library
	@if [ \! -f src/CoCoA-5/Main.C ]; \
	 then \
	   echo "--------------------------------------------"; \
	   echo ">>>>  CoCoA-5 not in this distribution  <<<<"; \
	   echo "--------------------------------------------"; \
	   exit; \
	 fi; \
	 if [ $(HAVE_BOOST) = "yes" ]; \
	 then \
	   cd src/CoCoA-5; \
	   $(MAKE) -s; \
	 else \
	   echo "*********************************************************"; \
	   echo "*** Skipping CoCoA-5: configured with --only-cocoalib ***"; \
	   echo "*********************************************************"; \
	 fi


.PHONY: server
server:  library
	@if [ \! -f src/server/CoCoAServer.C ]; \
	 then \
	   echo "------------------------------------------------"; \
	   echo ">>>>  CoCoAServer not in this distribution  <<<<"; \
	   echo "------------------------------------------------"; \
	 exit; \
	 fi; \
	 cd src/server; $(MAKE) -s all


.PHONY: benchmarks
benchmarks:  server
	@cd src/server/benchmarks; $(MAKE) -s benchmarks


.PHONY: doc
doc:
	@cd doc; $(MAKE) alldoc


.PHONY: examples
examples:  library
	@. configuration/shell-fns.sh; echounderline "Compiling the CoCoALib example programs..."
	@cd examples; $(MAKE) -s index.html executables
	@. configuration/shell-fns.sh; echounderline "Compilation of CoCoALib example programs completed; and index built."


.PHONY: clean clean-local clean-subdirs
clean:  clean-local  clean-subdirs
	@echo "Cleaned CoCoALib/"

clean-local:
	@/bin/rm -f  ./*~  ./.*~  ./.\#*
	@/bin/rm -rf lib/

clean-subdirs:
	@cd configuration; /bin/rm -f  ./*~  ./.*~  ./.\#*
	@cd doc;           $(MAKE) -s clean
	@cd examples;      $(MAKE) -s clean
	@cd include/CoCoA; $(MAKE) -s clean
	@cd src;           $(MAKE) -s clean


.PHONY: distclean
distclean:  veryclean
	@/bin/rm -rf configuration/last-config-cmd


.PHONY: veryclean veryclean-subdirs
veryclean:  clean-local  veryclean-subdirs
	@/bin/rm -rf SOURCE_RELEASE
	@/bin/rm -f configuration/autoconf.mk
	@echo "Verycleaned CoCoALib/"

veryclean-subdirs:
	@cd configuration; /bin/rm -f  ./*~  ./.*~  ./.\#*;  /bin/rm -rf  ExternalLibs/  ExternalLibs-CoCoA5/
	@cd doc;           $(MAKE) -s veryclean
	@cd examples;      $(MAKE) -s veryclean
	@cd include/CoCoA; $(MAKE) -s veryclean
	@cd src;           $(MAKE) -s veryclean
