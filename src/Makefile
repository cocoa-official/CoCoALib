# Makefile for CoCoALib/src/ directory.

COCOA_ROOT=..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/

.PHONY: default
default:
	@. $(COCOA_ROOT)/configuration/shell-fns.sh; \
	 echobox "ERROR: please run make from the CoCoALib root directory"

# The subdirectories containing library source code.
SRCDIRS=AlgebraicCore 
# Subdir parser/ is not yet ready to be built.
# Subdir tests/ is handled specially.
APPDIRS=server CoCoA-5

# This target will be built when called from "make" the CoCoALib root directory.
.PHONY: all
all: library
	@. $(COCOA_ROOT)/configuration/shell-fns.sh; \
	 for subdir in $(APPDIRS); \
	 do \
	   ( cd $$subdir; $(MAKE) -q all ); \
	   if [ $$? -ne 0 ] ; \
	   then \
	     echounderline "Compiling application subdirectory $$subdir/"; \
	     ( cd $$subdir; $(MAKE) -s all ); \
	     if [ $$? -ne 0 ] ;\
	     then \
	       echoerror "Compilation failed in application subdirectory $$subdir/"; \
	       exit 1; \
	     fi; \
	   fi; \
	 done; \
	 echounderline "Compiled all application subdirectories."

.PHONY: lib
lib: library

.PHONY: library
library:
	@. $(COCOA_ROOT)/configuration/shell-fns.sh; \
	 mkdir -p $(COCOA_ROOT)/lib; touch $(COCOA_ROOT)/lib/timestamp; \
	 if [ -f $(COCOA_LIB) -a $(COCOA_LIB) -ot $(COCOA_ROOT)/configuration/autoconf.mk ]; \
	 then \
	   /bin/rm -f $(COCOA_LIB); \
	 fi; \
	 for subdir in $(SRCDIRS); \
	 do \
	   ( cd $$subdir; $(MAKE) -q all ); \
	   if [ $$? -ne 0 ] ; \
	   then \
	     echounderline "Compiling CoCoALib source subdirectory $$subdir/"; \
	     ( cd $$subdir; $(MAKE) -s all ); \
	     if [ $$? -ne 0 ] ;\
	     then \
	       echoerror "Compilation failed in CoCoALib source subdirectory $$subdir/"; \
	       exit 1; \
	     fi; \
	   fi; \
	 done; \
	 if [ $(COCOA_ROOT)/lib/timestamp -nt $(COCOA_LIB) ]; \
	 then \
	   echo "[CoCoALib already compiled]"; \
	   echo; \
	 else \
	   ar -s $(COCOA_LIB); \
	   echobox "Compilation of CoCoALib completed."; \
	 fi


.PHONY: check
check: library
	@. $(COCOA_ROOT)/configuration/shell-fns.sh; \
	 echo; \
	 echounderline "Compiling and running all CoCoALib tests"
	 cd tests; $(MAKE) -s check


.PHONY: clean  clean-local  clean-subdirs
clean: clean-local  clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f  ./*~  ./.*~  ./.\#*

clean-subdirs:
	@for subdir in $(SRCDIRS) $(APPDIRS) tests; \
	 do \
	   ( cd $$subdir; $(MAKE) -s clean ); \
	 done

.PHONY: veryclean  veryclean-subdirs
veryclean: clean-local veryclean-subdirs
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-subdirs:
	@for subdir in $(SRCDIRS) $(APPDIRS) tests; \
	 do \
	   ( cd $$subdir; $(MAKE) -s veryclean ); \
	 done


.PHONY: dependencies
dependencies:
	@for subdir in $(SRCDIRS) $(APPDIRS) tests; \
	 do \
	   ( cd "$$subdir"; $(MAKE) -s dependencies ); \
	 done
