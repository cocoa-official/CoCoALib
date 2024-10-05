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



# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/Makefile,v 1.24 2022/02/22 09:36:26 abbott Exp $
# $Log: Makefile,v $
# Revision 1.24  2022/02/22 09:36:26  abbott
# Summary: Inserted some spaces
#
# Revision 1.23  2019/10/02 17:59:10  abbott
# Summary: Replaced source command by more portable . (dot)
#
# Revision 1.22  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.21  2016/05/18 12:18:52  abbott
# Summary: Added new target "lib: synonym for library
#
# Revision 1.20  2014/07/28 14:45:27  abbott
# Summary: Improved target clean-local
# Author: JAA
#
# Revision 1.19  2014/03/14 11:02:07  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.18  2011/07/20 09:04:20  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.17  2011/07/15 16:53:46  abbott
# dependencies target now recurses into both SRCDIRS and APPDIRS.
#
# Revision 1.16  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.15  2010/10/08 14:11:41  abbott
# Makefile cleaning:
#   (a) proper handling of recursive make,
#   (b) better organized targets (to make embedded shell scripts simpler)
#
# Revision 1.14  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" by calls to "$(MAKE)".
#
# Revision 1.13  2010/07/09 17:07:35  abbott
# Added missing dependency on the library.
#
# Revision 1.12  2010/02/09 10:36:52  abbott
# Added missing semicolon.
# Removed dependency on "library" for target "all" (because it caused
# libcocoa.a to be ranlib'd twice).
#
# Revision 1.11  2009/05/21 09:46:57  abbott
# Added two progress messages.
#
# Revision 1.10  2009/05/20 14:38:56  abbott
# Ooops, added missing backslash.
#
# Revision 1.9  2009/05/20 14:30:17  abbott
# Changed 2 error messages, and a progress message.
# Added indexing command for CoCoALib (i.e.  ar -s ...)
#
# Revision 1.8  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.7  2008/06/30 17:13:49  abbott
# Changed progress messages (now directory names end in a slash)
#
# Revision 1.6  2007/05/14 16:35:07  bigatti
# -- removed numerical directory (for Dortmund library)
#
# Revision 1.5  2007/03/23 18:38:43  abbott
# Separated the "convert" functions (and CheckedCast) into their own files.
# Many consequential changes.  Also corrected conversion to doubles.
#
# Revision 1.4  2007/03/16 17:42:41  abbott
# Improved a compilation message.
#
# Revision 1.3  2007/03/09 21:01:27  abbott
# Corrected clean target (previsouly skipped APPDIRS).
#
# Revision 1.2  2007/03/09 20:58:36  abbott
# Final retouches: modified ar flags (again), more important change to main Makefile.
#
# Revision 1.1.1.1  2007/03/09 15:16:11  abbott
# Imported files
#
