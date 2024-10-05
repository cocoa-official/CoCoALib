# Makefile for CoCoALib/src/CoCoA-5 directory.

CXXFLAGS_CUSTOM=-g

# installation directory  COCOA5_INSTALL_DIR  defined in autoconf.mk (set by configure script)

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
DEPEND_FILE=Makefile_dependencies

include version.mk
COCOA5_VERSION_DEFINES=-DCOCOA5_VER_MAJ="$(COCOA5_VER_MAJ)"  -DCOCOA5_VER_MIN="$(COCOA5_VER_MIN)"  -DCOCOA5_VER_MINMIN="$(COCOA5_VER_MINMIN)"
CWD=src/CoCoA-5/

COCOA5_DOC_DIR=CoCoAManual
COCOA5_DOC_XML=$(COCOA5_DOC_DIR)/CoCoAHelp.xml
COCOA5_WORDLIST=$(COCOA5_DOC_DIR)/wordlist.txt
COCOA5_DOC_AUX_DIR=$(COCOA5_DOC_DIR)/aux-files
COCOA5_DOC_HTML_DIR=$(COCOA5_DOC_DIR)/html
COCOA5_DOC_TEX_DIR=$(COCOA5_DOC_DIR)/tex

SRCS=AST.C  Lexer.C  Main.C  Interpreter.C  LineProviders.C  Parser.C \
     CoCoALibSupplement.C \
     BuiltInFunctions.C \
     BuiltInFunctions-CoCoALib.C \
     BuiltInFunctionsVarArgs-CoCoALib.C  BuiltInOneLiners-CoCoALib.C \
     BuiltInFunctions-Frobby.C \
     BuiltInFunctions-GFan.C \
     BuiltInFunctions-GSL.C \
     BuiltInFunctions-MathSAT.C \
     BuiltInFunctions-Normaliz.C \
     globals.C \
     OnlineHelp.C \
     VersionInfo.C Banner.C
#     CompilationDate.C


# Default target
.PHONY: all
all:  check  doc  wordlist
	@if [ $(BUILD_QT_GUI) = "yes" ]; \
	 then \
	   $(MAKE) Qt-gui; \
	 fi

OBJS=$(SRCS:.C=.o)

$(OBJS): Makefile  $(COCOA_ROOT)/configuration/autoconf.mk

CompilationDate.o: $(OBJS)  $(COCOA_LIB)  CompilationDate.C
	@echo "Compiling CompilationDate.o"
	$(COMPILE5) -c -o CompilationDate.o CompilationDate.C

COMPILE5=$(COMPILE) $(CXXFLAGS_CUSTOM)  $(COCOA5_CXX_DEFINES)


# Rule for compiling C++ code in *.C files into *.o object files
%.o: %.C
	@echo "Compiling `basename $@`"
	$(COMPILE5) -c -o $@ $<


CoCoAInterpreter: VersionInfo.o  CompilationDate.o
	@echo "CoCoAInterpreter: linking everything together"
	$(COMPILE5) $(OBJS) CompilationDate.o -o CoCoAInterpreter  $(LDLIBS)  $(COCOA5_LDLIBS)  -lpthread
	@#Next lines remove the *.dSYM directory which the Apple compiler creates
	@#but we remove it only if the compilation flag "-g" is absent.
	@AppleDir=CoCoAInterpreter.dSYM ; \
	echo " " $(CXXFLAGS) " " $(CXXFLAGS_CUSTOM) " " | grep -F " -g " >/dev/null; \
	if [ $$? -eq 1 -a -d "$$AppleDir" ] ; \
	then \
	  /bin/rm -rf "$$AppleDir"; \
	fi


VersionInfo.o: VersionInfo.C  check-version-defines  version.mk
	@echo "Compiling VersionInfo.o (with special preprocessor flags)"
	$(COMPILE) -c $(COCOA5_VERSION_DEFINES) -o VersionInfo.o VersionInfo.C


.PHONY: cocoa5
cocoa5: check-prerequisites-verbose  check-cocoalib
	@$(MAKE) CoCoAInterpreter

.PHONY: check-cocoalib
check-cocoalib:
	@(cd ../..; $(MAKE) -s lib)


# This target checks that the 3 fields for the version are purely alphanumeric
check-version-defines: check-version-defines.C  version.mk
	@$(COMPILE)  $(COCOA5_VERSION_DEFINES)  -o check-version-defines  check-version-defines.C
	@./check-version-defines; if [ $$? -ne 0 ]; \
	then \
	  echo "------------------------------------------------------------"; \
	  echo "ERROR: version fields not alphanumeric (see file version.mk)"; \
	  echo "------------------------------------------------------------"; \
	  exit 1; \
	fi

### 2022-11-11: this target no longer needed?
### .PHONY: check-prerequisites
### check-prerequisites:
###	@test "$(HAVE_BOOST)" = "yes"

.PHONY: check-prerequisites-verbose
check-prerequisites-verbose:
	@if [ "$(HAVE_BOOST)" = "no" ]; \
	 then \
	   echo "ERROR: Not building CoCoAInterpreter (& GUI):"; \
	   echo "ERROR: configuration option --only-cocoalib specified;"; \
	   echo "ERROR: Reconfigure *without* option --only-cocoalib (also requires BOOST libs)"; \
	   exit 1; \
	 fi



.PHONY: check
check: cocoa5
	@cd tests; $(MAKE) -s check


#------------------------------------------------------------------
# This part is relevant only if the QT-gui is present

.PHONY: Qt-gui
Qt-gui: check-prerequisites-verbose  $(OBJS)  $(COCOA_LIB)  QCodeEdit/Makefile  C5Makefile
	@echo "-----------------------------------"
	@echo ">>>>  Making the CoCoA Qt GUI  <<<<"
	@echo "-----------------------------------"
	@echo "[[CoCoA Qt-GUI: 1/3]]  Building QCodeEdit";  cd QCodeEdit;  CXX=$(CXX)  $(MAKE)  -s
	@echo "[[CoCoA Qt-GUI: 2/3]]  Building GUI";        CXX=$(CXX)  $(MAKE) -s -f C5Makefile
	@echo "[[CoCoA Qt-GUI: 3/3]]  Finishing GUI";       ./make-gui-finish.sh
	@echo
	@echo "---------------------------------------------------------"
	@echo " (-:  Amazing!!  CoCoA-5 Qt GUI successfully built!  :-) "
	@echo "---------------------------------------------------------"

QCodeEdit/Makefile:  $(COCOA_LIB)  QCodeEdit/QCodeEdit.pro.in  make-qcodeeditmakefile.sh
	if [ -f QCodeEdit/Makefile ]; then cd QCodeEdit; $(MAKE) -s clean; fi
	QMAKE=$(QMAKE)  ./make-qcodeeditmakefile.sh

C5Makefile:  $(COCOA_LIB)  C5.pro.in  make-c5makefile.sh
	if [ -f C5Makefile ]; then $(MAKE) -s -f C5Makefile distclean; fi
	QMAKE=$(QMAKE)  ./make-c5makefile.sh  $(CXXFLAGS_DEFINES)  $(CXXFLAGS)  $(COCOA5_VERSION_DEFINES)  $(LDLIBS)  $(COCOA5_LDLIBS)

# End of Qt-gui specific part
#------------------------------------------------------------------


.PHONY: install
install: cocoa5 doc wordlist
	@echo "===================================================================="
	@echo ">>> WARNING: CoCoA-5 installation procedure is still PRELIMINARY <<<"
	@echo "===================================================================="
	@echo
	@echo "Proceeding with installation after 5 seconds..."
	@echo "Target directory: $(COCOA5_INSTALL_DIR)"
	@echo
	@sleep 5
	@if [ \! -d "$(COCOA5_INSTALL_DIR)" ]; then echo "***** ERROR: $(COCOA5_INSTALL_DIR) does not exist or is not a directory *****"; exit 1; fi
	@if [ \! -w "$(COCOA5_INSTALL_DIR)" ]; then echo "***** ERROR: $(COCOA5_INSTALL_DIR) is not a writable directory *****"; echo "***** >>>>  Consider using \"sudo\" command  <<<<"; exit 1; fi
	@SUBDIR="$(COCOA5_INSTALL_DIR)/cocoa-$(COCOA5_VER_MAJ).$(COCOA5_VER_MIN).$(COCOA5_VER_MINMIN)"; \
	if [ -e "$$SUBDIR" ]; then /bin/chmod -R u+w "$$SUBDIR"; /bin/rm -rf "$$SUBDIR"; fi; \
	/bin/mkdir -m 755 "$$SUBDIR"; \
	/bin/cp -r packages emacs "$$SUBDIR"; \
	/bin/chmod 755 "$$SUBDIR/packages" "$$SUBDIR/emacs"; \
	$(INSTALL_CMD) -m 755 cocoa5 CoCoAInterpreter "$$SUBDIR"; \
	/bin/mkdir -m 755 "$$SUBDIR/CoCoAManual"; \
	/bin/cp -r CoCoAManual/html  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644 CoCoAManual/tex/CoCoAManual.pdf  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644 CoCoAManual/CoCoAHelp.xml  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644 CoCoAManual/wordlist.txt  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644 ../../doc/CoCoATranslationTable.html  "$$SUBDIR/CoCoAManual"; \
	/bin/rm -f "$(COCOA5_INSTALL_DIR)/cocoa5"; \
	/bin/ln -s "$$SUBDIR/cocoa5" "$(COCOA5_INSTALL_DIR)"
	@echo "CoCoA-5 installed in $(COCOA5_INSTALL_DIR)"


AllExamples.o: $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: ManExamples
ManExamples: cocoa5 OnlineHelp.o  PrintManExamples.o
	$(COMPILE5) OnlineHelp.o  PrintManExamples.o  -o PrintManExamples $(LDLIBS)
	/bin/rm -rf PrintManExamples.dSYM # Apple's compiler creates this dir
	@echo "creating ManExamples-in.cocoa5 ...";
	 ./PrintManExamples > $(COCOA_ROOT)/src/CoCoA-5/ManExamples-in.cocoa5
	@echo "...done.  Running examples...";
	./CoCoAInterpreter < ManExamples-in.cocoa5 > ManExamples-out;
	@echo "...done."
	@echo "===================================================="
	@echo ">>> NOW search for  \"ERROR:\"  in ManExamples-out <<<"
	@echo ">>> readline may truncate echo lines if too long <<<"
	@echo "===================================================="


wordlist.o: $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: wordlist
wordlist: $(COCOA5_WORDLIST)

$(COCOA5_WORDLIST):  OnlineHelp.C  wordlist.C  $(COCOA5_DOC_XML)
	$(MAKE)  -s  OnlineHelp.o  wordlist.o
	$(COMPILE5) OnlineHelp.o wordlist.o -o PrintWordlist $(LDLIBS)
	@/bin/rm -rf PrintWordlist.dSYM # Apple's compiler creates this dir
	@echo "creating wordlist..."; ./PrintWordlist | grep -v -F "OBSOLE" > $(COCOA5_WORDLIST)
	@/bin/rm PrintWordlist

.PHONY: doc
doc:
	@(cd $(COCOA5_DOC_DIR); $(MAKE))


# OptCoCoAInterpreter: $(SRCS)
# 	$(CXX) -o OptCoCoAInterpreter $(C_FLAGS) $(OPT_FLAGS) $(INCLUDE5) *.C $(LDLIBS)

.PHONY: lib
lib: library

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)


.PHONY: clean  clean-local  clean-subdirs
clean: clean-local  clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f $(OBJS)  CompilationDate.o  wordlist.o  "$(DEPEND_FILE).old"  "$(DEPEND_FILE).new"  core  gmon.out  ./*~  ./.*~  ./.\#*
	@/bin/rm -rf  C5.app  ./*.dSYM
	@/bin/rm -f  C5 .qmake.stash
	@/bin/rm -f  C5Makefile  C5.pro  moc_predefs.h
	@/bin/rm -f  moc_C5.cpp  qrc_C5.cpp  ui_Console.h  ui_Debugger.h  ui_MainWindow.h  ui_SourceEditor.h
	@/bin/rm -f  CoCoAInterpreter  OptCoCoAInterpreter  check-version-defines
	@/bin/rm -f  PrintWordlist  PrintManExamples  PrintManExamples.o
	@/bin/rm -f  AllExamples.o  ManExamples-in.cocoa5  ManExamples-out  my-test  my-file  PLOT-*  CoCoAPlot

clean-subdirs:
	cd $(COCOA5_DOC_DIR); /bin/rm -f  ./*~  ./.\#*
	cd emacs; /bin/rm -f ./*~  ./.\#*
	cd packages; /bin/rm -f ./*~  ./.\#*
	cd tests; /bin/rm -f ./*~  ./.\#*  ./*.found  ./*.cerr
	cd QCodeEdit/; if [ -f "Makefile" ]; then $(MAKE) -s -o  QCodeEdit.pro  distclean; fi
	/bin/rm -f QCodeEdit/QCodeEdit.pro  QCodeEdit/Makefile
	/bin/rm -rf QCodeEdit/.build/  GUI-debug/  GUI-release/


.PHONY: veryclean  veryclean-subdirs
veryclean:  clean-local  veryclean-subdirs
	@/bin/rm -f  .qmake.stash  QCodeEdit/libqcodeedit.a  QCodeEdit/Makefile
	@true >| "$(DEPEND_FILE)" # Empties $(DEPEND_FILE)
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-subdirs: clean-subdirs
	/bin/rm -f $(COCOA5_WORDLIST)
	/bin/rm -rf $(COCOA5_DOC_TEX_DIR)
	/bin/rm -rf $(COCOA5_DOC_HTML_DIR)
	if [ -d release-files ] ; \
        then \
	  cd release-files; \
	  /bin/rm -rf  MacCoCoA-5.*.app-*  CoCoA-5.*-*  cocoa-5.*-*; \
	fi


# This target checks whether we should really build DEPEND_FILE
.PHONY: dependencies
dependencies:
	@if [ \! -s "$(DEPEND_FILE)" -o \! "$(DEPEND_FILE)" -nt "$(COCOA_HDR)" ] ; \
	then \
	  $(MAKE) create-dependencies; \
	fi

# This is the target which really does rebuild the DEPEND_FILE
.PHONY: create-dependencies
create-dependencies:
	@/bin/rm -rf $(DEPEND_FILE).new
	@echo "Rebuilding dependencies in <CoCoALib>/$(CWD)"
	@echo "##################################################################"   >> $(DEPEND_FILE).new
	@echo "#        ---  CoCoALib version $(COCOALIB_VERSION)  ---"              >> $(DEPEND_FILE).new
	@echo "# Dependencies for directory $(CWD)"                                  >> $(DEPEND_FILE).new
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"   >> $(DEPEND_FILE).new
	@echo "# run \`make dependencies' in the CoCoA-5 directory."                 >> $(DEPEND_FILE).new
	@echo "##################################################################"   >> $(DEPEND_FILE).new
	@echo                                                                        >> $(DEPEND_FILE).new
	@if [ $(HAVE_BOOST) = "no" ] ; \
	then \
	  echo "***************************************"; \
	  echo "*** configured with --only-cocoalib ***   ==> not rebuilding CoCoA-5 dependencies"; \
	  echo "***************************************"; \
	  exit 0; \
	fi; \
	for srcfile in $(SRCS); \
	do \
	  objfile=`basename $$srcfile .C`.o; \
	  $(CXX) $(CXXFLAGS) $(CXXFLAGS_DEFINES) $(INCLUDE_PATHS) -MM "$$srcfile" -MQ "$$objfile" >> $(DEPEND_FILE).new; \
	  echo >> $(DEPEND_FILE).new; \
	done
	@if [ -s $(DEPEND_FILE) ] ; \
	then \
	  /bin/mv -f $(DEPEND_FILE) $(DEPEND_FILE).old ; \
	fi
	@/bin/mv $(DEPEND_FILE).new  $(DEPEND_FILE)


$(DEPEND_FILE):
	@touch "$(DEPEND_FILE)"

include $(DEPEND_FILE)


$(COCOA_ROOT)/configuration/autoconf.mk:
	@echo
	@echo "======================================================================="
	@echo ">>>  ERROR: you must first run configure script in directory $(COCOA_ROOT)  <<<"
	@echo "======================================================================="
	@echo
	@exit 1



# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/Makefile,v 1.166 2024/03/22 20:57:00 abbott Exp $
# $Log: Makefile,v $
# Revision 1.166  2024/03/22 20:57:00  abbott
# Summary: Changes to permit compilationon MinGW (redmine 1661)
#
# Revision 1.165  2024/03/21 20:10:34  abbott
# Summary: make lib is now "silent"
#
# Revision 1.164  2024/03/18 08:50:30  abbott
# Summary: Fixed re-linking problem (redmine 1795)
#
# Revision 1.163  2024/02/09 09:17:34  abbott
# Summary: Tried to insert a check that cocoalib is compiled
#
# Revision 1.162  2024/01/30 20:48:21  abbott
# Summary: Better error mesg if autoconf.mk does not exist
#
# Revision 1.161  2022/11/30 20:39:37  abbott
# Summary: Minor cleaning, improved layout
#
# Revision 1.160  2022/11/16 19:49:03  abbott
# Summary: Replaced fgrep (obsol) by grep; revised "no boost" err mesgs
#
# Revision 1.159  2022/11/04 15:26:14  abbott
# Summary: Several changes: "make" to $(MAKE); qmake to $(QMAKE); aldo added a BOOST_LIB_DIR arg for linker
#
# Revision 1.158  2022/04/10 16:19:32  abbott
# Summary: Corrected typo; removed some cruft (which had been put into CoCoAManual/Makefile)
#
# Revision 1.157  2022/03/17 16:17:18  abbott
# Summary: Aesthetic improvement
#
# Revision 1.156  2022/03/17 14:51:38  abbott
# Summary: Removed several superfluous include directives
#
# Revision 1.155  2022/03/07 16:56:28  bigatti
# Summary: added help line for ManExamples (about readline)
#
# Revision 1.154  2022/03/02 14:40:40  abbott
# Summary: Changed clean-local to remove new dirs (after changing GUI dir names)
#
# Revision 1.153  2022/02/25 10:40:00  abbott
# Summary: Revised so that Qt-gui builds & cleans properly (redmine 1393)
#
# Revision 1.152  2022/02/05 20:02:04  abbott
# Summary: Changed order of all target
#
# Revision 1.151  2022/02/04 20:07:50  abbott
# Summary: Changed order of target all: removed cocoa5 (implied by check)
#
# Revision 1.150  2021/04/01 20:03:52  abbott
# Summary: Now exclude obsolete wprds from wordlist (via fgrep -v)
#
# Revision 1.149  2020/12/12 20:32:55  abbott
# Summary: Now allow version to contain alphanumeric components
#
# Revision 1.148  2020/10/27 09:44:00  abbott
# Summary: Changed name of make var INCLUDE to INCLUDE_PATHS
#
# Revision 1.147  2020/03/07 14:45:05  abbott
# Summary: clean now removes PrintManExamples.o
#
# Revision 1.146  2020/03/06 20:35:15  abbott
# Summary: Improved behaviour of building Makefile_dependencies (even when that file does not exist)
#
# Revision 1.145  2020/02/12 15:45:24  abbott
# Summary: Added some helpful spaces
#
# Revision 1.144  2020/02/03 16:45:02  bigatti
# -- renamed AllExamples.C --> PrintManExamples.C
#
# Revision 1.143  2020/01/27 19:58:18  abbott
# Summary: Moved defn of DEPEND_FILE closer to top; removed some cruft
#
# Revision 1.142  2020/01/14 21:03:24  abbott
# Summary: Commented out dependency of Makefile_depend on configuration/version
#
# Revision 1.141  2020/01/09 18:36:19  abbott
# Summary: Removed non-existent PHONY target
#
# Revision 1.140  2019/12/18 08:21:52  abbott
# Summary: Renamed target "gui" to "Qt-gui"; now depends on flag BUILD_QT_GUI (instead of HAVE_QMAKE)
#
# Revision 1.139  2019/10/29 11:37:12  abbott
# Summary: Changed format of "skipping" message (more readable now)
#
# Revision 1.138  2019/10/09 12:56:52  abbott
# Summary: Remove also my-test and my-file (might be created by make ManExamples)
#
# Revision 1.137  2019/10/08 20:31:34  abbott
# Summary: Made clean-local remove also core and gmon.out
#
# Revision 1.136  2019/09/25 16:08:34  bigatti
# -- added BuiltInFunctionsVarArgs-CoCoALib
#
# Revision 1.135  2019/03/04 13:15:53  abbott
# Summary: Added new files globals.H/globals.C
#
# Revision 1.134  2018/08/06 16:23:19  abbott
# Summary: Improved install target
#
# Revision 1.133  2017/11/15 17:00:37  abbott
# Summary: Changed name of COCOA5_LDLIBS into LDLIBS_COCOA5
#
# Revision 1.132  2017/10/17 10:20:34  abbott
# Summary: Corrected badly placed space
#
# Revision 1.131  2017/10/04 10:30:24  abbott
# Summary: Removed useless rule; conseq change after separating LDLIBS for CoCoA-5 from those for CoCoALib
#
# Revision 1.130  2017/07/21 13:27:46  abbott
# Summary: Building of Makefile_dependencies should now be "clean" even when interrupted
#
# Revision 1.129  2017/07/14 09:32:04  bigatti
# -- added MathSAT
#
# Revision 1.128  2017/04/28 20:13:41  abbott
# Summary: Made for obvious the final mesg for ManExamples
#
# Revision 1.127  2017/04/28 15:48:19  abbott
# Summary: Also clean target now removes AllExamples.o
#
# Revision 1.126  2017/04/28 15:45:58  abbott
# Summary: Clean now removes ManExamples-out too
#
# Revision 1.125  2017/04/28 14:55:57  bigatti
# -- improved reminder for ManExamples
#
# Revision 1.124  2017/04/28 14:28:47  bigatti
# -- reminder for checkin ManExamples-out
#
# Revision 1.123  2017/03/01 17:26:18  abbott
# Summary: Removed bloody make veryclean-subdirs from QCodeEdit/Makefile rule
#
# Revision 1.122  2017/02/22 12:36:18  abbott
# Summary: Changed dependencies for wordlist
#
# Revision 1.121  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.120  2017/02/02 12:52:25  abbott
# Summary: Minor correction to installation script
#
# Revision 1.119  2016/12/13 17:16:22  abbott
# Summary: Now prints message after building the GUI
#
# Revision 1.118  2016/12/08 10:36:29  abbott
# Summary: Improved "install" target
#
# Revision 1.117  2016/12/01 10:21:30  abbott
# Summary: "make install" now creates symbolic link too
#
# Revision 1.116  2016/11/30 20:12:13  abbott
# Summary: Minor improvement to cocoa5 "install" target; dir is now specfied by configure script, and put in autoconf.mk
#
# Revision 1.115  2016/11/30 20:00:08  abbott
# Summary: Added first version of "install" target
#
# Revision 1.114  2016/10/12 13:49:11  bigatti
# -- fix for make ManExamples, and output redirected to ManExamples-out
#
# Revision 1.113  2016/05/18 12:28:08  abbott
# Summary: Added new target "lib": synonym for "library"
#
# Revision 1.112  2016/03/21 15:10:15  bigatti
# -- automatic version numbering for html doc
#
# Revision 1.111  2016/03/21 14:44:03  bigatti
# added automatic version numbering to TeX documentation
#
# Revision 1.110  2016/03/19 20:08:12  abbott
# Summary: Revised to work with new release-files/version
#
# Revision 1.109  2016/03/18 12:12:48  bigatti
# -- Renamed BuiltInOneLiners-CoCoALib
#
# Revision 1.108  2016/03/09 10:18:08  abbott
# Summary: Minor aesthetic change
#
# Revision 1.107  2016/02/18 11:29:11  abbott
# Summary: Now uses COCOA5_VERSION_DEFINES; also passes value to make-c5makefile
#
# Revision 1.106  2016/02/17 11:12:23  bigatti
# -- added version minmin
#
# Revision 1.105  2016/01/06 12:16:36  abbott
# Summary: Minor change to a progress message
#
# Revision 1.104  2015/12/02 17:07:50  abbott
# Summary: Corrected target clean-subdirs (remove *.found and *.cerr in tests subdir)
#
# Revision 1.103  2015/10/08 13:09:13  bigatti
# -- BuiltInFunctions.C split into BuiltInFunctions + BuiltInFunctions-CoCoALib
#
# Revision 1.102  2015/09/03 09:33:25  bigatti
# -- added BuiltInFunctions-GFan
#
# Revision 1.101  2015/09/02 14:47:54  bigatti
# -- make veryclean-subdirs also deletes old copies of release dirs
#   (but not the latest...)
#
# Revision 1.100  2015/04/16 20:23:41  abbott
# Summary: target clean-local now delete CompilationDate.o too
# Author: JAA
#
# Revision 1.99  2014/12/01 10:35:14  bigatti
# -- added "worlist" to "doc" target
#
# Revision 1.98  2014/10/30 13:23:42  abbott
# Summary: CompilationDate.o now depends on COCOA_LIB
# Author: JAA
#
# Revision 1.97  2014/09/05 16:12:24  abbott
# Summary: Missed a couple of changes CoCoAInterpreter --> cocoa5; now done
# Author: JAA
#
# Revision 1.96  2014/09/05 16:10:29  abbott
# Summary: Now rebuilds CompilationDate only when necessary; CoCoAInterpreter cleaner; new target cocoa5
# Author: JAA
#
# Revision 1.95  2014/09/01 12:09:27  abbott
# Summary: Strengthened clean-subirs; removed redundancy from veryclean-subdirs
# Author: JAA
#
# Revision 1.94  2014/08/01 11:54:08  bigatti
# -- fixed veryclean
#
# Revision 1.93  2014/08/01 11:51:43  bigatti
# == changed ManExamples (and veryclean)
#
# Revision 1.92  2014/08/01 11:32:02  bigatti
# -- changed ManExamples
#
# Revision 1.91  2014/07/30 16:56:43  bigatti
# ++ added "make doc"
#
# Revision 1.90  2014/07/28 14:55:07  abbott
# Summary: Improved targets (very)clean-local; new scheme for building dependencies
# Author: JAA
#
# Revision 1.89  2014/07/22 11:15:32  bigatti
# -- changed output file for target ManExamples
#
# Revision 1.88  2014/06/28 14:55:26  abbott
# Summary: Revised in line with new configuration scripts
# Author: JAA
#
# Revision 1.87  2014/06/17 10:17:57  abbott
# Summary: depend_file now depends on library.H (instead of version)
# Author: JAA
#
# Revision 1.86  2014/06/13 12:06:44  abbott
# Summary: Corrected name of "CoCoAManual/aux" subdir
# Author: JAA
#
# Revision 1.85  2014/05/14 13:16:21  bigatti
# -- ManExamples nowcreates ManExamples-out
#
# Revision 1.84  2014/05/14 07:37:29  bigatti
# -- commented out line for manual about old gui
#
# Revision 1.83  2014/04/22 12:47:47  abbott
# Summary: Added capability to increase max stack size when starting CoCoAInterpreter
# Author: JAA
#
# Revision 1.82  2014/04/17 12:22:04  abbott
# Summary: Changed clean-local, veryclean-local; new approach for rebuilding dependencies
# Author: JAA
#
# Revision 1.81  2014/04/04 16:41:12  bigatti
# Summary: makefile gui: gmpxx must be last!
#
# Revision 1.80  2014/03/21 12:01:39  abbott
# Summary: Corrected (& simplified) defn of LIBS5
# Author: JAA
#
# Revision 1.79  2014/03/14 11:03:31  abbott
# Summary: clean and veryclean targets now delete *~ and .*~ files too
# Author: JAA
#
# Revision 1.78  2014/03/07 14:17:57  abbott
# Summary: Corrected dependencies for C5Makefile and QCodeEdit/Makefile
# Author: JAA
#
# Revision 1.77  2014/03/05 11:15:23  bigatti
# -- undoing wrong commit
#
# Revision 1.75  2014/01/29 18:55:34  abbott
# Summary: Improved clean-subdirs target; wordlist.txt target now deletes PrintWordList executable.
# Author: JAA
#
# Revision 1.74  2014/01/29 18:00:12  abbott
# Summary: Target veryclean-local now deletes wordlist (generated file)
# Author: JAA
#
# Revision 1.73  2014/01/29 15:46:25  abbott
# Summary: Moved some rm command from veryclean-local to clean-local
# Author: JAA
#
# Revision 1.72  2014/01/28 09:49:23  abbott
# Corrected veryclean-local & veryclean-subdirs targets.
#
# Revision 1.71  2013/06/26 17:37:26  bigatti
# -- added BuiltInFunctions-Frobby.C
#
# Revision 1.70  2013/05/14 14:24:47  abbott
# Revised to be compatible with new BOOST libs script.
#
# Revision 1.69  2013/05/08 09:43:14  abbott
# Added dependency of all object files on the Makefile itself.
#
# Revision 1.68  2013/03/15 17:47:46  abbott
# Further change to make sure Qt dependent code is rebuilt when autoconf.mk changes.
#
# Revision 1.67  2013/03/15 11:03:36  abbott
# Added extra line to force rebuilding of QCodeEdit stuff when autoconf.mk changes.
#
# Revision 1.66  2013/03/12 15:16:57  abbott
# Attempt to fix a "wrong architecture" error message.
#
# Revision 1.65  2013/03/07 17:49:34  abbott
# Added (first impl of) new command SourceRegion.
# Had to hack uglily readAndExecute (in Interpreter.C).
#
# Revision 1.64  2012/10/24 14:05:59  abbott
# Made veryclean-local remove more files.
#
# Revision 1.63  2012/10/02 12:12:40  abbott
# Simplified defn of COMPILE5.
# Removed defn of INCLUDE5.
# Removed commented out code referring to ZIP library.
#
# Revision 1.62  2012/10/01 13:59:30  abbott
# Updated command for generating dependency file.
#
# Revision 1.61  2012/09/27 15:29:09  abbott
# Simplified and updated the rule for building C5Makefile.
#
# Revision 1.60  2012/08/06 09:45:41  abbott
# Removed an unwanted space in a printed message.
#
# Revision 1.59  2012/08/06 09:20:07  abbott
# Improved layout of message saying we're not building GUI.
# Added wordlist.o to files to remove when cleaning.
#
# Revision 1.58  2012/08/06 09:14:53  abbott
# Added a check for HAVE_QMAKE so we give the user a helpful message instead of
# an error when trying to builf GUI but QMAKE is absent.
#
# Revision 1.57  2012/08/05 12:26:03  abbott
# Added check-prerequisites and check-prerequisites-verbose targets
# (replaced the old targets boost-hdr-check and boost-lib-check)
#
# Revision 1.56  2012/08/02 16:35:35  abbott
# Added check-prerequisites target, so main Makefile can give a helpful mesg if there's a problem.
# Added a check for BOOST when rebuilding dependencies -- if absent, the dependency file is left
# (largely) empty.
#
# Revision 1.55  2012/07/19 17:03:23  bigatti
# -- passing GMPXX_LIB as last argument to make-c5makefile
#   (for linker problem with Normaliz)
#
# Revision 1.54  2012/07/16 10:20:24  bigatti
# -- added GSL and Normaliz arguments to "make-c5makefile" call
#
# Revision 1.53  2012/07/11 09:05:06  bigatti
# -- added BuiltInFunctions-Normaliz and BuiltInFunctions-GSL
#
# Revision 1.52  2012/07/05 14:24:25  bigatti
# -- added BuiltInOneLiners.C
#
# Revision 1.51  2012/07/03 12:43:20  abbott
# Added specific path for PrintWordlist
#
# Revision 1.50  2012/06/26 08:52:35  bigatti
# -- make target wordlist less verbose
#
# Revision 1.49  2012/06/26 06:52:49  bigatti
# -- removed experimental COCOA_5_OBJ_DIR
#
# Revision 1.48  2012/06/26 06:44:48  bigatti
# -- added texdoc target
# -- fixed CompilationDate proper compilation also for GUI
# -- added wordlist generation
#
# Revision 1.47  2012/06/18 10:04:22  bigatti
# -- added VersionInfo.C
#
# Revision 1.46  2012/06/07 16:07:05  bigatti
# -- added make for wordlist and htmldoc
#
# Revision 1.45  2012/04/11 09:55:59  abbott
# Makefile now performs boost-lib-check for CoCoAInterpreter
# (previsouly was only boost-hdr-check).
#
# Revision 1.44  2012/04/10 18:55:55  abbott
# Several changes to accomodate using BOOST filesystem fns in Main.C to check
# whether the default packages directory exists.  Needed to modify the script
# which generates C5Makefile as the extra BOOST libs must be passed as params.
#
# Revision 1.43  2012/03/16 17:04:53  abbott
# Corrected clean-subdirs target; it now leaves subdir Debug/ as empty but existing.
#
# Revision 1.42  2012/03/16 16:22:02  abbott
# Improved clean & clean-subdirs targets.
#
# Revision 1.41  2011/10/05 10:00:45  bigatti
# -- fixed make gui for non-standard installation of boost
#
# Revision 1.40  2011/10/03 16:07:57  abbott
# Corrected two minor bugs.
#
# Revision 1.39  2011/09/30 14:29:12  abbott
# Major revision to gui target recipe.
# Removed cruft from start of file.
#
# Revision 1.38  2011/09/30 13:02:33  bigatti
# -- modified target all: no longer fails if boost not present (just
#    prints a big warning)
# -- removed now useless configuration for boost
#
# Revision 1.37  2011/09/28 12:44:52  abbott
# Moved default target up -- correcting a bug I introduced :-(
#
# Revision 1.36  2011/09/21 14:48:11  abbott
# OBJS now depend on autoconf.mk ==> rebuild them if the compilation flags change.
#
# Revision 1.35  2011/09/12 16:49:36  bigatti
# -- added GMP_INC_DIR as 4th arg for make-gui.sh
#
# Revision 1.34  2011/08/17 12:00:48  abbott
# Changed the following targets clean-subdirs, veryclean, veryclean-local
# (cleaning works better in QCodeEdit directory).
#
# Revision 1.33  2011/07/21 17:34:35  abbott
# Consequential change for modified make-gui.sh script.
#
# Revision 1.32  2011/07/20 13:55:06  bigatti
# -- added some "@"
#
# Revision 1.31  2011/07/20 13:02:32  abbott
# Removed the word "version" from logging message when rebuilding dependencies.
#
# Revision 1.30  2011/07/20 09:55:51  abbott
# Corrected infinite recursion bug (when Makefile_dependencies doesn't exist).
#
# Revision 1.29  2011/07/20 09:04:19  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.28  2011/07/18 09:03:47  abbott
# Removed some pointless (& harmless?) double quotes.
# Added missing ./ before call to make-gui.sh
#
# Revision 1.27  2011/07/16 19:53:43  abbott
# New files and change to Makefile to allow C5 GUI to built automatically
# (as if it were a normal "make" target)
#
# Revision 1.26  2011/07/15 16:56:46  abbott
# Now gets BOOST location from autoconf.mk (i.e. as set by configure script).
#
# Revision 1.25  2011/05/18 15:22:42  bigatti
# -- fixed dependencies on cocoalib
#
# Revision 1.24  2011/05/16 19:02:20  lagorio
# Workaround per gcc bug
#
# Revision 1.23  2011/05/13 17:02:05  abbott
# Moved definition of OBJS; added two .PHONY declarations.
#
# Revision 1.22  2011/05/10 16:45:59  bigatti
# -- added  "make check"
#
# Revision 1.21  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.20  2011/02/15 09:39:34  bigatti
# -- added ManExamples target
#
# Revision 1.19  2011/01/27 15:52:10  bigatti
# -- added -DNDEBUG
#
# Revision 1.18  2011/01/19 16:36:58  bigatti
# -- removed dependency for BuiltInFunctions.o (I didn't do what I meant)
#
# Revision 1.17  2011/01/18 14:44:48  bigatti
# -- added #include CoCoALibSupplement.H
# -- added dependency of BuiltInFunctions from cocoalib.a
#
# Revision 1.16  2010/11/22 17:36:25  bigatti
# -- added "make library"
#
# Revision 1.15  2010/10/11 14:19:27  bigatti
# -- fixed boost path for bigatti
#
# Revision 1.14  2010/10/08 22:05:35  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.13  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" by calls to "$(MAKE)".
#
# Revision 1.12  2010/08/31 14:57:22  bigatti
# -- added OnlineHelp.C to sources
#
# Revision 1.11  2010/08/31 13:14:58  bigatti
# -- commented out compilation of zip library
#
# Revision 1.10  2010/08/23 20:19:22  lagorio
# *** empty log message ***
#
# Revision 1.9  2010/08/02 16:12:03  bigatti
# -- added code for zip support
#
# Revision 1.8  2010/05/07 14:19:01  bigatti
# -- added BuiltInFunctions.C
#
# Revision 1.7  2010/04/22 13:40:14  abbott
# Now inherits most compilation flags from autoconf.mk
# Make clean works properly.
# Object files now depend one the library.
#
# Revision 1.6  2010/04/22 09:30:17  bigatti
# -- added CompilationDate.C and forced its compilation for
#    CoCoAInterpreter target
#
# Revision 1.5  2010/04/21 16:08:06  bigatti
# -- more clever now: recompile only what is needed and with dependencies
#
# Revision 1.4  2010/04/13 15:24:44  bigatti
# -- updated after new "compile" file (with true interpreter starting!)
#
# Revision 1.3  2009/09/13 11:12:06  lagorio
# *** empty log message ***
#
# Revision 1.2  2009/09/13 08:41:13  lagorio
# *** empty log message ***
#
# Revision 1.1  2009/09/03 14:35:44  bigatti
# -- first import
#
