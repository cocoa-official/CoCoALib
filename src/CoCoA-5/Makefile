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


CoCoAInterpreter:  VersionInfo.o  CompilationDate.o
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


VersionInfo.o:  VersionInfo.C  check-version-defines  version.mk
	@echo "Compiling VersionInfo.o (with special preprocessor flags)"
	$(COMPILE) -c $(COCOA5_VERSION_DEFINES) -o VersionInfo.o VersionInfo.C


.PHONY: cocoa5
cocoa5:  check-prerequisites-verbose  check-cocoalib
	@$(MAKE) CoCoAInterpreter

.PHONY:  check-cocoalib
check-cocoalib:
	@(cd ../..; $(MAKE) -s lib)


# This target checks that the 3 fields for the version are purely alphanumeric
check-version-defines:  check-version-defines.C  version.mk
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
check:  cocoa5
	@cd tests; $(MAKE) -s check


#------------------------------------------------------------------
# This part is relevant only if the QT-gui is present

.PHONY: Qt-gui
Qt-gui:  check-prerequisites-verbose  $(OBJS)  $(COCOA_LIB)  QCodeEdit/Makefile  C5Makefile
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
install:  cocoa5  doc  wordlist
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
	$(INSTALL_CMD) -m 755  cocoa5  CoCoAInterpreter "$$SUBDIR"; \
	/bin/mkdir -m 755 "$$SUBDIR/CoCoAManual"; \
	/bin/cp -r CoCoAManual/html  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644  CoCoAManual/tex/CoCoAManual.pdf  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644  CoCoAManual/CoCoAHelp.xml  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644  CoCoAManual/wordlist.txt  "$$SUBDIR/CoCoAManual"; \
	$(INSTALL_CMD) -m 644  ../../doc/CoCoATranslationTable.html  "$$SUBDIR/CoCoAManual"; \
	/bin/rm -f "$(COCOA5_INSTALL_DIR)/cocoa5"; \
	/bin/ln -s "$$SUBDIR/cocoa5" "$(COCOA5_INSTALL_DIR)"
	@echo "CoCoA-5 installed in $(COCOA5_INSTALL_DIR)"


AllExamples.o: $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: ManExamples
ManExamples:  cocoa5  OnlineHelp.o  PrintManExamples.o
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


wordlist.o:  $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: wordlist
wordlist:  $(COCOA5_WORDLIST)

$(COCOA5_WORDLIST):  OnlineHelp.C  wordlist.C  $(COCOA5_DOC_XML)
	$(MAKE)  -s  OnlineHelp.o  wordlist.o
	$(COMPILE5) OnlineHelp.o wordlist.o -o PrintWordlist $(LDLIBS)
	@/bin/rm -rf PrintWordlist.dSYM # Apple's compiler creates this dir
	@echo "creating wordlist..."; ./PrintWordlist | grep -v -F "OBSOLE" > $(COCOA5_WORDLIST)
	@/bin/rm PrintWordlist

.PHONY: doc
doc:  cocoa5
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
