# Makefile for CoCoALib/src/CoCoA-5/tests/ directory

COCOA_ROOT=../../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/CoCoA-5/tests/

TESTS=exsegv.cocoa5 \
      exbugs.cocoa5 \
      bug-EvalTwice.cocoa5 \
      ErrMesg.cocoa5 \
      test-output.cocoa5 \
      quine.cocoa5 \
      AnonFunc.cocoa5 \
      SourceAnna.cocoa5 \
      demo-GeMiTo2011.cocoa5 \
      demo-Osaka2015.cocoa5 \
      test-manual.cocoa5 \
      whatiscocoa.cocoa5 \
      lecture-HF1.cocoa5 \
      lecture-HF2.cocoa5 \
      lecture-HF3.cocoa5 \
      lecture-HF4.cocoa5 \
      tut-CoCoLA1.cocoa5 \
      tut-CoCoLA2.cocoa5 \
      tut-CoCoLA3.cocoa5 \
      tut-CoCoLA4.cocoa5 \
      RealRoots.cocoa5 \
      ExtLibNORMALIZ.cocoa5 \
      test-ApproxSolve.cocoa5 \
      test-ArrAndPosets.cocoa5 \
      test-FactorAlgExt.cocoa5 \
      test-GB.cocoa5 \
      test-HomomorphismOps.cocoa5 \
      test-implicit.cocoa5 \
      test-PrimaryDecomposition0.cocoa5 \
      test-RingElems.cocoa5 \
      test-radical.cocoa5 \
      test-radical0dim.cocoa5 \
      test-saturate.cocoa5 \
      test-syz.cocoa5 \
      test-toric.cocoa5 \
      primary.cocoa5 \
      tricky-references.cocoa5 \
      TutHokkaido2.cocoa5 \
      TutHokkaido3.cocoa5 \
      TutHokkaido4.cocoa5 \
      TutHokkaido5.cocoa5 \
      tagging.cocoa5


.PHONY: all
all: check

.PHONY: lib
lib: library

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) -s library)

# This target should be made only after the CoCoA library has been compiled;
# normally it would be called by the command "make RunTests" in the CoCoALib root
# directory.  Most of the work is done by the RunTests.sh script.
.PHONY: check
check: ../CoCoAInterpreter
	@./RunTests.sh $(TESTS)

.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f  ./*.found  ./*.cerr  ./*~  ./.\#*

.PHONY: veryclean
veryclean: clean-local
	@echo "Verycleaned CoCoALib/$(CWD)"



$(EXECS): $(COCOA_LIB)
