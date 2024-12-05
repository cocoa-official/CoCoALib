# Makefile for CoCoALib/examples directory

COCOA_ROOT=..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=examples/

SRCS=ex-empty.C \
      ex-c++-basic.C  ex-c++-integers.C  ex-c++-arith.C  ex-c++-loop-for1.C  ex-c++-loop-for2.C  ex-c++-loop-while.C \
      ex-c++-bool.C  ex-c++-fn-defn.C  ex-c++-vector1.C  ex-c++-vector2.C  ex-c++-class.C \
      ex-00-intro.C \
      ex-AlexanderDual.C \
      ex-ApproxPts1.C \
      ex-BigInt1.C  ex-BigInt2.C  ex-BigInt3.C  ex-BigRat1.C  ex-BigRatInterval1.C \
      ex-bool3.C \
      ex-BuildInfo.C \
      ex-convert1.C  ex-CoprimeFactorBasis1.C  ex-CpuTimeLimit1.C  ex-CpuTimeLimit2.C \
      ex-cyclotomic.C  ex-combinatorics1.C \
      ex-DivMask1.C  ex-DivMask2.C \
      ex-DynamicBitset1.C  ex-DynamicBitset2.C \
      ex-ExternalLibs1.C \
      ex-factor1.C  ex-FloatApprox1.C  ex-frobby1.C \
      ex-geobucket1.C \
      ex-GFan1.C \
      ex-GMPAllocator1.C  ex-GMPAllocator2.C \
      ex-Graeffe.C \
      ex-hilbert1.C  ex-HomomorphismFns1.C \
      ex-ideal1.C  ex-ideal2.C  ex-IdealOfPoints1.C  ex-interrupt1.C  ex-interrupt2.C  ex-interrupt3.C  ex-IntegrationUIBCToSparsePolyRing.C \
      ex-Janet1.C  ex-Janet2.C \
      ex-LogStream1.C \
      ex-LRSDegeneracy.C \
      ex-matrix1.C  ex-matrix2.C  ex-matrix3.C  ex-matrix4.C  ex-module1.C  ex-module2.C \
      ex-MathSat1.C  ex-MathSat2.C \
      ex-MorseGraph.C \
      ex-MVT.C \
      ex-NF.C  ex-Normaliz1.C  ex-Normaliz2.C \
      ex-NumTheory1.C  ex-NumTheory2.C  ex-NumTheory3.C  ex-NumTheory4.C \
      ex-obsolescent.C  ex-OrderingGrading1.C  ex-OrthogPolys1.C \
      ex-PolyInput1.C  ex-PolyInput2.C  ex-PolyIterator1.C  ex-PolyIterator2.C \
      ex-PolyRing1.C  ex-PolyRing2.C  ex-PolyRing3.C  ex-PolyRing4.C \
      ex-Pommaret1.C  ex-Pommaret2.C \
      ex-PPMonoidElem1.C  ex-PPMonoidElem2.C \
      ex-PPMonoidHom1.C \
      ex-PPVector1.C \
      ex-PPWithMask1.C  ex-PPWithMask2.C \
      ex-ProgressReporter1.C \
      ex-QuotientBasis.C \
      ex-RandomSource1.C  ex-RandomSource2.C \
      ex-RandomBool1.C  ex-RandomLong1.C  ex-RandomBigInt1.C \
      ex-ring1.C  ex-ring2.C \
      ex-RingElem1.C  ex-RingElem2.C \
      ex-RingFp1.C  ex-RingFp2.C  ex-RingFq1.C \
      ex-RingHom1.C  ex-RingHom2.C  ex-RingHom3.C  ex-RingHom4.C  ex-RingHom5.C  ex-RingHom6.C \
      ex-RingQQ1.C \
      ex-RingTwinFloat1.C  ex-RingTwinFloat2.C  ex-RingTwinFloat3.C  ex-RingTwinFloat6.C \
      ex-RingWeyl1.C  ex-RingWeyl2.C  ex-RingWeyl3.C  ex-RingWeyl4.C  ex-RingWeyl5.C \
      ex-RingZZ1.C \
      ex-RootBound1.C \
      ex-SmallFp1.C  ex-SmallFp2.C  ex-SmallFp3.C \
      ex-SparsePolyOps1.C  ex-SparsePolyOps2.C  ex-SparsePolyOps3.C  ex-SparsePolyOps-eval.C \
      ex-SumBigRat.C  ex-symbol1.C  ex-symbol2.C \
      ex-ToString1.C  ex-ToString2.C \
      ex-UtilsTemplate1.C \
      ex-VectorOperations1.C  ex-verbose1.C  ex-VerificationLevel1.C


EXECS=$(SRCS:.C=)

# Next 3 lines are useful when developing: triggers recompilation of any *.C file when libcocoa changes
ALL_C=$(wildcard *.C)
ALL_EXECS=$(ALL_C:.C=)

.PHONY: default
default: 
	@echo "*** CoCoALib/examples/Makefile: default target ***"
	@(cd $(COCOA_ROOT); $(MAKE) examples)

$(ALL_EXECS): $(COCOA_LIB)



.PHONY: lib
lib:
	@(cd $(COCOA_ROOT); $(MAKE) lib)

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)

# This target should be made only after the CoCoA library has been compiled;
# normally it would be called by the command "make examples" in the CoCoALib root
# directory.
.PHONY: executables
executables: $(EXECS)


.PHONY: check
check: $(EXECS)
	@echo "Running all examples in CoCoALib/examples/.  Please wait."
	@failures=""; \
	 for prog in $(EXECS); \
	 do \
	   /bin/rm -rf errors; \
	   echo "Running $$prog"; \
	   if [ -f "$$prog.in" ]; \
	   then \
	     ./$$prog < $$prog.in > /dev/null 2> errors; \
	   else \
	     ./$$prog < /dev/null > /dev/null 2> errors; \
	   fi; \
	   if [ $$? -ne 0 -o -s errors ]; then failures="$$failures $$prog"; fi; \
	   /bin/rm -rf errors; \
	 done; \
	 . $(COCOA_ROOT)/configuration/shell-fns.sh; \
	 if [ -n "$$failures" ]; then echobox "These examples failed: $$failures"; exit 1; fi; \
	 echobox "Good news: all examples ran successfully."


.PHONY: valgrind
valgrind: executables
	@./ValgrindExamples.sh $(EXECS)


.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f $(EXECS)  ./*.o  core  a.out  gmon.out
	@/bin/rm -f  ./*~  ./.*~  ./.\#*
	@/bin/rm -rf  ./*.dSYM
	@$(MAKE) -f Makefile-indexhtml  clean-local

.PHONY: veryclean  veryclean-local
veryclean: veryclean-local
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-local: clean-local
	@/bin/rm -f $(ALL_EXECS)
	@$(MAKE) -f Makefile-indexhtml  veryclean


.PHONY: index.html
index.html:
	@$(MAKE) -f Makefile-indexhtml


$(COCOA_ROOT)/configuration/autoconf.mk:
	@echo
	@echo "===================================================================="
	@echo ">>>  ERROR: you must first run configure script in directory $(COCOA_ROOT)  <<<"
	@echo "===================================================================="
	@echo
	@exit 1
