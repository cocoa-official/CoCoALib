# Makefile for CoCoALib/src/tests/ directory

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/tests/

DEPEND_FILE=Makefile_dependencies

TESTS=test-empty.C \
      test-BigInt1.C  test-BigInt2.C  test-BigInt3.C \
      test-BigRat1.C  test-BigRat2.C  test-BigRat3.C  test-BigRat4.C \
      test-bool3.C \
      test-exbugs.C  \
      test-combinatorics1.C  test-convert1.C \
      test-cyclotomic.C \
      test-DenseUPolyRing1.C \
      test-DynamicBitset1.C \
      test-factor1.C \
      test-FloatApprox1.C \
      test-FractionField1.C \
      test-FreeModule1.C  test-FreeModule2.C \
      test-frobby1.C \
      test-GMPAllocator.C \
      test-GOperations1.C  test-GOperations2.C  test-GReductor1.C \
      test-Graeffe.C \
      test-gsl1.C \
      test-hilbert1.C  test-HomomorphismOps1.C \
      test-ideal1.C  test-ideal2.C  test-IdealOfPoints1.C \
      test-IntOperations1.C  test-IntOperations2.C  test-IntOperations3.C  test-IsInteger1.C \
      test-LRSDegeneracy1.C  test-LRSDegeneracy2.C \
      test-MachineInt1.C  test-MachineInt2.C \
      test-matrix1.C  test-matrix2.C  test-matrix3.C test-matrix4.C \
      test-MatrixForOrdering1.C  test-MatrixForOrdering2.C \
      test-MatrixSpecial1.C  test-MatrixSpecial2.C  test-MatrixOps1.C \
      test-MinPoly1.C test-MinPoly2.C \
      test-module1.C \
      test-MonomialIdeal1.C \
      test-GFan1.C \
      test-normaliz1.C test-normaliz2.C \
      test-NumTheory1.C  test-NumTheory2.C  test-NumTheory3.C  test-NumTheory4.C  test-NumTheory5.C \
      test-obsolescent.C  test-OpenMath1.C \
      test-OrderedDomain1.C  test-OrderedDomain2.C \
      test-output.C \
      test-PPMonoid1.C  test-PPMonoid2.C  test-PPMonoidSparse1.C \
      test-PPMonoidHom1.C \
      test-PPVector.C \
      test-QuotientRing1.C  test-QBGenerator1.C \
      test-radical1.C  test-RadicalMembership1.C \
      test-RandomSource1.C  test-RandomSource2.C  test-RandomSource3.C \
      test-RandomLong1.C test-RandomBigInt1.C \
      test-RingAssign1.C \
      test-RingElem1.C  test-RingElem2.C  test-RingElem3.C  test-RingElem4.C  test-RingElem5.C  test-RingElem6.C \
      test-RingFp1.C \
      test-RingFq1.C \
      test-RingHom1.C  test-RingHom2.C  test-RingHom3.C  test-RingHom4.C \
      test-RingQQ1.C \
      test-RingTwinFloat1.C  test-RingTwinFloat2.C   test-RingTwinFloat3.C \
      test-RingTwinFloat4.C  test-RingTwinFloat5.C \
      test-RingWeyl1.C \
      test-RingZZ1.C  test-RingZZ2.C  test-RingZZmodN1.C \
      test-RootBound1.C  test-RootBound2.C \
      test-SmallFp1.C  test-SmallFp2.C \
      test-SparsePolyRing1.C  test-SparsePolyRing2.C  \
      test-SqFreeFactor1.C \
      test-SturmSeq1.C \
      test-symbol1.C  test-symbol2.C \
      test-toric1.C \
      test-TmpIsTree.C  test-TmpMorseElement.C  test-TmpMorsePath.C  test-TmpMorseGraph.C  test-TmpResolutionMinimization.C  \
      test-utils1.C  test-utils2.C
#     test-RingTwinFloat3.C test-OpenMath2.C   -- skip these tests for the moment


EXECS=$(TESTS:.C=)


default:
	@echo "*** src/tests/Makefile: default target ***"
	@(cd $(COCOA_ROOT); $(MAKE) check)

.PHONY: lib
lib: library

.PHONY: library
library:
	@(cd $(COCOA_ROOT); $(MAKE) library)


# This target should be made only after the CoCoA library has been compiled;
# normally it would be called by the command "make check" in the CoCoALib root
# directory.
.PHONY: executables
executables: $(EXECS)

$(EXECS): $(COCOA_LIB)


# This target should be made only after the CoCoA library has been compiled;
# normally it would be called by the command "make RunTests" in the CoCoALib root
# directory.  Most of the work is done by the RunTests.sh script.
.PHONY: check
check: $(EXECS)
	@./RunTests.sh $(EXECS) 2> /dev/null  # discard err mesgs from the shell script


.PHONY: valgrind
valgrind: executables
	@./ValgrindTests.sh $(EXECS)


.PHONY: clean clean-local
clean: clean-local
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f  ./*.o  gmon.out  core  a.out  ./*~  ./.*~  ./.\#*
	@/bin/rm -f $(EXECS)  "$(DEPEND_FILE).new"  "$(DEPEND_FILE).old"  ./*.cout-found  ./*.cerr-found
	@/bin/rm -rf  ./*.dSYM

.PHONY: veryclean
veryclean: clean-local
	@true >| "$(DEPEND_FILE)" # Empties $(DEPEND_FILE)
	@echo "Verycleaned CoCoALib/$(CWD)"


# Morally the DEPEND_FILE depends on $(TESTS), but I do not tell make
# this as it tends to rebuild it too often.

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
	@echo "##################################################################"    >> $(DEPEND_FILE).new
	@echo "#        ---  CoCoALib version $(COCOALIB_VERSION)  ---"               >> $(DEPEND_FILE).new
	@echo "# Dependencies for directory CoCoALib/src/tests/"                      >> $(DEPEND_FILE).new
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"    >> $(DEPEND_FILE).new
	@echo "# run \`make dependencies' in the CoCoALib root directory."            >> $(DEPEND_FILE).new
	@echo "##################################################################"    >> $(DEPEND_FILE).new
	@echo                                                                         >> $(DEPEND_FILE).new
	@echo "#----------------------------------------------------------"           >> $(DEPEND_FILE).new
	@echo "# all executables also depend on libcocoa.a -- see Makefile"           >> $(DEPEND_FILE).new
	@echo "#----------------------------------------------------------"           >> $(DEPEND_FILE).new
	@echo                                                                         >> $(DEPEND_FILE).new
	@for srcfile in $(TESTS); \
	do \
	  execfile=`basename $$srcfile .C`; \
	  echo "$$execfile: $$srcfile" >> $(DEPEND_FILE).new; \
	  echo >> $(DEPEND_FILE).new; \
	done
	@if [ -s "$(DEPEND_FILE)" ] ; \
	then \
	  /bin/mv -f "$(DEPEND_FILE)" "$(DEPEND_FILE).old" ; \
	fi
	@/bin/mv $(DEPEND_FILE).new  $(DEPEND_FILE)


# Tests using external libs depend on those libs
test-frobby1: $(FROBBY_LIBS)
test-normaliz1: $(NORMALIZ_LIBS)

$(DEPEND_FILE):
	@touch "$(DEPEND_FILE)"

include $(DEPEND_FILE)
