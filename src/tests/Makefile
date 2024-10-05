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

#############################################################################
# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/Makefile,v 1.145 2024/08/01 08:46:32 bigatti Exp $
# $Log: Makefile,v $
# Revision 1.145  2024/08/01 08:46:32  bigatti
# Summary: added test-radical1 for radical(I)
#
# Revision 1.144  2024/02/23 11:08:09  abbott
# Summary: Changed actual output filenames to XYZ.cout-found and XYZ.cerr-found (easier to understand)
#
# Revision 1.143  2023/01/31 21:35:50  abbott
# Summary: Renamed test-Cyclotomicity to test-cyclotomic
#
# Revision 1.142  2023/01/18 17:16:04  abbott
# Summary: new tests for LRSDegeneracy
#
# Revision 1.141  2023/01/01 11:45:19  abbott
# Summary: Added new test-combinatorics1
#
# Revision 1.140  2022/11/22 21:23:45  abbott
# Summary: Added new tests for IsCyclotomic, IsLRSDegenerate, IsLRSDegenerateOrder, GraeffeN
#
# Revision 1.139  2021/09/16 12:02:23  abbott
# Summary: Added test for LinKerZZ (redmine 1590)
#
# Revision 1.138  2020/12/04 10:17:44  abbott
# Summary: Improved targets lib, library (now work as expected)
#
# Revision 1.137  2020/10/13 18:37:02  abbott
# Summary: Added test-RingTwinFloat3
#
# Revision 1.136  2020/09/28 11:22:27  abbott
# Summary: Replaced test-bugNNN by test-exbugs
#
# Revision 1.135  2020/03/06 20:35:15  abbott
# Summary: Improved behaviour of building Makefile_dependencies (even when that file does not exist)
#
# Revision 1.134  2020/02/12 15:41:51  abbott
# Summary: Merged test-bug4 and test-bug5 into test-bug2 and test-bug3
#
# Revision 1.133  2020/01/27 20:00:50  abbott
# Summary: Now automatically builds CoCoALib (if nec); minor tidying.
#
# Revision 1.132  2020/01/24 21:41:38  abbott
# Summary: Cleaned up; moved include DEPEND_FILE to bottom (works better there)
#
# Revision 1.131  2020/01/14 21:03:24  abbott
# Summary: Commented out dependency of Makefile_depend on configuration/version
#
# Revision 1.130  2019/12/18 08:16:49  abbott
# Summary: Revised configure; now must explicitly request to build Qt GUI; changed VERSION into COCOALIB_VERSION; conseq changes to Makefiles
#
# Revision 1.129  2019/10/17 20:22:58  abbott
# Summary: Reduced number of test files by putting several smaller tests into the same file (test-NumTheory, test-SparsePolyRing)
#
# Revision 1.128  2019/10/08 20:32:10  abbott
# Summary: Made clean-local remove also gmon.out
#
# Revision 1.127  2019/09/30 17:36:23  abbott
# Summary: Merged test-bug6, test-bug7, test-bug8, test-bug9 into test-bug1
#
# Revision 1.126  2019/03/04 11:42:12  abbott
# Summary: Added new test-NumTheory9
#
# Revision 1.125  2018/06/12 14:11:32  abbott
# Summary: Added new test-bug9 (from redmine 1185)
#
# Revision 1.124  2018/04/20 13:18:28  abbott
# Summary: Added new test-SturmSeq1.C
#
# Revision 1.123  2018/03/20 15:44:38  bigatti
# -- added test-MinPoly2 (with Alice Moallemy)
#
# Revision 1.122  2018/03/02 13:44:55  abbott
# Summary: Added test-NumTheory8
#
# Revision 1.121  2018/02/27 17:33:25  abbott
# Summary: Renamed two test files
#
# Revision 1.120  2017/12/18 13:10:32  bigatti
# -- HomomorphismOps renamed Ops from Fns
#
# Revision 1.119  2017/12/18 10:58:40  abbott
# Summary: Improved header inside Makefile_dependecies
#
# Revision 1.118  2017/12/15 16:33:34  abbott
# Summary: Removed useless ALL_EXECS (and ALL_C) variables
#
# Revision 1.117  2017/12/12 14:24:19  abbott
# Summary: Added test-RootBound2
#
# Revision 1.116  2017/12/01 21:40:58  abbott
# Summary: Added new tests
#
# Revision 1.115  2017/11/29 20:34:17  abbott
# Summary: Added SturmSeq and NumRealRoots
#
# Revision 1.114  2017/11/08 13:43:45  abbott
# Summary: Added new test for IsPrime of a RingElem (in ZZ): used to fail for negative values
#
# Revision 1.113  2017/09/14 15:54:55  abbott
# Summary: Added RootBound
#
# Revision 1.112  2017/09/12 15:06:33  abbott
# Summary: Added new test-MinPoly1.C  (see redmine 1101)
#
# Revision 1.111  2017/08/08 13:56:19  abbott
# Summary: Added new test for HomomorphismFns
#
# Revision 1.110  2017/07/21 13:27:46  abbott
# Summary: Building of Makefile_dependencies should now be "clean" even when interrupted
#
# Revision 1.109  2017/07/19 16:41:53  abbott
# Summary: Added test for RadicalMembership (IsInRadical, MinPowerInIdeal)
#
# Revision 1.108  2017/07/03 19:55:03  abbott
# Summary: Added new test-SparsePolyRing4
#
# Revision 1.107  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.106  2016/11/04 21:41:38  abbott
# Summary: Added new test for calling obsolescent fns
#
# * src/tests/Makefile:
#
# * src/tests/test-obsolescent.C:
#
# * src/tests/test-obsolescent.err: dded
#
# Revision 1.105  2016/11/03 12:26:34  abbott
# Summary: Added test-bug8
#
# Revision 1.104  2016/10/27 14:13:22  abbott
# Summary: Added new test-MatrixSpecial2 (for RandomUnimodularMat)
#
# Revision 1.103  2016/10/25 20:56:05  abbott
# Summary: Added new tests for IsSqFree
#
# Revision 1.102  2016/10/20 18:06:54  abbott
# Summary: Added test-NumTheory6 and test-SparsePolyRing3 (both for radical)
#
# Revision 1.101  2016/10/08 19:48:50  abbott
# Summary: Added new test-RingElem5
#
# Revision 1.100  2016/09/21 15:10:02  abbott
# Summary: Added test-RingAssign1
#
# Revision 1.99  2016/09/19 06:01:24  bigatti
# -- added test-GFan1
#
# Revision 1.98  2016/09/16 13:12:48  abbott
# Summary: Added test-bool3
#
# Revision 1.97  2016/07/21 14:37:13  abbott
# Summary: Added test-BigRat4 for reading BigRats
#
# Revision 1.96  2016/06/29 10:55:03  abbott
# Summary: Added new test test-bug7 (redmine 858)
#
# Revision 1.95  2016/05/18 12:30:29  abbott
# Summary: Added new target "lib": synonym for "library"
#
# Revision 1.94  2016/04/27 16:52:21  abbott
# Summary: new test related to redmine 870
#
# Revision 1.93  2016/03/30 09:48:03  abbott
# Summary: Removed test-RingTwinFloat6 (now it is an example)
#
# Revision 1.92  2016/03/25 20:19:21  abbott
# Summary: Added new tests: test-OrderedDomain1 & test-OrderedDomain2
#
# Revision 1.91  2016/01/27 14:04:22  abbott
# Summary: Added new test-RingFq1
#
# Revision 1.90  2015/12/08 13:56:10  abbott
# Summary: Updated Mario's code!  Very many changes!
#
# Revision 1.89  2015/11/25 13:00:54  abbott
# Summary: New test (see redmine issue 814)
#
# Revision 1.88  2015/11/24 12:51:08  abbott
# Summary: Removed test-F5 (as it was too troublesome)
#
# Revision 1.87  2015/11/04 12:17:57  abbott
# Summary: Added new tests test-SmallFp1 and test-SmallFp2
#
# Revision 1.86  2015/09/03 14:15:44  bigatti
# -- added test-normaliz2 (by Christof Soeger, in Aarhus)
#
# Revision 1.85  2015/06/25 16:10:29  abbott
# Summary: Added new test-NumTheory5
# Author: JAA
#
# Revision 1.84  2014/09/16 10:41:41  abbott
# Summary: Added new fn eratosthenes (with doc, example, test)
# Author: JAA
#
# Revision 1.83  2014/09/05 14:26:04  bigatti
# -- added test-gsl1
#
# Revision 1.82  2014/07/31 16:33:03  bigatti
# -- target ValgrindTests renamed --> valgrind
#
# Revision 1.81  2014/07/28 15:01:21  abbott
# Summary: Improved targets (very)clean-local; new scheme for building dependencies; added test-RingHom4
# Author: JAA
#
# Revision 1.80  2014/06/25 10:38:42  abbott
# Summary: Added new test-RingElem4
# Author: JAA
#
# Revision 1.79  2014/06/17 10:20:36  abbott
# Summary: Depend_file now depends on library.H (instead of version)
# Author: JAA
#
# Revision 1.78  2014/06/14 19:38:21  abbott
# Summary: Added new test-BigRat3 (CmpAbs for BigRat)
# Author: JAA
#
# Revision 1.77  2014/04/17 16:28:54  abbott
# Summary: Corrected name of test-SqFreeFactor1.C
# Author: JAA
#
# Revision 1.76  2014/04/17 06:19:53  bigatti
# -- added test-matrix4 for DeleteRow DeleteCol
#
# Revision 1.75  2014/04/10 15:37:43  abbott
# Summary: Added new test-FloatApprox1
# Author: JAA
#
# Revision 1.74  2014/04/03 15:36:46  abbott
# Summary: Added 2 new tests for utils
# Author: JAA
#
# Revision 1.73  2014/03/27 16:40:51  bigatti
# -- added test-ideal2
#
# Revision 1.72  2014/03/21 17:22:28  abbott
# Summary: New tests for inputting symbols
# Author: JAA
#
# Revision 1.71  2014/03/14 11:04:18  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.70  2014/03/06 16:38:27  abbott
# Summary: Cleaned up a bit  (two spaces between file names in list of test sources)
# Author: JAA
#
# Revision 1.69  2014/03/06 16:00:25  abbott
# Summary: Added new test-IntOperations2 (zero to power zero)
# Author: JAA
#
# Revision 1.68  2013/10/22 14:03:14  abbott
# Added test-SqfreeFactor1
#
# Revision 1.67  2013/06/27 16:57:46  abbott
# Added new tests for Mario Albert's resolution/morse impl.
#
# Revision 1.66  2013/06/17 09:01:32  abbott
# Added new test-RingElem3.
#
# Revision 1.65  2013/05/29 17:12:40  bigatti
# -- added test-QuotientRing1
#
# Revision 1.64  2013/05/28 16:12:46  abbott
# Added test-SparsePolyRing2.C
#
# Revision 1.63  2013/05/28 10:18:00  bigatti
# -- sorted
#
# Revision 1.62  2013/05/27 16:40:57  bigatti
# -- first import of test-FreeModule1, test-FreeModule2
# -- removed test-GradedFreeModule1
#
# Revision 1.61  2013/05/24 12:59:04  abbott
# Added test-FractionField1.
#
# Revision 1.60  2013/04/17 09:09:58  bigatti
# -- added test-JBMill1
#
# Revision 1.59  2013/02/14 17:37:12  bigatti
# -- added test for elimination matrices
#
# Revision 1.58  2013/01/21 15:30:36  abbott
# Renamed files called BuchbergerMoeller* to IdealOfPoints*.
#
# Revision 1.57  2013/01/21 13:40:55  abbott
# Added new test for Buchberger-Moeller & IdealOfPoints.
#
# Revision 1.56  2012/12/05 11:04:47  abbott
# Renamed existing tests for random sequence generators.
#
# Revision 1.55  2012/12/04 09:55:47  abbott
# Added new fns LeastNNegRemainder and SymmRemainder (with doc & new test).
# Some minor corrections to the doc (for operator/ and operator%).
#
# Revision 1.54  2012/10/03 15:21:35  abbott
# Added new test-RingElem2
#
# Revision 1.53  2012/06/29 15:19:57  abbott
# Added new test test-NumTheory3 (after fixing bug #199)
#
# Revision 1.52  2012/06/14 14:44:20  abbott
# Added new test for ring elems.
#
# Revision 1.51  2012/05/11 16:55:37  bigatti
# -- added test-GOperations2.C
#
# Revision 1.50  2012/04/13 16:26:14  abbott
# Added new test-matrix3.C for new fns solve/SolveByGauss.
#
# Revision 1.49  2012/04/02 15:06:10  bigatti
# -- added test-PPMonoidSparse1 (subset of test-PPMonoid1)
#
# Revision 1.48  2012/04/02 13:18:53  bigatti
# -- added test-MatrixForOrdering1
#
# Revision 1.47  2012/02/14 16:01:24  bigatti
# -- added test-PPMonoidHom1.C
#
# Revision 1.46  2011/12/05 16:57:27  bigatti
# -- changed: MachineInteger --> MachineInt (just in comment)
#
# Revision 1.45  2011/11/09 16:38:29  bigatti
# -- renamed MachineInteger --> MachineInt
#
# Revision 1.44  2011/08/27 23:08:32  abbott
# Added new test for "cmp" on MachineInts.
#
# Revision 1.43  2011/08/25 06:24:58  bigatti
# -- renamed test for ZZ and QQ
#
# Revision 1.42  2011/07/20 13:02:31  abbott
# Removed the word "version" from logging message when rebuilding dependencies.
#
# Revision 1.41  2011/07/20 09:55:51  abbott
# Corrected infinite recursion bug (when Makefile_dependencies doesn't exist).
#
# Revision 1.40  2011/07/20 09:04:19  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.39  2011/07/19 16:23:41  bigatti
# -- added test-normaliz1
#
# Revision 1.38  2011/05/09 14:50:06  bigatti
# -- added test-toric1.C
#
# Revision 1.37  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.36  2011/02/10 15:49:53  bigatti
# -- added test-MatrixSpecial1
#
# Revision 1.35  2011/01/14 17:21:43  abbott
# Added new test test-ZZ3
#
# Revision 1.34  2010/12/17 16:04:53  abbott
# Added test for new class RandomSource
#
# Revision 1.33  2010/11/22 17:48:46  bigatti
# -- added test-output
#
# Revision 1.32  2010/10/07 12:58:07  bigatti
# -- change: make --> $(MAKE)
#
# Revision 1.31  2010/06/29 15:18:29  abbott
# Added test for RandomZZStream
#
# Revision 1.30  2010/05/28 15:56:57  bigatti
# -- sorted alphabetically
# -- aligned structure with examples/Makefile
#
# Revision 1.29  2010/05/21 16:18:19  bigatti
# -- added test-DynamicBitset1
#
# Revision 1.28  2010/05/07 14:57:52  abbott
# Two main changes:
#   power(QQ,ZZ) now allows negative exponent
#   renamed QQ::AlreadyNormalized to QQ::AlreadyReduced
#           (and allowed denoms to be negative; the ctor then makes them positive).
#
# Revision 1.27  2010/03/23 15:50:46  bigatti
# -- commented out test-RingTwinFloat3
#
# Revision 1.26  2010/03/23 12:00:49  abbott
# Added test-ZZ2.
#
# Revision 1.25  2010/03/22 11:51:59  abbott
# Added test for ZZ.
#
# Revision 1.24  2010/03/05 21:38:54  abbott
# Added new NumTheory test.
#
# Revision 1.23  2010/02/16 10:19:29  abbott
# Added new class RandomLongStream with example and test.
#
# Revision 1.22  2010/02/03 18:01:22  bigatti
# -- added test-PPVector.C
#
# Revision 1.21  2010/02/02 16:46:28  abbott
# Added new test: test-PPMonoid2
#
# Revision 1.20  2010/02/01 22:45:30  abbott
# Redirected output of RunTests.sh to /dev/null (so
# that possibly confusing error messages do no appear
# on screen).
#
# Revision 1.19  2010/01/20 15:58:04  abbott
# Cleaned the test-RingTwinFloat tests (except number 3).
# Renamed the former test 7 to test 6 (since 7 was more
# thorough & covered what 6 did).
#
# Revision 1.18  2010/01/19 17:40:10  abbott
# Added new test-RingTwinFloat7
#
# Revision 1.17  2009/12/29 22:44:32  abbott
# Removed buggy proxy class ZZ::rtn.
# Consequent changes for function prototypes also in NumTheory.
# Corrected some minor buglets in NumTheory.
#
# Revision 1.16  2009/10/08 13:39:47  abbott
# Renamed "round" into "RoundDiv".
# Added some new versions of "RoundDiv".
# Added a test for "RoundDiv".
#
# Revision 1.15  2009/07/02 16:25:34  abbott
# Changed name of test-convert into test-convert1.
# Improved the clean target.
#
# Revision 1.14  2009/06/04 17:23:05  abbott
# Commented out the two RingFloat tests which currently fail.
# They should eventually be reactivated.
#
# Revision 1.13  2009/05/14 09:47:20  abbott
# Removed two spurious @ symbols.
#
# Revision 1.12  2009/02/11 15:11:36  bigatti
# -- added test-frobby1
#
# Revision 1.11  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.10  2008/12/12 11:32:01  abbott
# Updated Makefiles to make the new test/example for symbol visible.
#
# Revision 1.9  2008/10/09 15:10:04  bigatti
# -- added test-MonomialIdeal1
#
# Revision 1.8  2008/07/08 13:32:28  abbott
# Added missing .C suffix in list of sources.
#
# Revision 1.7  2008/06/05 14:20:37  bigatti
# -- added: make library
#
# Revision 1.6  2007/10/05 15:28:55  bigatti
# -- added abstract class DenseUPolyRing for representing dense
#    univariate polynomials
# -- added concrete class RingDenseUPolyClean, the cleanest
#    implementation
# -- just for testing, still horribly inefficient and incomplete
#
# Revision 1.5  2007/06/21 21:29:47  abbott
# Changed name of RingFloat into RingTwinFloat.
#
# Revision 1.4  2007/05/21 12:50:18  abbott
# Added new test for the convert function family.  Works now thanks to the
# new class MachineInt, and the consequent changes to ZZ.
#
# Revision 1.3  2007/05/14 16:37:37  bigatti
# -- removed tests for Dortmund library
#
# Revision 1.2  2007/03/27 16:57:32  bigatti
# -- added test-TmpIsTree
#
# Revision 1.1.1.1  2007/03/09 15:16:12  abbott
# Imported files
#
# Revision 1.16  2007/03/08 16:53:03  cocoa
# Minor improvement to some Makefiles (which were forgotten at last checkin).
#
# Revision 1.15  2007/03/06 12:34:54  bigatti
# -- added test for GMPAllocator (problems before/after GlobalManager)
#
# Revision 1.14  2007/03/05 21:36:52  cocoa
# Minor improvement to Makefiles etc. to avoid problems with outdated dependency files.
#
# Revision 1.13  2007/01/16 14:00:24  bigatti
# -- mainly a cleanup of TmpFGLM and TmpLESystemSolver (by S.Kaspar)
# -- added tests
#
# Revision 1.12  2007/01/08 16:56:17  cocoa
# -- added F5 code (by A.Arri)
#
# Revision 1.11  2006/12/04 13:59:30  cocoa
# -- added test-GReductor1
#
# Revision 1.10  2006/11/29 14:55:47  cocoa
# -- added ValgrindTests for a "quick" check on memory leaks
#    should return only "definitely lost 0" or "no memory leak possible"
#
# Revision 1.9  2006/11/16 18:16:10  cocoa
# -- added test for HilbertNumeratorMod
#
# Revision 1.8  2006/11/14 17:47:59  cocoa
# -- added first test for RingWeyl
#
# Revision 1.7  2006/11/11 14:51:43  cocoa
# Added simple test for the factorizer.
#
# Revision 1.6  2006/11/10 22:19:52  cocoa
# Revised makefiles to accommodate the factorizer code.
#
# Revision 1.5  2006/10/06 10:31:19  cocoa
# Added new test.
#
# Revision 1.4  2006/08/17 10:01:13  cocoa
# -- removed: test-bug2 (not yet checked in)
#
# Revision 1.3  2006/08/07 21:23:25  cocoa
# Removed almost all publicly visible references to SmallExponent_t;
# changed to long in all PPMonoid functions and SparsePolyRing functions.
# DivMask remains to sorted out.
#
# Revision 1.2  2006/07/19 07:18:05  cocoa
# -- added rudimentary test for GradedFreeModule (LPP, LPos, IsHomogeneous)
#
# Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
# Imported files
#
# Revision 1.7  2006/04/28 17:35:47  cocoa
# -- added test-GOperations1
#
# Revision 1.6  2006/04/27 15:16:07  cocoa
# -- added test-SparsePolyRing1
#
# Revision 1.5  2006/04/05 16:42:08  cocoa
# -- new test for homomorphisms on RingDistrMPolyInlFpPP
#
# Revision 1.4  2006/01/19 18:04:22  cocoa
# -- added test for  NF(f, I)
#
# Revision 1.3  2006/01/17 10:23:08  cocoa
# Updated DivMask; many consequential changes.
# A few other minor fixes.
#
# Revision 1.2  2005/10/18 12:06:36  cocoa
# Cleaned Makefiles, and fixed them so they should work wherever
# CoCoALib is unpacked.
#
# Replaced VERSION cpp macro with COCOA_VERSION.
#
# Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
# Imported files
#
# Revision 1.6  2005/10/14 15:25:07  cocoa
# Major tidying and cleaning to small prime finite fields.
# Several consequential changes.  Improved their documentation.
#
# Added Makefile and script to include/CoCoA/ directory to
# keep library.H up to date.
#
# Revision 1.5  2005/10/11 16:37:30  cocoa
# Added new small prime finite field class (see RingFpDouble).
#
# Cleaned makefiles and configuration script.
#
# Tidied PPMonoid code (to eliminate compiler warnings).
#
# Fixed bug in RingFloat::myIsInteger.
#
# Revision 1.4  2005/10/06 16:36:42  cocoa
# Added the capability find out build information at run-time.
# The Makefiles should be a little tidier too.
#
# Revision 1.3  2005/06/22 14:46:19  cocoa
# Improved configuration and build scripts with
# consequential changes to the Makefiles.
#
# Revision 1.2  2005/05/04 16:37:51  cocoa
# -- added tests: test-PPMonoid1 test-matrix2
#
# Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
# Imported files
#
# Revision 1.4  2005/04/01 16:18:20  cocoa
# Friday check-in.  Fixed a bug in the ctor for GeneralQuotientRingImpl.
#
# Revision 1.3  2005/02/11 16:45:24  cocoa
# Removed the useless and misleading functions myInit and myKill
# from the SmallFp*Impl classes; various consequential changes.
#
# Revision 1.2  2005/02/11 14:15:20  cocoa
# New style ring elements and references to ring elements;
# I hope I have finally got it right!
#
# Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
# Imported files
#
# Revision 1.10  2004/11/19 15:14:09  cocoa
# (a) Added new check to MemPool so that it can signal an
#     error if one tries to use a MemPool after it has been
#     destroyed.
# (b) Improved makefile in TEST/ so that it checks output,
#     and prints useful messages if the test fails.
# (c) Tidied ring.txt a bit (still more to do).
#
# Revision 1.9  2004/11/18 18:33:40  cocoa
# Now every ring know its own "one" element (as well as "zero").
# Several consequential changes.
#
# Revision 1.8  2004/11/12 15:49:29  cocoa
# Tidying prior to 0.90 release.
# (a) documentation improved (or marked as poor)
# (b) sundry minor improvements to the code
#
# Revision 1.7  2004/11/11 11:56:09  cocoa
# (1) Tidied makefiles, and introduced common.mki
# (2) Improved several tests, and cleaned them so that they
#     handle sanely any otherwise unhandled exceptions.
#
# Revision 1.6  2004/06/16 16:13:41  cocoa
# Improved I/O facilities with knock-on changes
#
# Revision 1.5  2004/05/24 15:52:13  cocoa
# Major update:
#   new error mechanism
#   many fixes
#   RingHoms almost work now
#   RingFloat much improved
#
# Revision 1.4  2004/03/20 17:46:10  cocoa
# Check in prior to departure to RWCA
#
# Revision 1.3  2004/02/03 16:30:13  cocoa
# makefile recognises more tests now.
#
# Revision 1.2  2003/10/17 10:53:17  cocoa
# No useful comment.
#
# Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
# Imported files
#
