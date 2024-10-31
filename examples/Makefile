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


#############################################################################
# Next few lines are for RCS header/log
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/examples/Makefile,v 1.134 2024/02/08 20:46:30 abbott Exp $
# $Log: Makefile,v $
# Revision 1.134  2024/02/08 20:46:30  abbott
# Summary: Added new class TupleIter (redmine 379)
#
# Revision 1.133  2024/01/30 20:47:05  abbott
# Summary: Better error mesg if autoconf.mk does not exist
#
# Revision 1.132  2024/01/21 19:50:33  abbott
# Summary: Added ex-interrupt3
#
# Revision 1.131  2023/12/09 21:39:11  abbott
# Summary: Added example for SumBigRat
#
# * examples/Makefile (Module):
#
# * examples/ex-SumBigRat.C (Module):
#
# Revision 1.130  2023/12/09 20:29:49  abbott
# Summary: Split "eval" fns out of SparsePolyOps-cyclotomic, into SparsePolyOps-eval; added example & doc
#
# Revision 1.129  2023/01/31 21:31:17  abbott
# Summary: Renamed ex-Cyclotomicity to ex-cyclotomic
#
# Revision 1.128  2022/11/22 21:19:27  abbott
# Summary: Added 3 new examples (LrsDegenerate, IsCyclotomic, Graeffe)
#
# Revision 1.127  2022/03/18 09:21:30  abbott
# Summary: Improved spacing
#
# Revision 1.126  2022/03/17 16:12:00  abbott
# Summary: 2 new examples: ex-00-intro, ex-c++-integers
#
# Revision 1.125  2021/10/07 14:27:01  abbott
# Summary: Removed ex-c++.C (superseded by ex-c++-basic.C)
#
# Revision 1.124  2021/09/16 11:50:33  bigatti
# Summary: renamed vector to vector1
#
# Revision 1.123  2021/09/13 13:56:39  abbott
# Summary: Renamed ex-RingElem1 to ex-RingElem2; renamed new ex-c++-RingElem to ex-RingElem1
#
# Revision 1.122  2021/09/03 12:57:47  abbott
# Summary: Added several new ex-c++-* examples (for Kassel CoCoA school Oct 2021)
#
# Revision 1.121  2020/02/14 12:19:56  abbott
# Summary: Added ex-c++-vector2.C
#
# Revision 1.120  2019/10/08 20:28:10  abbott
# Summary: Minor layout change
#
# Revision 1.119  2019/10/02 17:58:54  abbott
# Summary: Replaced source command by more portable . (dot)
#
# Revision 1.118  2019/03/27 13:45:30  bigatti
# (abbott) renamed ex-GCDFreeBasis1.C --> ex-CoprimeFactorBasis1.C
#
# Revision 1.117  2018/08/02 15:02:31  abbott
# Summary: Fixed small bug
#
# Revision 1.116  2018/08/02 14:56:31  abbott
# Summary: Improved so depends only files called ex-*.C (instead of *.C)
#
# Revision 1.115  2018/08/02 11:19:20  bigatti
# -- split Makefile and Makefile-indexhtml
#
# Revision 1.114  2018/07/31 15:40:34  abbott
# Summary: Moved "target" which says all execs depend on library
#
# Revision 1.113  2018/06/25 12:27:12  abbott
# Summary: Added ex-GCDFreeBasis1
#
# Revision 1.112  2018/04/18 15:38:41  abbott
# Summary: Added some spaces
#
# Revision 1.111  2018/04/18 14:12:57  abbott
# Summary: Added ex-BigRatInterval1.C
#
# Revision 1.110  2018/03/14 15:29:08  abbott
# Summary: Added ex-Verificationlevel1
#
# Revision 1.109  2018/02/12 14:47:18  abbott
# Summary: Changed back so that all executables depend on libcocoa.
#
# Revision 1.108  2017/12/18 10:57:19  abbott
# Summary: Changed spacing
#
# Revision 1.107  2017/12/15 16:36:28  bigatti
# -- added examples and doc
#
# Revision 1.106  2017/12/01 22:18:07  abbott
# Summary: Renamed ex-SturmSeq1.C to ex-SparsePolyOps3.C
#
# Revision 1.105  2017/12/01 21:54:47  abbott
# Summary: Removed CopyInfo.C from ALL_SRCS; corrected target CopyInfo
#
# Revision 1.104  2017/12/01 21:29:40  abbott
# Summary: Added new test
#
# Revision 1.103  2017/11/29 20:34:16  abbott
# Summary: Added SturmSeq and NumRealRoots
#
# Revision 1.102  2017/10/17 14:13:16  abbott
# Summary: Added ex-SparsePolyOps1.C
#
# Revision 1.101  2017/10/16 19:46:35  abbott
# Summary: Added example for new fns ChebyshevPoly, HermtePoly, LaguerrePoly
#
# Revision 1.100  2017/09/14 15:54:54  abbott
# Summary: Added RootBound
#
# Revision 1.99  2017/09/06 11:55:25  abbott
# Summary: Added ex-HomomorphismFns1
#
# Revision 1.98  2017/07/15 15:15:44  abbott
# Summary: Added examples for CpuTimeLimit
#
# Revision 1.97  2017/07/08 19:07:48  abbott
# Summary: updated example for interrupt; added new example too.
#
# Revision 1.96  2017/07/03 12:29:05  abbott
# Summary: Added ex-RingHom6.C
#
# Revision 1.95  2017/04/26 14:50:52  abbott
# Summary: New example: ex-ExternalLibs1
#
# Revision 1.94  2017/03/29 16:33:29  abbott
# Summary: Added ex_LogStream
#
# Revision 1.93  2017/02/22 12:31:03  abbott
# Summary: Added ex-c++-class.C
#
# Revision 1.92  2017/02/15 12:21:22  abbott
# Summary: Added 3 new ex-c++*.C examples
#
# Revision 1.91  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.90  2017/02/10 16:31:25  abbott
# Summary: Added new examples for C++: ex-c++-basic, ex-c++-arith, ex-c++-for-loop
#
# Revision 1.89  2017/02/08 17:00:22  abbott
# Summary: Added new examples ex-ideal1 ex-ideal2 ex-IdealOfPoints
#
# Revision 1.88  2017/01/19 16:53:13  abbott
# Summary: Moved target CopyInfo so that target "default" is first; renamed loop variable from F to SrcFile
#
# Revision 1.87  2016/11/11 13:20:37  abbott
# Summary: Added new example ex-verbose1
#
# Revision 1.86  2016/11/04 21:40:20  abbott
# Summary: Added example for obsolescent fns
#
# Revision 1.85  2016/10/14 09:41:31  abbott
# Summary: Changed target name "VelgrindExamples" into "valgrind"
#
# Revision 1.84  2016/09/21 15:09:22  abbott
# Summary: Added two spaces
#
# Revision 1.83  2016/06/10 13:54:14  abbott
# Summary: Added new simpler rule for compiling CopyInfo
#
# Revision 1.82  2016/05/18 12:18:20  abbott
# Summary: Added new lib target (synonym for library)
#
# Revision 1.81  2016/03/30 09:45:03  abbott
# Summary: Added new ex-RingTwinFloat6 (used to be a test)
#
# Revision 1.80  2016/01/27 14:48:14  abbott
# Summary: Added ex-RingFq1
#
# Revision 1.79  2015/12/08 13:56:08  abbott
# Summary: Updated Mario's code!  Very many changes!
#
# Revision 1.78  2015/11/21 19:50:56  abbott
# Summary: Added ex-FloatApprox1
#
# Revision 1.77  2015/11/04 10:15:07  abbott
# Summary: Added new ex-SmallFp2, ex-SmallFp3; removed duplicate ex-ToString*
#
# Revision 1.76  2015/09/02 11:10:17  bigatti
# -- added ex-GFan1
#
# Revision 1.75  2015/06/25 16:06:29  abbott
# Summary: Added ex-NumTheory4
# Author: JAA
#
# Revision 1.74  2015/05/20 16:53:01  abbott
# Summary: Added ex-interrupt1.C
# Author: JAA
#
# Revision 1.73  2014/09/16 10:41:41  abbott
# Summary: Added new fn eratosthenes (with doc, example, test)
# Author: JAA
#
# Revision 1.72  2014/07/31 15:59:59  abbott
# Summary: Replaced ex-io by ex-VectorOperations1
# Author: JAA
#
# Revision 1.71  2014/07/28 14:43:23  abbott
# Summary: Improved the cleaning targets
# Author: JAA
#
# Revision 1.70  2014/07/17 14:49:51  abbott
# Summary: Added ex-matrix4.C
# Author: JAA
#
# Revision 1.69  2014/07/14 15:32:29  abbott
# Summary: Added new example for UtilsTemplate
# Author: JAA
#
# Revision 1.68  2014/07/11 12:12:50  abbott
# Summary: Updated clean-local target
# Author: JAA
#
# Revision 1.67  2014/07/11 11:07:20  abbott
# Summary: Added better error handling to the generation of index.html
# Author: JAA
#
# Revision 1.66  2014/07/03 06:15:50  bigatti
# -- added ex-PPVector1
#
# Revision 1.65  2014/06/13 12:01:53  abbott
# Summary: Changed make variable EXAMPLES into SRCS; replaced phony "index" by real "index.html"
# Author: JAA
#
# Revision 1.64  2014/06/04 10:12:33  bigatti
# -- added ex-ToString1 ex-ToString2
#
# Revision 1.63  2014/05/20 09:02:48  abbott
# Summary: Renamed ex-decimal1 to ex-ToString1; added ex-ToString2
# Author: JAA
#
# Revision 1.62  2014/04/22 14:09:33  abbott
# Summary: Includes new example for ProgressReporter
# Author: JAA
#
# Revision 1.61  2014/04/17 08:39:23  bigatti
# -- added ex-PolyRing4.C (pseudo-constructors)
#
# Revision 1.60  2014/04/09 14:17:13  bigatti
# -- aggiunto ex-DynamicBitset2
#
# Revision 1.59  2014/03/14 11:01:16  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.58  2014/01/30 15:46:02  bigatti
# -- added ex-PolyInput2
#
# Revision 1.57  2013/06/27 16:46:37  abbott
# Added example for Mario Albert's resolution/morse code.
#
# Revision 1.56  2013/05/31 13:34:28  abbott
# Inserts index-preamble at the start of index.html
#
# Revision 1.55  2013/05/29 14:41:36  abbott
# Added ex-NumTheory2.C (which had been forgotten).
#
# Revision 1.54  2013/05/28 13:29:25  abbott
# Added ex-PolyRing3.C
#
# Revision 1.53  2013/05/28 07:08:54  bigatti
# -- added ex-module2.C
#
# Revision 1.52  2013/05/27 12:55:04  abbott
# Added new example ex-SmallFp1.C
#
# Revision 1.51  2013/03/26 14:54:33  abbott
# Added new example for conversion fns.
#
# Revision 1.50  2012/12/05 11:01:59  abbott
# Renamed the examples for random sequence generators.
#
# Revision 1.49  2012/10/19 16:29:10  bigatti
# -- added ex-Janet 1,2,3 (by Mario Albert)
#
# Revision 1.48  2012/10/05 06:46:21  bigatti
# -- added ex-geobucket1
#
# Revision 1.47  2012/07/04 12:16:22  abbott
# Added ex-matrix3.C
#
# Revision 1.46  2012/05/29 16:34:52  abbott
# Added ex-bool3.C
#
# Revision 1.45  2012/05/11 10:12:42  bigatti
# -- moved ex-C++ as second example (after ex-empty)
#
# Revision 1.44  2012/05/11 10:09:16  bigatti
# -- added ex-c++.C
#
# Revision 1.43  2012/05/10 14:42:17  abbott
# Added new example for anonymous symbols.
#
# Revision 1.42  2012/03/30 10:36:47  abbott
# Renamed ex-BigIntPrime1 to ex-BigInt3.
#
# Revision 1.41  2012/02/10 13:27:01  bigatti
# -- changed ex-RingQ/Z1  --> ex-RingQQ/ZZ1
#
# Revision 1.40  2012/02/03 10:36:26  bigatti
# -- added ex-RandomSource2
#
# Revision 1.39  2011/12/23 15:30:30  bigatti
# -- added ex-matrix2
#
# Revision 1.38  2011/08/26 10:20:11  bigatti
# -- renamed zz->BigInt, QQ->BigRat
#
# Revision 1.37  2011/08/02 13:27:47  abbott
# Added new example about use of decimal conversion fns.
#
# Revision 1.36  2011/07/19 16:21:23  bigatti
# -- added ex-frobby1
#
# Revision 1.35  2011/05/16 13:11:58  abbott
# Target veryclean target now deletes ALL_EXECS (and not just the "official" examples)
#
# Revision 1.34  2011/05/13 17:18:15  abbott
# Fixed clean-local target in Makefile
#
# Revision 1.33  2011/05/03 12:12:04  abbott
# Very minor "cosmetic" changes.
#
# Revision 1.32  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.31  2010/12/17 16:00:49  abbott
# Added new example for RandomSource.
#
# Revision 1.30  2010/10/22 09:13:39  abbott
# Added new GMPAllocator example (and renamed existing one).
# Added a check for output on cerr/clog -- any output means the example "failed".
#
# Revision 1.29  2010/10/08 22:05:36  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.28  2010/10/08 11:00:31  bigatti
# -- added "wildcard" for convenient compilation of non-official
#    examples (warning: working only with gmake?)
# -- added ex-Normaliz2.C
#
# Revision 1.27  2010/10/07 12:57:36  bigatti
# -- change: make --> $(MAKE)
#
# Revision 1.26  2010/07/16 15:47:34  abbott
# Added simple example for PPMonoidHoms.
#
# Revision 1.25  2010/06/29 15:17:25  abbott
# Added new example for RandomZZStream
#
# Revision 1.24  2010/05/28 16:02:02  bigatti
# -- added ex-DynamicBitset1.C
# -- aligned structure with tests/Makefile
# -- fixed sorting
#
# Revision 1.23  2010/03/30 15:23:23  bigatti
# -- removed indexes from default target
#
# Revision 1.22  2010/03/22 10:17:41  abbott
# Added new example for ZZ: actually renamed ex-ZZ1 to ex-ZZ2, & added new ex-ZZ1
#
# Revision 1.21  2010/03/18 17:36:12  bigatti
# -- added ex-PolyRing2.C
# -- fixed minor mis-sorting
#
# Revision 1.20  2010/03/11 16:33:07  bigatti
# -- fixed internal links in index.html
#
# Revision 1.19  2010/03/11 15:36:49  bigatti
# -- sorted examples alphabetically (except ex-empty.C)
#
# Revision 1.18  2010/03/11 15:32:47  bigatti
# -- new: creation of index.html
#
# Revision 1.17  2010/03/03 10:41:50  abbott
# Added example for basic number theory functions.
#
# Revision 1.16  2010/02/16 10:19:29  abbott
# Added new class RandomLongStream with example and test.
#
# Revision 1.15  2010/02/09 10:34:50  abbott
# Added missing semicolon (though it seemed to work OK even without it).
#
# Revision 1.14  2010/02/04 10:41:43  bigatti
# -- changes: default "make" makes indexes
#
# Revision 1.13  2010/02/01 22:37:43  abbott
# Added new examples for DivMask and PPWithMask
#
# Revision 1.12  2009/10/29 19:10:11  abbott
# Now uses echobox if an error occurs when running the examples.
#
# Revision 1.11  2009/07/08 12:26:53  abbott
# Added floor and ceil functions for QQs.
# Added example program for QQs.
# Minor correction to convert.C; minor cleaning to ex-ZZ1.C
#
# Revision 1.10  2009/07/02 16:22:42  abbott
# Better clean target.
# Added ex-factor1 example to build list.
#
# Revision 1.9  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.8  2008/12/12 11:32:02  abbott
# Updated Makefiles to make the new test/example for symbol visible.
#
# Revision 1.7  2008/11/19 09:21:47  bigatti
# -- added ex-PPMonoidElem2
#
# Revision 1.6  2008/07/16 10:03:54  bigatti
# -- added ex-MVT.C
#
# Revision 1.5  2007/06/21 21:29:47  abbott
# Changed name of RingFloat into RingTwinFloat.
#
# Revision 1.4  2007/06/06 15:41:08  abbott
# Fixed typo in name of ex-RandomBitStream1.C
#
# Revision 1.3  2007/06/06 15:16:48  abbott
# Added new RandomBitStream class (now based upon GMP's random generator).
# Consequential changes to Makefiles, etc.  There's even doc and an example!
#
# Revision 1.2  2007/05/14 16:35:36  bigatti
# -- removed examples for Dortmund library
#
# Revision 1.1.1.1  2007/03/09 15:16:11  abbott
# Imported files
#
