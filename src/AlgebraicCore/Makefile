# Makefile for CoCoALib/src/AlgebraicCore/ directory.

COCOA_ROOT=../..
include $(COCOA_ROOT)/configuration/autoconf.mk
CWD=src/AlgebraicCore/
DEPEND_FILE=Makefile_dependencies

SUBDIRS=TmpFactorDir TmpHilbertDir

FROBBY_FILES=ExternalLibs-Frobby.C
GFAN_FILES=ExternalLibs-GFan.C
GSL_FILES=ExternalLibs-GSL.C
MATHSAT_FILES=ExternalLibs-MathSAT.C
NORMALIZ_FILES=ExternalLibs-Normaliz.C

SRCS=time.C assert.C utils.C  utils-gmp.C  VectorOps.C  ProgressReporter.C  BuildInfo.C \
     obsolescent.C  exception.C  error.C  MemPool.C  GlobalManager.C \
     SignalWatcher.C  interrupt.C  CpuTimeLimit.C  LogStream.C  verbose.C \
     bool3.C  VerificationLevel.C MachineInt.C  BigInt.C  BigIntOps.C  BigRat.C  BigRatOps.C \
     convert.C  ToString.C  FloatApprox.C \
     NumTheory-misc.C  NumTheory-CRT.C  NumTheory-SimplestRat.C  NumTheory-gcd.C  NumTheory-modular.C \
     NumTheory-ContFrac.C   NumTheory-factor.C  NumTheory-CoprimeFactorBasis.C  NumTheory-prime.C \
     NumTheory-RatReconstruct.C  NumTheory-root.C \
     random.C  combinatorics.C \
     BigRatInterval.C \
     RingElemInput.C  symbol.C  OpenMath.C  OpenMathXML.C \
     SmallFpImpl.C  SmallFpLogImpl.C  SmallFpDoubleImpl.C  DUPFp.C \
     SmallFqVecImpl.C  SmallFqLogImpl.C  SmallFqUtils.C \
     ring.C  ring-AutomaticConversion.C  ideal.C  FieldIdeal.C  RingHom.C  CanonicalHom.C   HomomorphismOps.C \
     RingFp.C RingFpLog.C RingFpDouble.C RingZZ.C RingQQ.C RingTwinFloat.C \
     RingFqLog.C  RingFqVec.C \
     QuotientRing.C  FractionField.C  PolyRing.C  PolyRing-content.C \
     SparsePolyRing.C  SparsePolyIter.C  SparsePolyOps-RingElem.C  SparsePolyOps-RingElem-homog.C \
     SparsePolyOps-EvalInterval.C  SparsePolyOps-RingHom.C \
     SparsePolyOps-ideal.C  SparsePolyOps-ideal-monomial.C  \
     SparsePolyOps-ideal-points.C  SparsePolyOps-ideal-ZeroDim.C \
     SparsePolyOps-ideal-radical.C      SparsePolyOps-ideal-FGLM.C \
     SparsePolyOps-ideal-RadicalMembership.C  SparsePolyOps-implicit.C \
     SparsePolyOps-RealRadical.C  SparsePolyOps-SturmSeq.C  \
     SparsePolyOps-vector.C  SparsePolyOps-involutive.C  \
     SparsePolyOps-IsPrimitivePoly.C  SparsePolyOps-cyclotomic.C  SparsePolyOps-eval.C  SparsePolyOps-Graeffe.C  \
     NumTheory-PrimeModSeq.C  SparsePolyOps-LRSDegeneracy.C  \
     TmpPPVector.C \
     RootBound.C \
     degree.C  PPOrdering.C  PPMonoid.C  PPMonoidHom.C \
     DivMask.C  OrdvArith.C  PPWithMask.C \
     PPMonoidEv.C  PPMonoidEvOv.C  PPMonoidOv.C  PPMonoidSparse.C \
     QBGenerator.C  QuasiPoly.C  PolyFamilies.C \
     matrix.C  MatrixFp.C  MatrixOps.C  MatrixOps-arith.C  MatrixOps-det.C  MatrixOps-rank.C  MatrixOps-flatten.C \
     MatrixOps-adj.C  MatrixOps-norms.C  MatrixOps-rref.C  MatrixOps-LinKerZZ.C  MatrixOps-KroneckerProd.C \
     MatrixOps-IsZeroDet.C  DenseMatrix.C  MatrixView.C \
     LinDepMill.C  MatrixForOrdering.C  MatrixSpecial.C  RandomSparseNonSing01Mat.C \
     SparsePolyOps-resultant.C  SparsePolyOps-MinPoly.C \
     DistrMPolyClean.C  RingDistrMPolyClean.C \
     DistrMPolyInlPP.C  RingDistrMPolyInlPP.C \
     DistrMPolyInlFpPP.C  RingDistrMPolyInlFpPP.C \
     UPoly.C  DenseUPolyRing.C  RingDenseUPolyClean.C  DenseUPolyClean.C \
     GBEnv.C \
     TmpGPoly.C  TmpGPair.C  TmpGReductor.C  TmpGRStats.C  TmpGOperations.C \
     SugarDegree.C  SparsePolyOps-reduce.C \
     ReductionCog.C  geobucket.C \
     module.C  FGModule.C  FreeModule.C  submodule.C  QuotientModule.C \
     ModuleOrdering.C \
     ApproxPts.C  ApproxPts2.C  TmpLDP.C \
     RingExtAlg.C  RingWeyl.C \
     SparsePolyOps-hilbert.C TmpToric.C \
     factor.C  factor-AlgExtn.C  RingElemOps-CoprimeFactorBasis.C \
     TmpF5.C  TmpF5Mat.C  TmpIsTree.C \
     TmpMayerVietorisTree.C  TmpMonomialFunctions.C \
     DynamicBitset.C \
     TmpJBDatastructure.C  TmpJBSets.C \
     TmpMorseElement.C  TmpMorsePaths.C  TmpMorseGraph.C  TmpResolutionMinimization.C \
     TmpJBAlgorithm.C  TmpJBMill.C  TmpPBMill.C \
     TmpPartialMorseBetti.C  TmpMorseBetti.C  TmpMorseResolution.C \
     TmpStabilityAlgorithm.C  TmpUniversalInvolutiveBasisContainer.C \
     BenchmarkToolkit.C \
     ExternalLibs.C \
     $(FROBBY_FILES) $(GFAN_FILES) $(GSL_FILES) \
     $(MATHSAT_FILES) $(NORMALIZ_FILES)

# (see below) OBJS depend on $(COCOA_ROOT)/configuration/autoconf.mk
OBJS=$(SRCS:.C=.o)

# The default target compiles debug_new.C and leak_checker.C as
# well as all the files listed in SRCS.
default:
	@echo "*** $(CWD)Makefile: default target ***"
	@cd $(COCOA_ROOT); $(MAKE) library

.PHONY: all
all: check-source-files $(COCOA_LIB) debug_new.o leak_checker

$(COCOA_LIB): $(OBJS) $(SUBDIRS)
	mkdir -p `dirname $(COCOA_LIB)`
	for subdir in $(SUBDIRS); \
	do \
	  echo "Compiling directory $(CWD)$$subdir/"; \
	  touch "$$subdir"; \
	  ( cd "$$subdir"; $(MAKE) -s ); \
	  if [ $$? -ne 0 ]; then echo "***** Compilation of $(CWD)$$subdir/ FAILED *****"; exit 1; fi; \
	done
	echo "doing ar for AlgebraicCore"
	ar $(ARFLAGS) $(COCOA_LIB) $(OBJS)

$(OBJS):  $(COCOA_ROOT)/configuration/autoconf.mk

.PHONY: check-source-files
check-source-files:
	@./CheckSourceFiles.sh $(SRCS)


#######################################################
# Special compilation targets with extra CPP symbols.
BuildInfo.o: BuildInfo.C $(COCOA_ROOT)/include/CoCoA/BuildInfo.H
	@echo "Compiling BuildInfo.o (with special preprocessor flags)"
	PLATFORM=$$(uname -snr); \
	$(COMPILE) -c  -DCOCOA_VERSION="\"$(COCOALIB_VERSION)\""  -DCOCOA_CXX="\"$(CXX)\""  -DCOCOA_CXXFLAGS="\"$(CXXFLAGS)\""  -DCOCOA_PLATFORM="\"$(PLATFORM)\""  -o BuildInfo.o BuildInfo.C

IdealOfPoints.o: IdealOfPoints.C $(COCOA_ROOT)/include/CoCoA/IdealOfPoints.H
	@echo "Compiling IdealOfPoints --  keep your fingers crossed"
	$(COMPILE) -ITmpFactorDir/ -c IdealOfPoints.C

factor.o: factor.C $(COCOA_ROOT)/include/CoCoA/factor.H
	@echo "Compiling factor --  keep your fingers crossed"
	$(COMPILE) -ITmpFactorDir/ -ITmpFactorDir/multivariate/ -c factor.C

TmpHilbert.o: TmpHilbert.C $(COCOA_ROOT)/include/CoCoA/TmpHilbert.H
	@echo "Compiling TmpHilbert --  keep your fingers crossed"
	$(COMPILE) -ITmpHilbertDir/ -c TmpHilbert.C

debug_new.o: debug_new.C $(COCOA_ROOT)/include/CoCoA/debug_new.H
	@echo "Compiling debug_new.o (not part of CoCoALib)"
	$(COMPILE) -c -o debug_new.o debug_new.C

leak_checker: leak_checker.C
	@echo "Compiling leak_checker (not part of CoCoALib)"
	$(COMPILE) -o leak_checker leak_checker.C
	/bin/rm -rf leak_checker.dSYM


.PHONY: clean clean-local clean-subdirs
clean: clean-local clean-subdirs
	@echo "Cleaned CoCoALib/$(CWD)"

clean-local:
	@/bin/rm -f  $(OBJS)  a.out  core  "$(DEPEND_FILE).old"  "$(DEPEND_FILE).new"  leak_checker  debug_new.o
	@/bin/rm -f  ./*~  ./.*~  ./.\#*
	@/bin/rm -rf  ./*.dSYM

clean-subdirs:
	for subdir in $(SUBDIRS); do ( cd $$subdir; $(MAKE) -s clean ); done

.PHONY: veryclean  veryclean-subdirs
veryclean:  clean-local  veryclean-subdirs
	@/bin/rm -f "$(DEPEND_FILE)"
	@true >| "$(DEPEND_FILE)" # Empties $(DEPEND_FILE)
	@echo "Verycleaned CoCoALib/$(CWD)"

veryclean-subdirs:
	for subdir in $(SUBDIRS); do ( cd $$subdir; $(MAKE) -s veryclean ); done


# This target checks whether we should really build DEPEND_FILE...
# (only if COCOA_HDR is newer than DEPEND_FILE).
.PHONY: dependencies
dependencies:
	@if [ \! -s "$(DEPEND_FILE)" -o \! "$(DEPEND_FILE)" -nt "$(COCOA_HDR)" ] ; \
	then \
	  $(MAKE) create-dependencies; \
	fi
	for subdir in $(SUBDIRS); do ( cd $$subdir; $(MAKE) -s dependencies ); done

# This is the target which really does rebuild the DEPEND_FILE
.PHONY: create-dependencies
create-dependencies:
	@/bin/rm -rf $(DEPEND_FILE).new
	@echo "Rebuilding dependencies in <CoCoALib>/$(CWD)"
	@echo "##################################################################"   >> $(DEPEND_FILE).new
	@echo "#        ---  CoCoALib version $(COCOALIB_VERSION)  ---"              >> $(DEPEND_FILE).new
	@echo "# Dependencies for directory $(CWD)"                                  >> $(DEPEND_FILE).new
	@echo "# DO NOT EDIT: this file is automatically generated; to regenerate"   >> $(DEPEND_FILE).new
	@echo "# run \`make dependencies' in the CoCoALib root directory."           >> $(DEPEND_FILE).new
	@echo "##################################################################"   >> $(DEPEND_FILE).new
	@echo                                                                        >> $(DEPEND_FILE).new
	@for srcfile in $(SRCS); \
	do \
	  objfile=`basename $$srcfile .C`.o; \
	  $(CXX) $(CXXFLAGS) $(INCLUDE_PATHS) -ITmpFactorDir -ITmpFactorDir/linalg -ITmpFactorDir/multivariate -MM "$$srcfile" -MQ "$$objfile" >> $(DEPEND_FILE).new; \
	  echo >> $(DEPEND_FILE).new; \
	done
	@if [ -s "$(DEPEND_FILE)" ] ; \
	then \
	  /bin/mv -f "$(DEPEND_FILE)" "$(DEPEND_FILE).old" ; \
	fi
	@/bin/mv -f $(DEPEND_FILE).new  $(DEPEND_FILE)


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
# $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/Makefile,v 1.191 2024/05/13 13:10:41 bigatti Exp $
# $Log: Makefile,v $
# Revision 1.191  2024/05/13 13:10:41  bigatti
# Summary: added FGLM
#
# Revision 1.190  2024/02/05 11:26:40  bigatti
# Summary: added SparsePolyOps-ideal-radical
#
# Revision 1.189  2024/01/30 20:47:57  abbott
# Summary: Better error mesg if autoconf.mk does not exist
#
# Revision 1.188  2023/12/09 20:29:49  abbott
# Summary: Split "eval" fns out of SparsePolyOps-cyclotomic, into SparsePolyOps-eval; added example & doc
#
# Revision 1.187  2023/11/21 22:06:32  abbott
# Summary: Added CompilationPlatform to BuildInfo (not on redmine)
#
# Revision 1.186  2023/03/15 21:33:55  abbott
# Summary: New file SparsePolyOps-IsPrimitivePoly
#
# Revision 1.185  2023/01/31 21:33:35  abbott
# Summary: SparsePolyOps-cyclotomic now contains everything
#
# Revision 1.184  2023/01/01 11:48:22  abbott
# Summary: INserted new file SparsePolyOps-cyclotomic (revised from Nico's impl)
#
# Revision 1.183  2022/12/05 21:43:42  abbott
# Summary: New file MatrixOps-IsZeroDet.C
#
# Revision 1.182  2022/11/16 19:43:51  abbott
# Summary: Added cyclotmicity, graeffe, PrimeModSeq, LRSDegeracy (from Mexis)
#
# Revision 1.181  2022/02/14 14:47:52  bigatti
# Summary: Removed SparsePolyOps-interreduce, added SparsePolyOps-vector
#
# Revision 1.180  2022/02/05 20:00:00  abbott
# Summary: Renamed RadicalMembership files (redmine 1612)
#
# Revision 1.179  2022/01/20 19:16:41  abbott
# Summary: Added new file NumTheory-root.C
#
# Revision 1.178  2021/10/04 08:53:26  abbott
# Summary: Removed apply.H/.C (redmine 1467, 1598)
#
# Revision 1.177  2021/03/01 19:40:48  abbott
# Summary: Added new file PolyRing-content.C
#
# Revision 1.176  2020/10/27 09:42:18  abbott
# Summary: Changed name of make var INCLUDE to INCLUDE_PATHS; removed ref to obsolete var COCOA_CXX_DEFINES
#
# Revision 1.175  2020/10/14 20:01:54  abbott
# Summary: Renamed SparsePolyOps-interreduce to SparsePolyOps-interreduced
#
# Revision 1.174  2020/10/02 18:40:37  abbott
# Summary: Added new fns FlattenByRows & FlattenByCols
#
# Revision 1.173  2020/09/28 11:12:40  abbott
# Summary: Further splitting of MatrixOps (redmine 1196)
#
# Revision 1.172  2020/09/22 18:19:32  abbott
# Summary: Added new files MatrisOps-det, SparsePolyOps-interreduce
#
# Revision 1.171  2020/06/19 14:59:55  abbott
# Summary: Added new fn AutomaticConversionHom
#
# Revision 1.170  2020/05/26 12:06:18  abbott
# Summary: Renamed TensorMat to KroneckerProd; doc & tests updated
#
# Revision 1.169  2020/03/11 17:00:26  abbott
# Summary: Added new fn HomogCompt; split of fns to do with homog polys into new file.  Cleaned includes.
#
# Revision 1.168  2020/03/06 20:35:15  abbott
# Summary: Improved behaviour of building Makefile_dependencies (even when that file does not exist)
#
# Revision 1.167  2020/02/28 08:56:18  abbott
# Summary: Added LinKerZZ
#
# Revision 1.166  2020/02/12 15:38:54  bigatti
# -- made new file SparsePolyIter.C for BeginIter and EndIter
#
# Revision 1.165  2020/01/27 19:54:37  abbott
# Summary: Flags to ar from ARFLAGS
#
# Revision 1.164  2020/01/26 14:30:03  abbott
# Summary: Split NumTheory into several pieces (redmine 1161)
#
# Revision 1.163  2020/01/14 21:03:24  abbott
# Summary: Commented out dependency of Makefile_depend on configuration/version
#
# Revision 1.162  2020/01/09 18:35:41  abbott
# Summary: Added rref (Row reduced echelon form)
#
# Revision 1.161  2019/12/18 08:16:49  abbott
# Summary: Revised configure; now must explicitly request to build Qt GUI; changed VERSION into COCOALIB_VERSION; conseq changes to Makefiles
#
# Revision 1.160  2019/10/15 12:57:55  bigatti
# -- renamed files for ideals
#
# Revision 1.159  2019/10/08 20:30:05  abbott
# Summary: Renamed subresultant to SparsePolyOps-resultant
#
# Revision 1.158  2019/09/23 08:10:24  abbott
# Summary: Renamed OrthogonalPolys to PolyFamilies
#
# Revision 1.157  2019/03/27 14:21:31  bigatti
# (abbott) renamed GCDFreeBasis --> CoprimeFactorBasis
# (abbott) removed ThreadsafeCounter --> using atomic
#
# Revision 1.156  2019/03/18 11:24:58  abbott
# Summary: Split NumTheory into several smaller files
#
# Revision 1.155  2018/08/06 14:09:11  bigatti
# -- removed SparsePolyOps (merged into SparsePolyOps-RingElem)
#
# Revision 1.154  2018/06/25 12:29:17  abbott
# Summary: Added NumTheory-GCDFreeBasis
#
# Revision 1.153  2018/05/22 14:16:39  abbott
# Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
#
# Revision 1.152  2018/05/18 16:42:42  bigatti
# -- added SparsePolyOps-RingElem.C
#
# Revision 1.151  2018/05/18 12:13:36  bigatti
# -- renamed IntOperations --> BigIntOps
#
# Revision 1.150  2018/05/17 16:04:34  bigatti
# -- renamed MatrixOperations --> MatrixOps
# -- renamed VectorOperations --> VectorOps
#
# Revision 1.149  2018/04/20 13:12:15  abbott
# Summary: Added new fn eval(poly, interval) -- still only a prototype
#
# Revision 1.148  2018/04/19 13:38:40  bigatti
# -- added factor-AlgExtn
#
# Revision 1.147  2018/04/18 14:26:40  abbott
# Summary: Added BigRatInterval
#
# Revision 1.146  2018/04/16 21:47:36  bigatti
# - added SparsePolyOps-RingHom
#
# Revision 1.145  2018/04/09 16:35:43  bigatti
# -- settling new names
#
# Revision 1.144  2018/04/06 15:23:38  bigatti
# -- renamed many files with prefix SparsePolyOps
#
# Revision 1.143  2018/03/13 17:34:07  abbott
# Summary: Added new files VerificationLevel
#
# Revision 1.142  2018/03/12 14:45:11  bigatti
# -- Renamed SparsePoly-ideal/involutive into SparsePolyOps-ideal/involutive
#
# Revision 1.141  2018/02/27 17:28:34  abbott
# Summary: Renamed NumTheory_prime to NumTheory-prime
#
# Revision 1.140  2018/02/27 17:26:23  bigatti
# -- added SparsePoly-involutive.C
#
# Revision 1.139  2018/02/27 17:12:37  bigatti
# -- Renamed SparsePolyRing_ideal.C --> SparsePoly-ideal.C
#
# Revision 1.138  2018/02/27 16:36:49  bigatti
# -- split SparsePolyRing.C part on ideals into SparsePolyRing_ideal.C
#
# Revision 1.137  2018/02/27 10:59:41  abbott
# Summary: Added NumTheory_prime; removed RandomSmallPrime (absorbed into NumTheory_prime)
#
# Revision 1.136  2018/01/31 10:01:29  abbott
# Summary: Added new files RingExtAlg (exterior algebra)
#
# Revision 1.135  2018/01/17 11:19:38  abbott
# Summary: Added LinDepMill
#
# Revision 1.134  2017/12/18 13:12:14  abbott
# Summary: Renamed files
#
# Revision 1.133  2017/12/18 13:09:37  bigatti
# -- HomomorphismOps renamed Ops from Fns
#
# Revision 1.132  2017/12/01 21:34:04  abbott
# Summary: Added RealRadical
#
# Revision 1.131  2017/11/29 20:34:17  abbott
# Summary: Added SturmSeq and NumRealRoots
#
# Revision 1.130  2017/11/15 17:01:27  abbott
# Summary: Added 2 new files
#
# Revision 1.129  2017/10/16 19:53:50  abbott
# Summary: Added new fns ChebyshevPoly, HermitePoly, LaguerrePoly
#
# Revision 1.128  2017/09/14 15:55:28  abbott
# Summary: Added RootBound
#
# Revision 1.127  2017/07/21 13:21:23  abbott
# Summary: Split olf interrupt into two ==> new file SignalWatcher; refactored interrupt and CpuTimeLimit
#
# Revision 1.126  2017/07/19 16:39:02  abbott
# Summary: Added RadicalMembership
#
# Revision 1.125  2017/07/15 15:17:48  abbott
# Summary: Added CpuTimeLimit
#
# Revision 1.124  2017/04/26 14:53:39  abbott
# Summary: Added new file ExternalLibs.C
#
# Revision 1.123  2017/03/13 12:20:55  abbott
# Summary: Removed CPPFLAGS_check.C
#
# Revision 1.122  2017/03/01 17:16:24  abbott
# Summary: Added automatic check for some CPP flags (THREADSAFE_HACK)
#
# Revision 1.121  2017/02/24 08:23:01  bigatti
# -- added MATHSAT_FILES
#
# Revision 1.120  2017/02/20 12:08:53  abbott
# Summary: Added new files MatrixFp.C and subresultant.C
#
# Revision 1.119  2017/02/16 10:12:23  bigatti
# -- tolto MatrixFp.C
#
# Revision 1.118  2017/02/14 17:06:29  abbott
# Summary: Updated clean/veryclean targets in all Makefiles
#
# Revision 1.117  2017/02/01 10:37:05  abbott
# Summary: Added GCDFreeBasis
#
# Revision 1.116  2017/02/01 10:16:14  bigatti
# -- added LogStream.C
#
# Revision 1.115  2016/11/11 13:24:33  abbott
# Summary: Added verbose.C
#
# Revision 1.114  2016/11/07 12:21:44  bigatti
# -- removed MatrixFp (work in progress for John Abbott)
#
# Revision 1.113  2016/11/03 12:29:58  abbott
# Summary: Added file for obsolescent fns; also there is a global flag saying whether to give error if calling one.
#
# Revision 1.112  2016/05/23 12:57:39  bigatti
# -- removed MatrixFp (work in progress)
#
# Revision 1.111  2016/05/03 13:29:35  abbott
# Summary: Major update to Mario's Janet/Pommaret code
#
# Revision 1.110  2016/03/25 20:41:15  abbott
# Summary: Renamed utils_gmp to utils-gmp
#
# Revision 1.109  2016/03/18 12:14:26  bigatti
# -- added MinPoly.C
#
# Revision 1.108  2015/12/18 15:25:28  abbott
# Summary: Added source for non-prime finite fields
#
# Revision 1.107  2015/12/08 13:56:09  abbott
# Summary: Updated Mario's code!  Very many changes!
#
# Revision 1.106  2015/12/02 17:06:40  abbott
# Summary: Corrected include dirs for dependencies
#
# Revision 1.105  2015/09/02 09:29:45  bigatti
# -- added GFAN_FILES
#
# Revision 1.104  2015/06/26 14:58:39  abbott
# Summary: Added new file exception.C
# Author: JAA
#
# Revision 1.103  2015/05/20 13:40:01  abbott
# Summary: Added 2 new files: interrupt & combinatorics
# Author: JAA
#
# Revision 1.102  2014/11/04 18:09:15  abbott
# Summary: Added TmpImplicit
# Author: JAA
#
# Revision 1.101  2014/10/30 13:22:50  abbott
# Summary: moved check-source-files to all target from COCOA_LIB (triggered too many rebuilds)
# Author: JAA
#
# Revision 1.100  2014/07/28 15:43:50  abbott
# Summary: Corrected printed message
# Author: JAA
#
# Revision 1.99  2014/07/28 14:47:49  abbott
# Summary: Improved targets (very)clean-local; new scheme for building dependencies
# Author: JAA
#
# Revision 1.98  2014/07/14 10:03:49  abbott
# Summary: Added new QuasiPoly.C file
# Author: JAA
#
# Revision 1.97  2014/06/25 10:35:08  abbott
# Summary: Added missing CXXFLAGS_DEFINES when building dependency file
# Author: JAA
#
# Revision 1.96  2014/06/17 10:08:41  abbott
# Summary: depend_file now depends on library.H (instead of version)
# Author: JAA
#
# Revision 1.95  2014/06/13 12:05:14  abbott
# Summary: Added new (phony) target check-source-files; added two new source files
# Author: JAA
#
# Revision 1.94  2014/05/15 12:32:38  abbott
# Summary: Removed file io.C (no longer needed)
# Author: JAA
#
# Revision 1.93  2014/04/22 13:26:20  abbott
# Summary: Added new file ProgressReporter
# Author: JAA
#
# Revision 1.92  2014/04/17 13:41:16  abbott
# Summary: renamed MatrixViews to MatrixView
# Author: JAA
#
# Revision 1.91  2014/04/11 15:16:32  abbott
# Summary: Renamed TmpFactor.C to factor.C
# Author: JAA
#
# Revision 1.90  2014/04/11 13:33:28  abbott
# Summary: Replaced decimal.C by ToString.C
# Author: JAA
#
# Revision 1.89  2014/04/10 15:45:20  abbott
# Summary: Added FloatApprox.C to sources; updated Makefile_dependencies (& veryclean-local)
# Author: JAA
#
# Revision 1.88  2014/03/14 11:02:24  abbott
# Summary: clean and veryclean targets now delete .*~ files too
# Author: JAA
#
# Revision 1.87  2014/03/07 14:20:33  abbott
# Summary: Added UPoly.C
# Author: JAA
#
# Revision 1.86  2014/03/05 11:13:12  bigatti
# -- undoing wrong commit
#
# Revision 1.84  2014/01/30 15:51:29  bigatti
# -- added RingElemInput.C
#
# Revision 1.83  2013/06/27 16:52:54  abbott
# Added Mario Albert's new resolution/morse impl (still temporary code).
#
# Revision 1.82  2013/06/05 17:45:57  abbott
# Added TmpLDP.C
#
# Revision 1.81  2013/06/03 09:13:35  bigatti
# renamed ModuleTermOrdering into ModuleOrdering
#
# Revision 1.80  2013/05/27 13:09:05  abbott
# Consequential changes from moving all server-related files to src/server/.
#
# Revision 1.79  2013/04/11 15:08:00  abbott
# Major update to IdealOfPoints: first steps towards revitalizing the old BM code.
#
# Revision 1.78  2013/03/15 11:00:50  abbott
# Added check for exponent overflow when powering a PP.
# Merged PPMonoidEv and PPMonoidEvZZ implementations into a single file.
# Implemented new interface for pseudo-ctors for PPMonoidEv which uses a "flag"
# to say whether exponents are big or not.
#
# Revision 1.77  2013/01/23 15:52:05  bigatti
# added BenchmarkToolkit.[CH]
#
# Revision 1.76  2013/01/21 15:30:36  abbott
# Renamed files called BuchbergerMoeller* to IdealOfPoints*.
#
# Revision 1.75  2012/12/04 20:11:18  abbott
# Replaced the 4 source files for random generators, by the single new source file.
#
# Revision 1.74  2012/11/23 17:27:34  abbott
# Added DUPFp and BuchbergerMoeller to standard sources.
#
# Revision 1.73  2012/10/19 16:24:05  bigatti
# -- added      TmpJBDatastructure.C TmpJBSets.C TmpJBEnv.C   (by Mario Albert)
#
# Revision 1.72  2012/10/02 12:56:55  abbott
# Corrected special rule for compiling BuildInfo.
#
# Revision 1.71  2012/10/01 13:58:41  abbott
# Updated command for generating dependency file.
#
# Revision 1.70  2012/05/29 14:56:15  abbott
# Separated ThreadsafeCounter from symbol; also employed it in ring.C.
#
# Revision 1.69  2012/05/28 12:35:28  bigatti
# -- removed personal *.C file
#
# Revision 1.68  2012/05/28 09:18:21  abbott
# Created IntOperations which gathers together all operations on
# integers (both big and small).  Many consequential changes.
#
# Revision 1.67  2012/04/23 09:47:24  abbott
# Added CXXFLAGS and INCLUDE to the command for generating the dependencies file.
#
# Revision 1.66  2012/03/16 15:40:12  abbott
# Merged contents of NumTheoryQQ (continued fraction functions) into NumTheory.
# Merged the doc too.
#
# Revision 1.65  2012/02/10 17:08:44  abbott
# Changed RingZ.C and RingQ.C into RingZZ.C and RingQQ.C
#
# Revision 1.64  2011/11/09 16:27:34  bigatti
# -- renamed MachineInteger --> MachineInt
#
# Revision 1.63  2011/09/23 15:37:24  abbott
# Moved dependency declaration for $(OBJS)
#
# Revision 1.62  2011/09/23 13:20:35  bigatti
# -- QQ.C renamed into BigRat.C
#
# Revision 1.61  2011/08/14 15:52:17  abbott
# Changed ZZ into BigInt (phase 1: just the library sources).
#
# Revision 1.60  2011/08/12 15:21:26  abbott
# Added BigInt impl (derived from ZZ); not used by anyone yet.
#
# Revision 1.59  2011/07/21 09:31:22  bigatti
# -- updated inclusion of NORMALIZ
#
# Revision 1.58  2011/07/20 13:02:32  abbott
# Removed the word "version" from logging message when rebuilding dependencies.
#
# Revision 1.57  2011/07/20 09:55:51  abbott
# Corrected infinite recursion bug (when Makefile_dependencies doesn't exist).
#
# Revision 1.56  2011/07/20 09:04:20  abbott
# New approach to generating Makefile_dependencies: affects almost all Makefiles.
#
# Revision 1.55  2011/06/27 13:44:21  bigatti
# -- added TmpMonomialIdeal
#
# Revision 1.54  2011/06/23 16:07:13  abbott
# Added incomplete but compiling PPMonoidSparse: first prototype,
# simple rather than fast!
#
# Revision 1.53  2011/05/09 14:27:48  bigatti
# -- added TmpToric.C
#
# Revision 1.52  2011/05/03 09:39:24  abbott
# All Makefiles now offer both "clean" and "veryclean" as targets.
#
# Revision 1.51  2011/03/22 20:17:46  abbott
# Removed SmallPrime.C from SRCS.
#
# Revision 1.50  2011/03/21 17:36:39  bigatti
# -- added ExternalLibs-GSL
#
# Revision 1.49  2011/03/01 15:23:43  abbott
# New version to mark addition of MantissaAndExponent, FloatStr.
#
# Revision 1.48  2011/02/10 15:30:59  bigatti
# -- added MatrixSpecial
#
# Revision 1.47  2010/12/17 16:03:15  abbott
# Added new class RandomSource.
#
# Revision 1.46  2010/11/17 15:52:53  abbott
# Removed GMPAllocator.C from list of source files.
#
# Revision 1.45  2010/11/05 14:44:45  bigatti
# -- renamed TmpFrobby.H into ExternalLibs-Frobby.H
# -- renamed TmpNormaliz.H into ExternalLibs-Normaliz.H
#
# Revision 1.44  2010/10/08 22:05:35  abbott
# Removed pointless CoCoALib prefix in definitions of CWD in Makefiles.
#
# Revision 1.43  2010/10/08 14:11:41  abbott
# Makefile cleaning:
#   (a) proper handling of recursive make,
#   (b) better organized targets (to make embedded shell scripts simpler)
#
# Revision 1.42  2010/10/08 13:38:18  bigatti
# -- added TmpNormalizTypes.C
#
# Revision 1.41  2010/10/07 15:41:22  abbott
# Replaced explicit recursive calls to "make" by calls to "$(MAKE)".
#
# Revision 1.40  2010/10/01 16:49:39  bigatti
# -- added NORMALIZ_FILES
#
# Revision 1.39  2010/07/09 17:04:02  abbott
# Added PPMonoid homs.
#
# Revision 1.38  2010/06/29 15:18:08  abbott
# Added RandomZZStream.
#
# Revision 1.37  2010/03/30 15:33:46  bigatti
# -- added DynamicBitset
#
# Revision 1.36  2010/03/23 15:37:00  bigatti
# -- added GBEnv.C
# -- moved "ar" for fasted compilation on Anna's PPC
#
# Revision 1.35  2010/03/23 14:59:55  bigatti
# -- removed RingTwinFloatTZ.C
#
# Revision 1.34  2010/03/18 16:35:08  abbott
# Added NumTheoryQQ files for number theory operations involving QQs.
#
# Revision 1.33  2010/02/16 10:40:09  abbott
# Removed something that shouldn't have been there (well, not yet).
#
# Revision 1.32  2010/02/16 10:19:29  abbott
# Added new class RandomLongStream with example and test.
#
# Revision 1.31  2010/02/03 16:55:19  bigatti
# -- added TmpMonomialFunctions
#
# Revision 1.30  2010/02/03 12:40:15  bigatti
# -- added RegisterServerOpsUser.C
#
# Revision 1.29  2009/12/03 17:26:34  abbott
# Renamed EFGModule to FGModule.
# Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
# Removed some cruft (old code that was not used by anyone).
#
# Revision 1.28  2009/07/02 16:29:42  abbott
# Added new class QQ to represent rational numbers.
# Consequent change to the Makefile.
#
# Revision 1.27  2009/06/05 12:14:55  abbott
# Major change:
#   created new files NumTheory.H/C  which contain basic number theory operations
#   removed several basic number theory operations from ZZ.H/C
#   removed gcd from MachineInteger.H/C
#   changed names of some basic fns:
#      IsPPrime -> IsProbPrime
#      invmod -> InvMod    (changed signature too)
#      powermod -> PowerMod  (changed signature too)
#   added new fns
#      NextProbPrime & PrevProbPrime
#   consequent changes to other code and tests and examples
#
# Revision 1.26  2009/05/20 14:29:16  abbott
# Moved indexing of CoCoALib into the Makefile of the parent directory (namely src/).
#
# Revision 1.25  2009/05/14 09:47:02  abbott
# Changed "ar" flags to make library construction a bit faster (I hope).
#
# Revision 1.24  2009/02/20 13:27:52  bigatti
# -- renamed sugar.[CH] --> SugarDegree.[CH]
#
# Revision 1.23  2009/02/09 13:59:59  bigatti
# -- new abstract class for computing sugar (tested but not yet used)
#
# Revision 1.22  2009/01/30 17:11:14  bigatti
# -- frobby integration
#
# Revision 1.21  2008/12/16 10:14:13  bigatti
# -- changed makefiles for compatibility with Solaris make (no "-C" option)
#
# Revision 1.20  2008/11/23 18:58:32  abbott
# Major overhaul to preprocessing and SOI/NBM code.
# Split SOI/NBM off into a separate file.
# Preprocessing is now "rational" (but internally guided by C++ doubles).
# SOI/NBM now each have 3 similar interfaces: one purely rational, one for
# input which is represented as doubles, and one which converts the input
# to RingTwinFloat values and produces a result which is over some RingTwinFloat
# (the precision is increased automatically until an answer is obtained).
#
# Revision 1.19  2008/05/30 12:49:43  abbott
# Added an aesthetic space.
#
# Revision 1.18  2008/05/27 16:23:54  bigatti
# -- added TmpPPVector and TmpMayerVietorisTree
#
# Revision 1.17  2008/04/21 11:23:11  abbott
# Separated functions dealing with matrices and PPOrderings into a new file.
# Added matrix norms, and completed adjoint.
#
# Revision 1.16  2008/04/18 15:35:57  abbott
# (long overdue) Major revision to matrices
#
# Revision 1.15  2007/11/20 10:01:26  bigatti
# -- change: TmpF5.C --> TmpF5Mat.C  (by Alberto Arri)
# -- updated and improved test-F5.C
#
# Revision 1.14  2007/10/05 15:28:55  bigatti
# -- added abstract class DenseUPolyRing for representing dense
#    univariate polynomials
# -- added concrete class RingDenseUPolyClean, the cleanest
#    implementation
# -- just for testing, still horribly inefficient and incomplete
#
# Revision 1.13  2007/09/24 14:05:31  abbott
# Added check that subdirectory compilation succeeded; previously ignored subdir errors.
#
# Revision 1.12  2007/06/21 21:29:47  abbott
# Changed name of RingFloat into RingTwinFloat.
#
# Revision 1.11  2007/06/06 15:16:48  abbott
# Added new RandomBitStream class (now based upon GMP's random generator).
# Consequential changes to Makefiles, etc.  There's even doc and an example!
#
# Revision 1.10  2007/05/21 12:57:28  abbott
# New class for passing machine integers as args; includes some simple
# operations on machine integers (cmp, gcd, IsNegative,...).  Operations
# between ZZ and machine integers modified to use the new class.  Inexact
# integer division (of a ZZ) by a negative value now triggers an error;
# new error for reporting inexact integer division by a negative value.
#
# Revision 1.9  2007/05/14 16:36:34  bigatti
# -- removed .C files for Dortmund library (Kaspar)
#
# Revision 1.8  2007/05/03 10:37:13  abbott
# Updated Makefile.
#
# Revision 1.7  2007/04/27 14:58:05  bigatti
# -- content of CoCoAServer.C split into dedicated files
# -- new registration mechanism
#
# Revision 1.6  2007/03/23 18:38:43  abbott
# Separated the "convert" functions (and CheckedCast) into their own files.
# Many consequential changes.  Also corrected conversion to doubles.
#
# Revision 1.5  2007/03/09 20:58:36  abbott
# Final retouches: modified ar flags (again), more important change to main Makefile.
#
# Revision 1.4  2007/03/09 19:09:19  abbott
# Changed flags to ar -- still have trouble with numerical stuff vanishing
# from libcocoa.a when something in AlgebraicCore has to be recompiled.
#
# Revision 1.3  2007/03/09 18:56:56  bigatti
# -- added Tmp prefix to Groebner related files
#
# Revision 1.2  2007/03/09 18:21:56  abbott
# Changed ar flags: use r instead of q (which causes trouble if you partially
# recompile).
#
# Revision 1.1.1.1  2007/03/09 15:16:11  abbott
# Imported files
#
#
