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
     module.C  FinGenModule.C  FreeModule.C  submodule.C  QuotientModule.C \
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

# Get current Git revision hash (if possible)
COCOALIB_GIT_REV=$(shell git rev-parse --short HEAD 2> /dev/null)
ifneq ($(COCOALIB_GIT_REV),)
  # Check if the working directory is clean
  GIT_STATUS=$(shell git status --porcelain 2> /dev/null)
  ifneq ($(GIT_STATUS),)
    COCOALIB_GIT_REV+=[modified]
  endif
else
  COCOALIB_GIT_REV=unknown
endif

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
# Saves the current Git revision for BuildInfo.C
.PHONY: update-git-rev
update-git-rev: ;

.gitrev: update-git-rev
	@current_rev="$$(cat $@ 2> /dev/null)"; \
    if [ "$$current_rev" != "$(COCOALIB_GIT_REV)" ]; then \
	  echo "Git revision changed --  updating .gitrev"; \
      echo "$(COCOALIB_GIT_REV)" > $@; \
    fi


#######################################################
# Special compilation targets with extra CPP symbols.
BuildInfo.o: BuildInfo.C $(COCOA_ROOT)/include/CoCoA/BuildInfo.H .gitrev
	@echo "Compiling BuildInfo.o (with special preprocessor flags)"
	PLATFORM=$$(uname -snr); \
	$(COMPILE) -c  -DCOCOA_VERSION="\"$(COCOALIB_VERSION)\""  -DCOCOA_GIT_REVISION="\"$(COCOALIB_GIT_REV)\""  -DCOCOA_CXX="\"$(CXX)\""  -DCOCOA_CXXFLAGS="\"$(CXXFLAGS)\""  -DCOCOA_PLATFORM="\"$(PLATFORM)\""  -o BuildInfo.o BuildInfo.C

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
	@/bin/rm -f  $(OBJS)  a.out  core  "$(DEPEND_FILE).old"  "$(DEPEND_FILE).new"  leak_checker  debug_new.o .gitrev
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
