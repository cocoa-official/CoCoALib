#ifndef CoCoA_GReductor_H
#define CoCoA_GReductor_H

//   Copyright (c)  2005-2017  John Abbott, Anna M. Bigatti
//   Author2:  2005-2010  Massimo Caboara, 2010-2017 Anna M. Bigatti 

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


//#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/GBEnv.H"
#include "CoCoA/TmpGPoly.H"
#include "CoCoA/TmpGPair.H"
#include "CoCoA/TmpGRStats.H"
#include "CoCoA/utils.H"
//#include "CoCoA/CoCoA4io.H"


namespace CoCoA
{

  //  void ReadInt(std::istream&,int&,SkipTagType);

  class FreeModule;     // forward declaration -- defined in FreeModule.H
  class SparsePolyRing; // forward declaration -- defined in SparsePolyRing.H

class GBCriteria
{
public:
  enum CoprimeFlag   { UseCoprime, DontUseCoprime };
  enum GMFlag        { UseGM, DontUseGM };
  enum BackFlag      { UseBack, DontUseBack };
  enum DivFlag       { UseDiv, DontUseDiv }; ///< remove poly if its LPP is divisible by LPP of new poly; true except for RingWeyl
  enum AllSetMarker  { AllSet };

public:
  GBCriteria(AllSetMarker) {myBack=true; myCoprime=true; myDiv=true; myGM=true;}
  GBCriteria(CoprimeFlag c, GMFlag gm, BackFlag b, DivFlag d)
  {
    myCoprime = (c == UseCoprime);
    myBack = (b == UseBack);
    myDiv = (d == UseDiv);
    myGM = (gm == UseGM);
  }
public:
  bool myBack;
  bool myCoprime;
  bool myDiv;
  bool myGM;
};

  class GReductor
  {
  public:

    enum UseDynamicAlgFlag    { UseDynamicAlg, DontUseDynamicAlg };
    enum BuchbergerOpTypeFlag { HomogeneousAlg, SaturatingAlg, AffineAlg };


    GReductor(const GRingInfo&,
              const PolyList&,
//              const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit(),
              const BuchbergerOpTypeFlag theBuchbergerOpType=HomogeneousAlg,
              const Reductors::UseBorelFlag=Reductors::DontUseBorel,
              const UseDynamicAlgFlag IsDynamic=DontUseDynamicAlg,
              const GBCriteria criteria=GBCriteria(GBCriteria::AllSet));

    GReductor(const GRingInfo&,
              const GPolyList&,
//              const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit(),
              const BuchbergerOpTypeFlag theBuchbergerOpType=HomogeneousAlg,
              const Reductors::UseBorelFlag=Reductors::DontUseBorel,
              const UseDynamicAlgFlag IsDynamic=DontUseDynamicAlg,
              const GBCriteria criteria=GBCriteria(GBCriteria::AllSet));
    GReductor(const GReductor&);// copy ctor not working, to be fixed
    ~GReductor(){};
    void myCtorAux(const BuchbergerOpTypeFlag theBuchbergerOpType,
                   const UseDynamicAlgFlag IsDynamic);
    GReductor(const GRingInfo&,
              GPolyList&,
              const ClearMarker,
//              const CpuTimeLimit& CheckForTimeout=NoCpuTimeLimit(),
              const BuchbergerOpTypeFlag theBuchbergerOpType=HomogeneousAlg,
              const Reductors::UseBorelFlag=Reductors::DontUseBorel,
              const UseDynamicAlgFlag IsDynamic=DontUseDynamicAlg,
              const GBCriteria criteria=GBCriteria(GBCriteria::AllSet));

    bool IsDynamic() const {return IsDynamicAlgorithm;}
    BuchbergerOpTypeFlag myGetBuchbergerOpType() const {return myBuchbergerOpType;}
    bool myPreparationDone() const {return myPrepared;}
    bool WrongLPPFound() const {return myWrongLPPFoundValue;}
    const Stats& myStats() const {return myStat;}

    long myReductorsLen() const {return myTrueReductors.myLen();}
    long myGBasisLen() const {return len(myGB);}
    long myPairsLen() const {return len(myPairs);}

    void mySetBuchbergerOpType(const BuchbergerOpTypeFlag theBuchbergerOpType) 
        {myBuchbergerOpType=theBuchbergerOpType;}

    void myCopyGBasis(PolyList&);// GB output
    void myCopyMinGens(PolyList&);// GB output
    void myCopyGBasis(GPolyList&);// GB output in GPoly form
    void myGBasisClear(GPolyList&); ///< GB output in GPoly form
    void myCopyGBasis(VectorList& outGBasis);
    void myCopyMinGens(VectorList& outMinGens);
    void myDoAFFGBasis(); // I am working on sugar selection strategies
    void myDoGBasis(); ///< flag may be set for using Borel reductors
    void myDoGBasisTEST(); ///< flag may be set for using Borel reductors
    RingElem myDoGBasisElimFirst(ConstRefPPMonoidElem ElimIndsProd);
    void myReduceCurrentSPoly();
    void myPrepareGBasis(); ///< flag may be set for using Borel reductors
    void myFinalizeGBasis();//Last operations (stats at the moment) immediately before shutting down a GB computation
    void myDoGBasis(const int ReductionNumber); // Performs ReductionNumber reductions, -1 means unlimited
    void myReduceUntilNonZeroRedSPoly(); // Reduces until an SPoly does not reduces to zero. Stop immediately when that happens (no updates)
     // DYNAMIC ALG: reduces until a non zero SPoly has a LPP with is different from the best possibile LPP w.r.t. HPoly.
     // In that case, settles the WrongLPP filed in the reductor.
    void myReduceUntilWrongLPPFound(RefPPMonoidElem,
                                    std::vector<RingElem>& );
//A24    GPoly GetSPoly()const{return mySPoly;};// Reads the SPoly
//A24    void GetCandidateGBasis(PolyList&)const;// the polys computed up to now and the non processed generators
    long myAge() const {return myAgeValue;}
    long myTruncDeg() const {return myTruncDegValue;}
//A24    long GetNReductions() const {return myStat.myNReductions;};
    void mySetTruncDeg(long TruncDeg) { myTruncDegValue = TruncDeg; }
//////    void mySetSPoly(GPoly& p){mySPoly=p;};// Sets the SPoly  not used
    void myDoGBasisRealSolve();
    void myDoGBasisSelfSatCore(); ///< dehomog algorithm
    ////    void _myDoSATGBasis(); ///< dehomog algorithm
    ////    void myDoSATMixGBasis(); ///< dehomog mix algorithm
    void myStampaPPGB(std::ostream&)const; ///< print ?
    void myStampaGB(std::ostream&)const; ///< print ?
    void myStampaPairs(std::ostream&)const; ///< print ?
    void myStampaReductors(std::ostream&)const; ///< print ?
    void myStampaStats(std::ostream& out) const{myStat.myStampa(out);}; ///< print ?
    friend std::ostream& operator<<(std::ostream& out, const GReductor& GR);
    const SparsePolyRing& myPRing()const{return myGRingInfoValue.myNewSPR();};
    const GRingInfo& myGRingInfo()const{return myGRingInfoValue;};
//A24    void Rebuild(const PolyList&);///< rebuild the GReductor initliazinig it with the PL
/// not used?    void myUpdateBasisOnly();///<Updates the Basis only.
    void myCreatePairs();// this should be const, but this requires
                                   //the GPair ctor to be const, and I have to think about that
    void myCreateInputPolyPairs(GPolyList&);
///    void myDoGBasisByBatch();//Tmp, DYNAMIC
///    void myPrepareGBasisPairsExcluded();//Tmp, DYNAMIC
///    void myBuildNewPairsAll(GPairList&);//Tmp, DYNAMIC

  private:
    long myGMInsert(GPairList&,GPair);
    void myBuildNewPairs(GPairList&);
    void myUpdateBasisAndPairs();
    void myUpdateBasisAndPairs(const GPoly& inPoly);
    void myUpdateIncomingWDeg();
    void myApplyBCriterion();

// These two should go in GPoly (interface) with the
// real body in PolyRing
//void smart_dehomog(GPoly&,long);
//void smart_dehomog_DRL(GPoly&,long);

  private:  // member fields
    const GRingInfo myGRingInfoValue;
    GPairList myPairs;// here the polys are ptrs to myPolys
    Reductors myTrueReductors;// the true reductors.
    GPolyPtrList myGB;// the candidate Gbasis - NB polys are ptrs to myPolys
    GPolyList myPolys;// the REAL Polys, the other are ptrs
    GPoly mySPoly;
    degree myOldDeg;// used for flow control and stats. The degree it refers to is the degree of the pair
    degree myIncomingWDeg; // used for flow control and stats
    long myTruncDegValue;
    Stats myStat;// the statistics
    bool myPrepared; // Default false. True after a myPrepareGBasis has been performed
    long myAgeValue;
    bool IsDynamicAlgorithm;
    bool myWrongLPPFoundValue; // DYNAMIC ALGORITHM
    GBCriteria myCriteria;
    BuchbergerOpTypeFlag myBuchbergerOpType; // Type of operation performed on the reductor
    static int ourDefaultStatLevel; ///< default verbosity level for statistics
//    CpuTimeLimit myCheckForTimeout;


  };// End class GReductor

   void monic(PolyList&);

   const GRingInfo& GetGRI(const GPolyList& theGPL);
   FreeModule owner(const VectorList& theVL);

    enum ModOrdTypeForcing {NoForcing,PosWDegTO,WDegTOPos,WDegPosTO};
    SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM);
    SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM,ModOrdTypeForcing MOType);
    SparsePolyRing MakeNewPRingForSimpleEmbedding(const SparsePolyRing& theOldP);
    SparsePolyRing MakeNewPRingForSimpleEmbedding(const SparsePolyRing& theOldP,ModOrdTypeForcing MOType);
    FreeModule MakeNewFreeModuleForSyz(const GPolyList& GPL);
    FreeModule MakeNewFreeModuleForSyz(const VectorList& VL);
    FreeModule MakeNewFreeModuleForSyz(const PolyList& PL);
    SparsePolyRing MakeElimRingFromOld(const SparsePolyRing& theOldP,
                                       const std::vector<long>& theElimVars,
				       const bool IsHomog);
    SparsePolyRing MakeNewPRingFromModulePosFirst(const FreeModule& FM,
					          bool HomogInput);

    SparsePolyRing MakeNewPRingForSimpleEmbeddingPosFirst(const SparsePolyRing& OldP,
                                                          bool HomogInput);
    // Embed p in the CompIndex component, grading given by the embedding
    GPoly EmbedPoly(ConstRefRingElem p,
                    const GRingInfo& theGRI,
                    const long CompIndex);
    // Embed p in the CompIndex component, grading given by the_d
    GPoly EmbedPoly(ConstRefRingElem p,
                    const GRingInfo& theGRI,
                    const degree& the_d,
                    const long CompIndex);

   // Embed v
    GPoly EmbedVector(const ModuleElem& v,
    		      const GRingInfo& theGRI);

   // Embed v starting from the StartingFromCompIndex component
    GPoly EmbedVector(const ModuleElem& v,
                      const GRingInfo& theGRI,
                      const long StartingFromCompIndex);

   // Embed theVL
    GPolyList EmbedVectorList(const VectorList& theVL,
                              const GRingInfo& theGRI);
    // Embed theVL starting from the StartingFromCompIndex component
    GPolyList EmbedVectorList(const VectorList& theVL,
                              const GRingInfo& theGRI,
                              const long StartingFromCompIndex);

    // Just transform the PolyList in a GPolyList
    GPolyList EmbedPolyList(PolyList& thePL,
                            const GRingInfo& theGRI);


    GPolyList EmbedPolyList(const PolyList& thePL,
                            const GRingInfo& theGRI,
                            const degree& the_d,
                            const long CompIndex);

    // Embed with the shifts Shifts
    GPolyList EmbedPolyList(const PolyList& InputPolyList,
                            const GRingInfo& theGRI,
                            const long CompIndex);
			
    // The embedding used in syzygy computations
    GPolyList SyzEmbedVectorList(const VectorList& InputVectorList,
                                 const GRingInfo& theGRI);
  // The embedding used in syzygy computations
    GPolyList SyzEmbedPolyList(const PolyList& InputPolyList,
                                 const GRingInfo& theGRI);
				
    // The embedding used in intersection computations
    GPolyList IntEmbedVectorLists(const VectorList& VL1,
                                  const VectorList& VL2,
                                  const GRingInfo& theGRI);
				
    // The embedding used in intersection computations
    GPolyList IntEmbedPolyLists(const PolyList& PL1,
                                const PolyList& PL2,
                                const GRingInfo& theGRI);
 			
  // The special poly embedding used in colon computations
    GPolyList ColonEmbedVectorLists(const VectorList& VL1,
                                    const VectorList& VL2,
                                    const GRingInfo& theGRI);
				
    GPolyList ColonEmbedPolyLists(const PolyList& PL1,
                                  const PolyList& PL2,
                                  const GRingInfo& theGRI);

    void SyzEmbedGPolyList(GPolyList& theGPL);

    void IntEmbedGPolyList(GPolyList& theGPL1, GPolyList& theGPL2);

    void ColonEmbedGPolyList(GPolyList& theGPL, GPoly& the_gp);

  ModuleElem DeEmbedPoly(ConstRefRingElem p,
                         const GRingInfo& theGRI,
                         const long ComponentsLimit); // the component in p that goes to the 0 component of the output vector v. Lesser components of p go to higher component of v
		
    ModuleElem DeEmbedPoly(ConstRefRingElem p,
                           const GRingInfo& theGRI);

    VectorList DeEmbedPolyList(const PolyList& PL,
                               const GRingInfo& theGRI);

    VectorList DeEmbedPolyList(const PolyList& PL,
                               const GRingInfo& theGRI,
                               const long ComponentsLimit); // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding



    void DeEmbedPoly(ModuleElem& theOutputP,
                     const GPoly& the_gp,
                     const long ComponentsLimit); // the component in p that goes to the 0 component of the output vector v. Lesser components of p go to higher component of v
		
    void DeEmbedPoly(ModuleElem& theOutputP,
                     GPoly& the_gp);

    void DeEmbedPolyList(VectorList& theOutputVL,
                         GPolyList& theGPL);

    void DeEmbedPolyList(VectorList& theOutputVL,
                         const GPolyList& theGPL,
                         const long ComponentsLimit); // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding

    void DeEmbedPoly(RingElem& theOutputP,
                      GPoly& the_gp); // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding

    void DeEmbedPolyList(PolyList& theOutputPL,
                         GPolyList& theGPL,
                         const long ComponentsLimit); // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding

    // input are embedded polys, output true (OldP) Polys
    // this is done directly and not passing through a VectorList to avoid copying
    RingElem DeEmbedPolyToP(ConstRefRingElem the_p,
                            const GRingInfo& theGRI);

    PolyList DeEmbedPolyListToPL(const PolyList& PL,
                                 const GRingInfo& theGRI,
                                 const long ComponentsLimit); // Poly whose LPP has last var degree bigger than this number disappear on DeEmbedding


    void PolyList2GPolyList(const PolyList&,GPolyList&,const GRingInfo&);
    void GPolyList2PolyList(const GPolyList&,PolyList&);
    void PolyList2GPolyListClear(PolyList&,GPolyList&,const GRingInfo&);
    void GPolyList2PolyListClear(GPolyList&,PolyList&);
    RingElem homog(ConstRefRingElem the_p, const std::vector<RingElem>& the_Y); ///< hp = hom(p,y) in the ring with the y's
    void homogenized(ModuleElem& the_hv,
                     const ModuleElem& the_v,
                     const GRingInfo& theGRI);
    std::vector<long> PolyList2IndexList(const PolyList&);
    PPMonoidElem IndexList2PPMonoidElem(const PPMonoid&,
                                        const std::vector<long>&);
    std::vector<long> PPMonoidElem2IndexList(ConstRefPPMonoidElem);
    bool IsHomog(const PolyList&);
    bool IsHomog(const VectorList&);
    std::vector<degree> DegStructure(ConstRefRingElem);
    std::vector<std::vector<degree> > DegStructure(const ModuleElem&);
    std::vector<std::vector<degree> > DegStructure(const PolyList&);
    std::vector<std::vector<std::vector<degree> > > DegStructure(const VectorList&);
    PolyList MakePolyList(ConstRefRingElem);
    VectorList MakeVectorList(const ModuleElem&);
    RingElem CoeffCommonDenominator(ConstRefRingElem f);
    PolyList WithoutDenominators(const PolyList& PL, SparsePolyRing Rx);
    PolyList WithDenominator1Hom(const PolyList& PL, SparsePolyRing P);


}// end namespace cocoa



// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/TmpGReductor.H,v 1.41 2024/05/30 14:46:56 bigatti Exp $
// $Log: TmpGReductor.H,v $
// Revision 1.41  2024/05/30 14:46:56  bigatti
// Summary: Renamed myFirstPairWDeg into myIncomingWDeg, and added myUpdateIncomingWDeg
//
// Revision 1.40  2024/05/28 13:04:23  bigatti
// Summary: added field for truncation
// some rearranging in view of major cleaning
//
// Revision 1.39  2023/07/04 09:16:42  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.38  2022/02/18 14:11:50  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.37  2018/06/27 12:15:18  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.36  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.35  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.34  2017/11/29 17:48:15  bigatti
// -- added GBasisRealSolveCore
// -- added myUpdateBasisAndPairs with one arg
// -- some comments
//
// Revision 1.33  2017/11/24 17:46:40  bigatti
// -- renamed GBasisSelfSat --> GBasisSelfSatCore
// -- added GBasisSelfSat in cpkg5
//
// Revision 1.32  2017/11/23 12:37:40  bigatti
// -- added GBasisSelfSat
//
// Revision 1.31  2017/05/02 12:16:02  bigatti
// -- removed NewSparsePolyRing (now NewPolyRing does the same)
//
// Revision 1.30  2017/04/28 13:54:19  bigatti
// -- minor cleaning
// -- renamed AllGraded --> HOMOG, AllAffine --> NOGRADING
//
// Revision 1.29  2017/04/18 09:28:03  bigatti
// -- removed StatLevel member field (now using VerbosityLevel)
//
// Revision 1.28  2015/03/04 10:17:34  bigatti
// -- added (..)ElimFirst
//
// Revision 1.27  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.26  2013/10/28 13:16:47  bigatti
// -- myCreateSpecialPairs --> myCreateInputPolyPairs
//
// Revision 1.25  2013/06/12 08:51:43  bigatti
// -- added mySetMinimalGen
// -- commented out unused code
//
// Revision 1.24  2013/01/31 11:42:14  bigatti
// -- added myStats()
//
// Revision 1.23  2012/10/16 09:54:27  abbott
// Replaced  RefRingElem  by  RingElem&
// Modified many calls to  myAppendClear  (because a temporary value will
// not match a non-const reference in a fn call).
//
// Revision 1.22  2012/10/03 12:26:08  bigatti
// -- changed homogenized --> homog: now returns a RingElem (instead of void)
//
// Revision 1.21  2011/12/05 16:54:40  bigatti
// -- changed: MachineInteger --> MachineInt (just in comment)
//
// Revision 1.20  2011/03/22 11:25:49  bigatti
// -- unsigned long --> long
//
// Revision 1.19  2011/03/11 17:42:09  bigatti
// -- changed  unsigned int --> long
//
// Revision 1.18  2011/03/11 15:46:50  bigatti
// -- only whitespaces
//
// Revision 1.17  2011/03/11 15:31:58  bigatti
// -- indentation
//
// Revision 1.16  2011/03/11 14:50:39  abbott
// Changed some size_t into long.
//
// Revision 1.15  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.14  2010/07/27 07:37:13  bigatti
// -- new class GBCriteria, simplified GReductor ctor
//
// Revision 1.13  2010/07/16 09:29:58  bigatti
// -- minor cleaning and coding conventions
//
// Revision 1.12  2010/05/14 09:53:09  bigatti
// -- removed empty ctor for SugarDegree
// -- added marker for SugarDegree(uninitialized)
// -- SugarDegree for GBasis input is initialized by myPrepareGBasis
//
// Revision 1.11  2010/03/23 14:43:07  bigatti
// -- class GRingInfo estracted from TmpGPoly
//
// Revision 1.10  2009/04/27 13:16:13  bigatti
// -- added BuchbergerOpTypeMarker
//
// Revision 1.9  2009/03/20 17:37:18  bigatti
// -- added IsCopCriterionActive
//
// Revision 1.8  2009/01/30 13:41:50  bigatti
// -- enum instead of bool arguments
//
// Revision 1.7  2008/09/19 14:08:16  bigatti
// -- modified GRStats (M.Caboara)
//
// Revision 1.6  2008/07/09 16:08:40  abbott
// Added missing ampersand (for C++ const reference).
//
// Revision 1.5  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.4  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.3  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.13  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.12  2007/03/07 22:24:03  bigatti
// -- reintroduced TmpGTypes.H (waiting for a better solution)
//
// Revision 1.11  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.10  2006/12/06 17:11:16  cocoa
// -- removed #include "config.H"
//
// Revision 1.9  2006/11/24 17:20:07  cocoa
// -- reorganized includes of header files
//
// Revision 1.8  2006/11/22 15:33:39  cocoa
// -- changed handling of timings and number of reductions (M.Caboara)
//
// Revision 1.7  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.6  2006/10/06 16:46:17  cocoa
// -- syzygies for non-homogenous polynomials (Max)
// -- wip: evolution of Groebner Framework (Max)
//
// Revision 1.5  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.4  2006/08/17 09:32:26  cocoa
// -- added: flags for homogeneous input
//
// Revision 1.3  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.2  2006/07/19 07:13:00  cocoa
// -- moved "IsHomogeneous" to SparsePolyRing/FreeModule
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.16  2006/05/16 09:01:26  cocoa
// -- added "const" to arg of homogenized
// -- changed ordering matrix for modules (compatibility with GradedFreeModule)
//
// Revision 1.15  2006/05/12 17:03:16  cocoa
// -- swapped arguments in homogenized
//
// Revision 1.14  2006/05/11 16:00:22  cocoa
// -- fixed spelling of "homogenize"
//
// Revision 1.13  2006/05/04 14:25:16  cocoa
// -- major cleaning of FreeModule: created GradedFreeModule and moved
//    some code around
//
// Revision 1.12  2006/04/27 14:03:27  cocoa
// -- tidied up include files
//
// Revision 1.11  2006/04/26 09:50:40  cocoa
// -- added GReductor::ourDefaultStatLevel variable to allow CoCoAServer
//    to set statistics level
//
// Revision 1.10  2006/04/21 16:45:17  cocoa
// -- new functions by Max
//
// Revision 1.9  2006/04/11 16:42:31  cocoa
// -- added functions for Elim
//
// Revision 1.8  2006/03/13 16:56:53  cocoa
// -- changed: member data  myGRingInfoValue  is no longer a reference
//
// Revision 1.7  2006/03/02 13:45:57  cocoa
// -- changes by Max
//
// Revision 1.6  2006/02/13 14:46:45  cocoa
// -- changes by Max
//
// Revision 1.5  2006/02/13 13:45:04  cocoa
// -- changes by Max (GRingInfo)
//
// Revision 1.4  2006/01/18 15:58:20  cocoa
// -- new changes my Max
//
// Revision 1.3  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.2  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.5  2004/10/29 16:02:52  cocoa
// -- added field myDivMaskImplPtr for creating GPolys with the same
//    DivMask::base
// -- calls to LCMwMask (instead of LCM) for fast divisibility tests
// -- function IsDivisible had wrong semantics --> swapped arguments everywhere
//
// Revision 1.4  2004/06/16 16:13:41  cocoa
// Improved I/O facilities with knock-on changes
//
// Revision 1.3  2004/03/04 11:37:17  cocoa
// -- updated code for Borel reductors:
//    ReductorData fields: myGPoly->myGPolyPtr;  new: myCount, IamBorelUpdated
//    myBorelReductors is now in Reductors (was in GReductor)
//    Reductors: field: IhaveBorelReductors;  type: enum UseBorelMarker
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.11  2003/06/23 17:09:54  abbott
// Minor cleaning prior to public release.
// Improved the include directives,
//
// Revision 1.10  2003/05/29 16:49:58  bigatti
// - change: myRingSpecialIndex is now an int
// - added:  flag for criterion (to be disables for non-comm Groebner)
//
// Revision 1.9  2003/05/28 14:23:09  bigatti
// - new code for modules
//
// Revision 1.8  2003/05/14 17:06:07  bigatti
// - new functions for "BorelReductorsPolys" and saturating algorithm
//
// Revision 1.7  2002/09/19 17:25:28  bigatti
// - Cleaner code based on PolyRing
//
// Revision 1.6  2002/05/13 11:47:34  bigatti
// - new data structure for "Reductors"
//
// Revision 1.5  2002/04/15 17:17:30  bigatti
// - Max's new code
//
// Revision 1.4  2002/04/09 14:21:57  bigatti
// - CKR criterion
//
// Revision 1.3  2002/01/31 17:29:47  bigatti
// - new function: Stampa_Reductors_Stats
//
// Revision 1.2  2001/12/12 18:19:55  bigatti
// - new structure of reduction
//
// Revision 1.1  2001/12/05 13:27:48  bigatti
// Initial revision
//

#endif
