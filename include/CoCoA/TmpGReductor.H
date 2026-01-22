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

    enum BuchbergerOpTypeFlag { HomogeneousAlg, SaturatingAlg, AffineAlg };

    void myCtorAux(const BuchbergerOpTypeFlag theBuchbergerOpType);
    
    GReductor(const GRingInfo&,
              const PolyList&,
              const BuchbergerOpTypeFlag theBuchbergerOpType=HomogeneousAlg,
              const GBCriteria criteria=GBCriteria(GBCriteria::AllSet));
    GReductor(const GRingInfo&,
              const GPolyList&,
              const BuchbergerOpTypeFlag theBuchbergerOpType=HomogeneousAlg,
              const GBCriteria criteria=GBCriteria(GBCriteria::AllSet));
    GReductor(const GReductor&);// copy ctor not working, to be fixed
    ~GReductor(){};

    // GReductor(const GRingInfo&,
    //           GPolyList&,
    //           const ClearMarker,
    //           const BuchbergerOpTypeFlag theBuchbergerOpType=HomogeneousAlg,
    //           const GBCriteria criteria=GBCriteria(GBCriteria::AllSet));

    BuchbergerOpTypeFlag myGetBuchbergerOpType() const {return myBuchbergerOpType;}
    bool myPreparationDone() const {return myPrepared;}
    //    bool WrongLPPFound() const {return myWrongLPPFoundValue;}
    const Stats& myStats() const {return myStat;}

    long myReductorsLen() const {return myTrueReductors.myLen();}
    long myGBasisLen() const {return len(myGB);}
    long myPairsLen() const {return len(myPairs);}

    void mySetBuchbergerOpType(const BuchbergerOpTypeFlag theBuchbergerOpType) 
        {myBuchbergerOpType=theBuchbergerOpType;}

    std::vector<RingElem> myExportGBasis();
    std::vector<RingElem> myExportMinGens();
    void myCopyGBasis_GPoly(GPolyList&);// GB output in GPoly form
    void myGBasisClear(GPolyList&); ///< GB output in GPoly form
    void myCopyGBasis_module(VectorList& outGBasis);
    void myCopyMinGens_module(VectorList& outMinGens);
    void myDoAFFGBasis(); // I am working on sugar selection strategies
    void myDoGBasis();
    void myDoGBasisTEST();
    RingElem myDoGBasisElimFirst(ConstRefPPMonoidElem ElimIndsProd);
    void myReduceCurrentSPoly();
    void myPrepareGBasis();
    void myFinalizeGBasis();//Last operations (stats at the moment) immediately before shutting down a GB computation
    void myDoGBasis(const int ReductionNumber); // Performs ReductionNumber reductions, -1 means unlimited
    void myReduceUntilNonZeroRedSPoly(); // Reduces until an SPoly does not reduces to zero. Stop immediately when that happens (no updates)
    // void myReduceUntilWrongLPPFound(RefPPMonoidElem,
    //                                 std::vector<RingElem>& );
//A24    GPoly GetSPoly()const{return mySPoly;};// Reads the SPoly
//A24    void GetCandidateGBasis(PolyList&)const;// the polys computed up to now and the non processed generators
    long myAge() const {return myAgeValue;}
    long myTruncDeg() const {return myTruncDegValue;}
//A24    long GetNReductions() const {return myStat.myNReductions;};
    void mySetTruncDeg(long TruncDeg) { myTruncDegValue = TruncDeg; }
    static constexpr long ourNoTruncValue = -1;  // no GB truncation
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
    //    bool IsDynamicAlgorithm;
    //    bool myWrongLPPFoundValue; // DYNAMIC ALGORITHM
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
    SparsePolyRing MakeElimRing(const SparsePolyRing& P,
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


}// end namespace cocoa

#endif
