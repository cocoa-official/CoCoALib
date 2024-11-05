//   Copyright (c)  2005-2017 John Abbott and Anna M. Bigatti
//   Author:  2005  Massimo Caboara, 2010-2023 Anna M. Bigatti

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

#include "CoCoA/TmpGOperations.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseMatrix.H" // for NewDenseMat
#include "CoCoA/FractionField.H"
#include "CoCoA/FreeModule.H"
#include "CoCoA/MatrixForOrdering.H" // for MakeTermOrdMat
#include "CoCoA/MatrixView.H" // for ConcatVer
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingQQ.H" // for IsQQ in ComputeSaturationByPrincipal
#include "CoCoA/RingZZ.H" // for RingZZ
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H" // for GradingMat
#include "CoCoA/VectorOps.H"
#include "CoCoA/factor.H" // for factor in ComputeSaturationByPrincipal
#include "CoCoA/matrix.H" // for ConstMatrixView
#include "CoCoA/symbol.H"
#include "CoCoA/verbose.H"
// #include "CoCoA/TmpGReductor.H"  // included in TmpGOperations.H 

#include <algorithm>
using std::stable_sort;
#include <list>
using std::list;
#include <functional>
using std::less;
using std::binary_function;
using std::mem_fun_ref; // for calling GPair::complete on GPairList
#include <utility>
using std::make_pair;
#include <iostream>
using std::ostream;
using std::endl;
using std::flush;

namespace CoCoA
{


  namespace // anonymous
  {
    bool IsHomogGrD0(const PolyList& L)
    {
      if (L.empty()) return true;
      if (GradingDim(owner(L[0]))==0) return true;
      return IsHomog(L);
    }

    bool IsHomogGrD0(const VectorList& L)
    {
      if (L.empty()) return true;
      if (GradingDim(RingOf(owner(L[0])))==0) return true;
      return IsHomog(L);
    }

    bool IsHomogGrD0(ConstRefRingElem f)
    {
      if (GradingDim(owner(f))==0) return true;
      return IsHomog(f);
    }

  }
      

// GBasis

  void ComputeGBasis(VectorList& outGB, VectorList& outMinGens, const VectorList& inGens, const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    const FreeModule FM=owner(inGens);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM));
    const SparsePolyRing OldP(RingOf(FM));
    if (!IsField(CoeffRing(OldP)))
      CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool IsSatAlg=false;
    const GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(inGens),
                        IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=EmbedVectorList(inGens,GRI);
    if (EmbeddedPolys.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    VectorList TmpGB;
    VectorList TmpMinGens;
    GBR.myCopyGBasis(TmpGB);
    if (GradingDim(OldP)>0 && IsHomog(inGens)) GBR.myCopyMinGens(TmpMinGens);
    outGB.clear(); // just to remember to clean this up
    outMinGens.clear(); // just to remember to clean this up
    //    outGB =      DeEmbedPolyList(TmpGB, GRI);
    //    outMinGens = DeEmbedPolyList(TmpMinGens, GRI);
    outGB =      TmpGB;
    outMinGens = TmpMinGens;
  }//ComputeGBasis


  void ComputeGBasis(PolyList& outGB, PolyList& outMinGens, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    SparsePolyRing SPR(owner(inGens));
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool IsSatAlg=false;
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(inGens, Rx));
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList TmpGB;
      PolyList TmpMinGens;
      GBR.myCopyGBasis(TmpGB);
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myCopyMinGens(TmpMinGens);
      outGB.clear();//just to remember to clean this up
      outMinGens.clear();//just to remember to clean this up
      outGB = WithDenominator1Hom(TmpGB, SPR);
      monic(outGB);  //  2016-11-22: make monic
      outMinGens = WithDenominator1Hom(TmpMinGens, SPR);
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, inGens);
      GBR.myDoGBasis();// homog input standard alg interred
      GBR.myCopyGBasis(outGB);
      outMinGens.clear();//just to remember to clean this up
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myCopyMinGens(outMinGens);
    }          
  }//ComputeGBasis
  

  namespace { // anonymous
    bool NoWDegGt(const PolyList& F, long D)
    {
      for (auto& f:F)  if (wdeg(f)[0] > D) return false;
      return true;
    }
  }
  
    
  void ComputeGBasisTrunc(PolyList& outGB, PolyList& outMinGens, long& TruncDeg, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    SparsePolyRing SPR(owner(inGens));
    CoCoA_ASSERT_ALWAYS(TruncDeg > 0);
    if (!IsField(CoeffRing(SPR)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (GradingDim(SPR)!=1)  CoCoA_THROW_ERROR1(ERR::ReqGradingDim1);
    if (!IsHomog(inGens))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    bool IsSatAlg=false;
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(inGens, Rx));
      GBR.mySetTruncDeg(TruncDeg);
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList TmpGB;
      PolyList TmpMinGens;
      GBR.myCopyGBasis(TmpGB);
      if (NoWDegGt(inGens, TruncDeg)) GBR.myCopyMinGens(TmpMinGens);
      outGB = WithDenominator1Hom(TmpGB, SPR);
      monic(outGB);  //  2016-11-22: make monic
      outMinGens = WithDenominator1Hom(TmpMinGens, SPR);
      if (GBR.myTruncDeg() == 0)  TruncDeg = 0;
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, inGens);
      GBR.mySetTruncDeg(TruncDeg);
      GBR.myDoGBasis();// homog input standard alg interred
      GBR.myCopyGBasis(outGB);
      outMinGens.clear();//just to remember to clean this up
      if (NoWDegGt(inGens, TruncDeg)) GBR.myCopyMinGens(outMinGens);
      if (GBR.myTruncDeg() == 0)  TruncDeg = 0;
    }
  }//ComputeGBasisTrunc
  
  
  /*
  void ComputeSATMixGBasis(PolyList& theGB, const PolyList& thePL)
  {      
    if (thePL.empty())
    {
      theGB.clear();
      return;
    }
    bool IsSatAlg=false;
    SparsePolyRing SPR(owner(thePL));
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR))) CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(thePL, Rx));
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList TmpGB;
      GBR.myGBasis(TmpGB);
      theGB.clear();//just to remember to clean this up
      theGB=WithDenominator1Hom(TmpGB, SPR);
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(thePL),IsSatAlg,NewDivMaskEvenPowers());
      GReductor GBR(GRI, thePL);
      GBR.myDoGBasis();// homog input standard alg interred
      GBR.myGBasis(theGB);
    }
  }//ComputeSATMixGBasis
  */

  void ComputeGBasisSelfSatCore(PolyList& outGB, const PolyList& inPL, const CpuTimeLimit& CheckForTimeout)
  {
    if (inPL.empty())
    {
      outGB.clear();
      return;
    }
    //    bool IsSatAlg=true;
    SparsePolyRing SPR(owner(inPL));
    PolyList TmpGB;
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(inPL), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(inPL, Rx),
                    GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      GBR.myCopyGBasis(TmpGB);
      TmpGB = WithDenominator1Hom(TmpGB, SPR);
    }
    else
    {
      GRingInfo GRI(SPR, IsHomogGrD0(inPL), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, inPL, GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      GBR.myCopyGBasis(TmpGB);
    }
    swap(outGB, TmpGB);
  } // ComputeGBasisSelfSatCore


  void ComputeGBasisRealSolve(PolyList& outGB, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGB.clear();
      return;
    }
    SparsePolyRing SPR(owner(inGens));
    if (!IsField(CoeffRing(SPR)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (characteristic(SPR)!=0)  CoCoA_THROW_ERROR1(ERR::ReqChar0);
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(inGens),false/*IsSatAlg*/,NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(inGens, Rx));
      GBR.myDoGBasisRealSolve();// homog input standard alg interred
      PolyList TmpGB;
      PolyList TmpMinGens;
      GBR.myCopyGBasis(TmpGB);
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myCopyMinGens(TmpMinGens);
      outGB.clear();//just to remember to clean this up
      outGB = WithDenominator1Hom(TmpGB, SPR);
      monic(outGB);  //  2016-11-22: make monic
    }
    else
    {
      GRingInfo GRI(SPR,IsHomogGrD0(inGens),false/*IsSatAlg*/,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, inGens);
      GBR.myDoGBasisRealSolve();
      GBR.myCopyGBasis(outGB);
    }          
  } // ComputeGBasisRealSolve


/*
void ComputeGBasisFrameWork(PolyList& theGB, const PolyList& theInputPolyList)

 {
   theGB.clear();
   if (theInputPolyList.empty())
     return;
   const SparsePolyRing SPR(owner(theInputPolyList));
   GRingInfo GRI(SPR,IsHomogGrD0(theInputPolyList),NewDivMaskEvenPowers());
   GReductor GBR(GRI, theInputPolyList,GReductor::ourDefaultStatLevel);


   if (false)
   {
     GBR.myDoGBasis();// homog input standard alg interred
     GBR.myGBasis(theGB);
   }


   if (false)
   {
     GBR.myPrepareGBasis();
     while (GBR.myPairsLen()!=0)
     {
       GBR.myDoGBasis(1);
     }
     GBR.myGBasis(theGB);
   }

  if (false)
  {
   //non null Spoly by non null Spoly
   GBR.myPrepareGBasis();
   //GPoly SP(GRI);
   while (GBR.myPairsLen()!=0)
   {
     GBR.myReduceUntilNonZeroRedSP();
     //SP=GBR.GetSPoly();
     //.......Things.......
     //PolyList PL;
     //GBR.GetCandidateGBasis(PL);
     //GBR.SetSPoly(SP);
     if (GBR.myPairsLen()!=0)
     {
       myTrueReductors.interreduce(mySPoly);
       myUpdateBasisAndPairs();
       if (myCurrentPairDeg!=myOldDeg&&myTrueReductors.IhaveBorelReductors())
         myTrueReductors.myBorelReductorsUpdateInNextDegree(myCurrentPairDeg);
     }
    }
   GBR.myCopyGBasis(theGB);
  }

 }//ComputeGBasisFrameWork
*/



  void ComputeElim(PolyList& theElimResult, const PolyList& thePL, ConstRefPPMonoidElem inds,
                   const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeElim");
    VERBOSE(99) << "-- called --" << std::endl;
    if (thePL.empty())
    {
      theElimResult.clear();
      return;
    }
    PolyList ElimResult;
    const SparsePolyRing OldSPR=owner(thePL);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    bool IsHomogGrD0PL = false;
    if  (GradingDim(OldSPR)!=0 && IsHomogGrD0(thePL)) IsHomogGrD0PL = true;
    SparsePolyRing NewSPR=MakeElimRingFromOld(OldSPR,IndexList, IsHomogGrD0PL);
    RingHom OldToNew = PolyAlgebraHom(OldSPR, NewSPR, indets(NewSPR));
    RingHom NewToOld = PolyAlgebraHom(NewSPR, OldSPR, indets(OldSPR));
    PolyList NewPL;
    for(const RingElem& g: thePL)
      NewPL.push_back(OldToNew(g));
///    for (PolyList::const_iterator it=thePL.begin();it!=thePL.end();++it)
///      NewPL.push_back(OldToNew(*it));
    PPMonoidElem ElimIndsProd(LPP(OldToNew(monomial(OldSPR,inds))));
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis(GB, MinGens, NewPL, CheckForTimeout);
    for (const RingElem& g: GB)
    {
      if ( IsConstant(g) ) // redmine #1647
      {
        ElimResult.clear();
        VERBOSE(99) << "g is constant" << std::endl; /////////////////
        ElimResult.push_back(OldSPR->myOne());
        break;
      }
      if (IsCoprime(LPP(g), ElimIndsProd))
        ElimResult.push_back(NewToOld(g));
///    for (PolyList::iterator it=GB.begin();it!=GB.end();++it)
///      if (IsCoprime(LPP(*it), ElimIndsProd))
///        ElimResult.push_back(NewToOld(*it));
    }
    swap(theElimResult,ElimResult);
    return;
  }//ComputeElim


  RingElem ComputeElimFirst(const PolyList& inPL, ConstRefPPMonoidElem inds,
                            const CpuTimeLimit& CheckForTimeout)
  {
    const SparsePolyRing OldSPR=owner(inPL);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    bool IsHomogGrD0PL = (GradingDim(OldSPR)==0 ? false : IsHomogGrD0(inPL));
    if (!IsHomogGrD0PL)  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    bool IsSatAlg = false;
    SparsePolyRing NewSPR=MakeElimRingFromOld(OldSPR,IndexList, IsHomogGrD0PL);
    RingHom OldToNew = PolyAlgebraHom(OldSPR, NewSPR, indets(NewSPR));
    RingHom NewToOld = PolyAlgebraHom(NewSPR, OldSPR, indets(OldSPR));
    PPMonoidElem ElimIndsProd = LPP(OldToNew(monomial(OldSPR,inds)));
    if (!IsField(CoeffRing(OldSPR)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(OldSPR)))
      CoCoA_THROW_ERROR2(ERR::NYI, "Only for FFp");
    else
    {
      GRingInfo GRI(NewSPR,IsHomogGrD0(inPL),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, OldToNew(inPL));
      return NewToOld(GBR.myDoGBasisElimFirst(ElimIndsProd));
    }
    return zero(OldSPR); // just to keep the compiler quiet
  }


  void ComputeElim(VectorList& /*theElimResult*/,
                   const VectorList& /*theVL*/,
                   ConstRefPPMonoidElem /*inds*/,
                   const CpuTimeLimit& /*CheckForTimeout*/)
  {
    // Remember to check if ordering is TOPOS ordering here
    // if not, use hom
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeElim


  void ComputeSyz(VectorList& theSyzResult, const FreeModule& SyzFM, const VectorList& theVL,
            const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsHomogGrD0(theVL))
      CoCoA_THROW_ERROR2(ERR::NYI, "Not Yet Tested for non-homog input");
    if (theVL.empty())
    {
      theSyzResult.clear();
      return;
    }
    bool IsSatAlg=false;
    VectorList SyzResult;
    ModOrdTypeForcing MOType;
    if (IsHomogGrD0(theVL))
      MOType=WDegPosTO;
    else
      MOType=PosWDegTO;
    const FreeModule FM=owner(theVL);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM,MOType));
    const SparsePolyRing OldP(RingOf(FM));
    // Note: the GRI should build itself SyzFM and NewP from the data and deduce FM and OldP.
    //       All the embedding/deembedding functions should be members of GRI.
    GRingInfo GRI(NewP,OldP,FM,SyzFM,IsHomogGrD0(theVL),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=SyzEmbedVectorList(theVL,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    SyzResult=DeEmbedPolyList(OutputPolyList,GRI,NumCompts(FM));
    swap(theSyzResult,SyzResult);
    return;
  }//ComputeSyz


  void ComputeSyz(VectorList& theSyzResult, const FreeModule& SyzFM, const PolyList& inPL,
            const CpuTimeLimit& CheckForTimeout)
  {
    //if (!IsHomogGrD0(inPL))
    //  CoCoA_THROW_ERROR2("Not Yet Tested for non-homogeneous input", "ComputeSyz(const VectorList&)");
    if (inPL.empty())
    {
      theSyzResult.clear();
      return;
    }
    bool IsSatAlg=false;
    VectorList SyzResult;
    ModOrdTypeForcing MOType;
    if (IsHomogGrD0(inPL))
      MOType=WDegPosTO;
    else
      MOType=PosWDegTO;
    const SparsePolyRing OldP(owner(inPL));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbedding(OldP,MOType));
    GRingInfo GRI(NewP,OldP,SyzFM,IsHomogGrD0(inPL),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=SyzEmbedPolyList(inPL,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    SyzResult=DeEmbedPolyList(OutputPolyList,GRI,1);
    swap(theSyzResult,SyzResult);
    return;
  }//ComputeSyz


  void ComputeIntersection(VectorList& theIntersectionResult,
                           const VectorList& theVL1,
                           const VectorList& theVL2,
                           const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    VectorList IntersectionResult;
    const FreeModule FM=owner(theVL1);
    //const SparsePolyRing NewP(MakeNewPRingFromModule(FM,WDegPosTO));
    const SparsePolyRing NewP(MakeNewPRingFromModulePosFirst(FM,IsHomogGrD0(theVL1)&&IsHomogGrD0(theVL2)));
    const SparsePolyRing OldP(RingOf(FM));
    GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(theVL1)&&IsHomogGrD0(theVL2),
                  IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=IntEmbedVectorLists(theVL1,theVL2,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    IntersectionResult=DeEmbedPolyList(OutputPolyList,GRI,NumCompts(FM));
    swap(IntersectionResult,theIntersectionResult);
    return;
  }//ComputeIntersection


  void ComputeIntersection(PolyList& theIntersectionResult,
                           const PolyList& inPL1,
                           const PolyList& inPL2,
                           const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    PolyList IntersectionResult;
    const SparsePolyRing OldP(owner(inPL1));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbeddingPosFirst(OldP,IsHomogGrD0(inPL1)&&IsHomogGrD0(inPL2)));
    GRingInfo GRI(NewP,OldP,IsHomogGrD0(inPL1)&&IsHomogGrD0(inPL2),IsSatAlg,
                  NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=IntEmbedPolyLists(inPL1,inPL2,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    IntersectionResult=DeEmbedPolyListToPL(OutputPolyList,GRI,1);
    swap(theIntersectionResult,IntersectionResult);
    return;
  }//ComputeIntersection


// Colon by a single vector
  void ComputeColonByPrincipal(PolyList& theColonResult, const VectorList& theVL1, const VectorList& theVL2,
            const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    //    for (VectorList::const_iterator it=theVL1.begin();it!=theVL1.end();++it)
    //      if (!IsHomogGrD0(*it))
      if (!IsHomogGrD0(theVL1))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
      //    for (VectorList::const_iterator it=theVL2.begin();it!=theVL2.end();++it)
      //      if (!IsHomogGrD0(*it))
      if (!IsHomogGrD0(theVL2))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
    PolyList ColonResult;
    const FreeModule FM=owner(theVL1);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM,WDegPosTO));
    const SparsePolyRing OldP(RingOf(FM));
    GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(theVL1)&&IsHomogGrD0(theVL2),
                  IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=ColonEmbedVectorLists(theVL1,theVL2,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    ColonResult=DeEmbedPolyListToPL(OutputPolyList,GRI,NumCompts(FM));
    swap(ColonResult,theColonResult);
    return;
  }//ComputeColonByPrincipal

// Colon module:module
  void ComputeColon(PolyList& theResult,
                    const VectorList& theVL1,
                    const VectorList& theVL2,
                    const CpuTimeLimit& CheckForTimeout)
  {
    if (theVL1.empty() && theVL2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    PolyList ColonResult;
    if (theVL2.empty())
    {
      ColonResult.push_back(one(RingOf(owner(theVL1))));
      swap(ColonResult, theResult);
      return;
    }
    if (theVL1.empty())
    {
      swap(ColonResult,theResult);
      return;
    }
    PolyList tmp;
    VectorList::const_iterator it=theVL2.begin();
    ComputeColonByPrincipal(ColonResult,theVL1,MakeVectorList(*it), CheckForTimeout);
    ++it;
    for (;it!=theVL2.end();++it)
    {
      ComputeColonByPrincipal(tmp,theVL1,MakeVectorList(*it), CheckForTimeout);
      ComputeIntersection(ColonResult,ColonResult,tmp, CheckForTimeout);
    }
    swap(ColonResult,theResult);
    return;
  }//ComputeColon


  void ComputeColonByPrincipal(PolyList& theResult,
                               const PolyList& inPL1,
                               ConstRefRingElem f,
                               const CpuTimeLimit& CheckForTimeout)
  {
    PolyList TmpColonResult;
    if (IsZero(f))
    {
      TmpColonResult.push_back(one(owner(f)));
      swap(theResult, TmpColonResult);
      return;
    }
    if (inPL1.empty())
    {
      theResult.clear();
      return;
    }
    bool IsSatAlg=false;
    const SparsePolyRing OldP(owner(inPL1));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbeddingPosFirst(OldP,IsHomogGrD0(inPL1)&&IsHomogGrD0(f)));
    GRingInfo GRI(NewP,OldP,IsHomogGrD0(inPL1)&&IsHomogGrD0(f),
                  IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=ColonEmbedPolyLists(inPL1, std::vector<RingElem>(1,f), GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    TmpColonResult = DeEmbedPolyListToPL(OutputPolyList, GRI, 1);
    swap(theResult, TmpColonResult);
    return;
  }//ComputeColonByPrincipal


// Colon ideal:ideal
  void ComputeColon(PolyList& theResult,
                    const PolyList& inPL1,
                    const PolyList& inPL2,
                    const CpuTimeLimit& CheckForTimeout)
  {
    if (inPL1.empty() && inPL2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    PolyList TmpColonResult;
    if (inPL2.empty())
    {
      TmpColonResult.push_back(one(owner(inPL1)));
      swap(theResult, TmpColonResult);
      return;
    }
    if (inPL1.empty()) // ComputeColonByPrincipal checks this again
    {
      theResult.clear();
      return;
    }
    PolyList tmp;
    PolyList::const_iterator it=inPL2.begin();
    ComputeColonByPrincipal(TmpColonResult, inPL1, *it, CheckForTimeout);
    for (++it; it!=inPL2.end(); ++it)
    {
      ComputeColonByPrincipal(tmp, inPL1, *it, CheckForTimeout);
      ComputeIntersection(TmpColonResult, TmpColonResult, tmp, CheckForTimeout);
      if (TmpColonResult.empty()) break;
    }
    swap(theResult, TmpColonResult);
    return;
  }//ComputeColon


  void ComputeColonByPrincipal(VectorList& /*theResult*/,
                               const VectorList& /*theVL*/,
                               const PolyList& /*inPL*/,
                               const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeColonByPrincipal
  

  namespace // anonymous namespace for file local functions and definitions
  {
  // These procedures are for EqualLPPs used by Saturation
  
  // PolyList::const_iterators are ordered according to LPP of their polys
  bool ByLPP(const PolyList::const_iterator& it,
             const PolyList::const_iterator& it1)
  {
    return (SparsePolyRingPtr(owner(*it))->myCmpLPP(raw(*it),raw(*it1)) == -1);
  }


  bool AreEqualLPPs(std::vector<PolyList::const_iterator>& theV1,
                    std::vector<PolyList::const_iterator>& theV2)
  {
    const long lenV1 = len(theV1);
    const long lenV2 = len(theV2);
    if (lenV1 != lenV2)  return false;
    stable_sort(theV1.begin(), theV1.end(), ByLPP);
    stable_sort(theV2.begin(), theV2.end(), ByLPP);
    for (int i=0; i!=lenV1; ++i)
      if (LPP(*(theV1[i])) != LPP(*(theV2[i])))  return false;
    return true;
  }
    
    
//     bool AreEqualLPPs(const std::vector<VectorList::const_iterator>& theV1,
//                       const std::vector<VectorList::const_iterator>& theV2)
//     {
    // const long lenV1 = len(V1);
    // const long lenV2 = len(V2);
//       if (lenV1 != lenV2)  return false;
//       stable_sort(theV1.begin(), theV1.end(), ByLPP);
//       stable_sort(theV2.begin(), theV2.end(), ByLPP);
//       for (int i=0; i!=lenV1; ++i)
//         //if (*(theV1[i])!=*(theV2[i]))
//         // when LT exists if (LT(*(theV1[i]))!=LT(*(theV2[i])))
//         return false;
//       return true;
//     }

    // Useful when you have I\subset J and you want to check I==J
    bool AreEqualLPPs(const PolyList& inPL1, const PolyList& inPL2)
    {
      std::vector<PolyList::const_iterator> V1,V2;
      for (PolyList::const_iterator it=inPL1.begin(); it!=inPL1.end(); ++it)
        V1.push_back(it);
      for (PolyList::const_iterator it=inPL2.begin(); it!=inPL2.end(); ++it)
        V2.push_back(it);
      return AreEqualLPPs(V1,V2);
    }
    
//     // Useful when you have M\subset N and you want to check M==N
//     bool AreEqualLPPs(const VectorList& theVL1, const VectorList& theVL2)
//     {
//       std::vector<VectorList::const_iterator> V1,V2;
//       for (VectorList::const_iterator it=theVL1.begin();it!=theVL1.end();++it)
//         V1.push_back(it);
//       for (VectorList::const_iterator it=theVL2.begin();it!=theVL2.end();++it)
//         V2.push_back(it);
//       return AreEqualLPPs(V1,V2);
//     }
  
  } // anonymous namespace

   
//---------------------------------------------------------


  void ComputeSaturationByIrred(/*const*/ PolyList& TmpPL,
                                ConstRefRingElem f,
                                const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturationByIrred");
    VERBOSE(99) << "-- called --" << std::endl;
    const RingElem h = (indets(owner(f))).back();
    TmpPL.push_back(h*f-1);
    ComputeElim(TmpPL, TmpPL, LPP(h), CheckForTimeout);
  } // ComputeSaturationByIrred

  
  namespace { // anonymous

    RingElem SatByIndet(const RingElem& f, long IndetIndex)
    {
      const ring& P = owner(f);
      long pwr = exponent(LPP(f), IndetIndex);
      if (pwr == 0) return f;
      const PPMonoidElem t = power(indet(PPM(P),IndetIndex), pwr);
      RingElem ans(P);
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        PushBack(ans, coeff(it), PP(it)/t);
      }
      return ans;
    }

  }
  

  void ComputeSaturationHomogByPP(PolyList& outPL,
                                  const PolyList& inPL,
                                  ConstRefPPMonoidElem t,
                                  const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturationHomogByPP");
    VERBOSE(99) << "-- called --" << std::endl;
//    if (IsZero(I))  return ideal(zero(Ph));
    const ring& P = owner(inPL);
    const ring& K = CoeffRing(P); // CoCoA_ASSERT(IsField(k));
    const long GrDim = GradingDim(P);
    const long NumX = NumIndets(P);
    std::vector<symbol> names;  for (long i=0; i<NumX; ++i) { names.push_back(IndetSymbol(P,i)); }
    VERBOSE(90) << "GrDim = " << GrDim << "   W = " << GradingMat(P) << std::endl;
    matrix WRL0 = NewDenseMat(ConcatVer(GradingMat(P), ZeroMat(RingZZ(), 1, NumIndets(P))));
    PolyList tmpPL = inPL;
    PolyList MinGens; // unused
    std::vector<long> expv = exponents(t);
    for (long i=0; i<NumIndets(P); ++i)
      if (expv[i] != 0)
      {
        VERBOSE(95) << "doing indet " << i << " = " << indet(P,i) << std::endl;
        SetEntry(WRL0, GrDim,i, -1);
        ring PzDegRevLex = NewPolyRing(K, names, NewMatrixOrdering(MakeTermOrdMat(WRL0), GrDim));
        SetEntry(WRL0, GrDim,i, 0);
        RingHom phi = PolyAlgebraHom(owner(tmpPL[0]), PzDegRevLex, indets(PzDegRevLex));
        ComputeGBasis(tmpPL, MinGens, phi(tmpPL), CheckForTimeout);
        for (RingElem& g: tmpPL)  g = SatByIndet(g,i);
      }
    RingHom phi = PolyAlgebraHom(owner(tmpPL[0]), P, indets(P));
//    tmpPL = phi(tmpPL);
    for (auto& g:tmpPL)
    {
      if ( IsZero(g) )  CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
      if ( IsConstant(g) ) // redmine #1647
      {
        tmpPL.clear();
        VERBOSE(99) << "g is constant" << std::endl; /////////////////
        tmpPL.push_back(one(P));
        break;
      }
      g = phi(g);
    }
    swap(outPL, tmpPL);
  } // ComputeSaturationHomogByPP


  void ComputeSaturationByPrincipal(PolyList& outPL,
                                    const PolyList& inPL,
                                    ConstRefRingElem f,
                                    const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturationByPrincipal");
    VERBOSE(99) << "-- called --" << std::endl;
    // non-homogeneous
    if (IsInvertible(f))
    {
      PolyList TmpPL(inPL);
      swap(outPL, TmpPL);
      return;
    }
    if (IsMonomial(f) /*IsPP*/ && IsHomog(inPL))
    {
      ComputeSaturationHomogByPP(outPL, inPL, LPP(f), CheckForTimeout);
      return;
    }
    const SparsePolyRing P = owner(inPL);
    const SparsePolyRing NewP = NewPolyRing(CoeffRing(P),
                                            NewSymbols(NumIndets(P)+1));
    const std::vector<RingElem>& x_new = indets(NewP);
    RingHom PToNewP = PolyAlgebraHom(P, NewP, std::vector<RingElem>(x_new.begin(), x_new.end()-1));
    std::vector<RingElem> x = indets(P);
    x.push_back(one(P));
    RingHom NewPToP = PolyAlgebraHom(NewP, P, x);
    PolyList TmpPL = PToNewP(inPL);
    // 2023-05-02: JAA: would be better to detect & handle specially when f is monomial!!
    // If f is not monomial & we can factorize, do a succession of saturations...
    //    if (!IsMonomial(f) && IsQQ(CoeffRing(P)))
    if (IsQQ(CoeffRing(P)))
    {
      const std::vector<RingElem> F = factor(f).myFactors();//const factorization<RingElem> F=factor(f);
      ComputeSaturationByIrred(TmpPL, PToNewP(F[0]), CheckForTimeout);
      for (long i=1; i<len(F); ++i)
        ComputeSaturationByIrred(TmpPL, PToNewP(F[i]), CheckForTimeout);
    }
    else // 2023-05-02: JAA: would be better to use radical(f) below.
      ComputeSaturationByIrred(TmpPL, PToNewP(f), CheckForTimeout);
    TmpPL = NewPToP(TmpPL);
    swap(outPL, TmpPL);
  } // ComputeSaturationByPrincipal


  void ComputeSaturation(PolyList& theSaturationResult,
                         const PolyList& inPL1,
                         const PolyList& inPL2,
                         const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturation");
    VERBOSE(99) << "-- called --" << std::endl;
    if (inPL1.empty() && inPL2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    if (inPL2.empty())
    {
      theSaturationResult.clear();// this or swap? this look better
      theSaturationResult.push_back(one(owner(inPL1)));
      return;
    }
    if (inPL1.empty())
    {
      theSaturationResult.clear();
      return;
    }

    PolyList TmpPL2;
    if (len(inPL2)==1)
      ComputeSaturationByPrincipal(TmpPL2, inPL1, inPL2.front(), CheckForTimeout);
    else
    {
      PolyList TmpPL1;
      ComputeColon(TmpPL1, inPL1, inPL2, CheckForTimeout);
      ComputeColon(TmpPL2, TmpPL1, inPL2, CheckForTimeout);
      while (!AreEqualLPPs(TmpPL1, TmpPL2))
      {
        swap(TmpPL1,TmpPL2);
        ComputeColon(TmpPL2, TmpPL1, inPL2, CheckForTimeout);
      }
    }
    swap(theSaturationResult,TmpPL2);
  }//ComputeSaturation


  void ComputeSaturationByPrincipal(VectorList& /*theSaturation*/,
                                    const VectorList& /*theVL*/,
                                    const PolyList& /*inPL*/,
                                    const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeSaturationByPrincipal


  void ComputeSaturation(VectorList& /*theSaturation*/,
                         const VectorList& /*theVL*/,
                         const PolyList& /*inPL*/,
                         const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeSaturation


  void ComputeHomogenization(VectorList& /*theHomogResult*/,
                             const VectorList& /*theVL*/,
                             const PolyList& /*inPL*/,
                             const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeHomogenization


  void ComputeHomogenization(PolyList& theHomogResult,
                             const PolyList& inPL1,
                             const PolyList& theIndets,
                             const CpuTimeLimit& CheckForTimeout)
  {
    if (inPL1.empty())
    {
      theHomogResult.clear();
      return;
    }
    PolyList HomogResult;
    const SparsePolyRing SPR=owner(inPL1);
    RingElem IndetProduct(product(theIndets));
    RingElem tmp(SPR);
    for (const RingElem& g: inPL1)
      HomogResult.push_back(homog(g, theIndets));
//    for (PolyList::const_iterator it=inPL1.begin();it!=inPL1.end();++it)
//      HomogResult.push_back(homog(*it, theIndets));
    ComputeSaturationByPrincipal(HomogResult,HomogResult,IndetProduct, CheckForTimeout);
    swap(theHomogResult,HomogResult);
    return;
  }//ComputeHomogenization


// WARN: it supposes ComputeSaturationByPrincipal returns a GB
  bool RadicalMembership(const PolyList& PL, ConstRefRingElem the_f,
                         const CpuTimeLimit& CheckForTimeout)
  {
    PolyList PL1;
    ComputeSaturationByPrincipal(PL1, PL, the_f, CheckForTimeout);
    if (len(PL) != 1) return false;
    //    monic(PL1);
    //    if (IsOne(PL1.front()))
    //      return true;
    //    else
    //      return false;
    return IsInvertible(PL1.front());
  }//RadicalMembership

  void ComputeLT(VectorList& /*theLTResult*/,
                 const VectorList& /*theVL*/,
                 const CpuTimeLimit& /*CheckForTimeout*/)
  {
    // Waiting for LT of a ModuleElem
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeLT


  void ComputeLT(PolyList& theLTResult, const PolyList& inPL,
                 const CpuTimeLimit& CheckForTimeout)
  {
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis(GB, MinGens, inPL, CheckForTimeout);
    if (GB.empty())
    {
      swap(theLTResult,GB);
      return;
    }
    PolyList L;
    SparsePolyRing P(owner(*(GB.begin())));
    for (const RingElem& g: GB)
    {
      L.push_back(monomial(P,LPP(g)));
//??      SparsePolyIter it_p=BeginIter(g);
//??      L.push_back(monomial(P, PP(it_p)));
    }
    // for (PolyList::const_iterator it=GB.begin();it!=GB.end();++it)
    // {
    //   SparsePolyIter it_p=BeginIter(*it);
    //   L.push_back(monomial(P, PP(it_p)));
    // }
    swap(theLTResult,L);
  }//ComputeLT

}// end namespace cocoa


// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/TmpGOperations.C,v 1.70 2024/09/04 15:03:34 bigatti Exp $
// $Log: TmpGOperations.C,v $
// Revision 1.70  2024/09/04 15:03:34  bigatti
// Summary: Error codes: removed MixedSizes (into IncompatDims),
// removed BLASFailed, added ExternalLibs and ReqGradingDim1,
// changed ExpectedCoeffsInField into ReqCoeffsInField
// Now using CoCoA_THROW_ERROR1/2
//
// Revision 1.69  2024/05/28 12:57:46  bigatti
// Summary: added ComputeGBasisTrunc;
// changed myGBasis(G)/myMinGens(G)  into myCopyGBasis(G)/myCopyMinGens(G)
//
// Revision 1.68  2024/03/22 13:50:07  bigatti
// Summary: ComputeSaturationByPrincipal now cycles also on factors of monomial
//
// Revision 1.67  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.66  2024/02/15 15:32:57  bigatti
// Summary: added verbosity; added management of 1's in output of
// ComputeElim and ComputeSaturationHomogByPP
//
// Revision 1.65  2023/08/01 10:10:43  bigatti
// Summary: added ComputeSaturationHomogByPP
//
// Revision 1.64  2023/07/04 09:28:02  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.63  2023/07/03 13:51:44  bigatti
// Summary: intercepting gbasis of submodule with only 0 gens
//
// Revision 1.62  2023/03/27 14:00:25  abbott
// Summary: Commented out unused CheckForTimeout arg (in NYU fns)
//
// Revision 1.61  2023/03/14 22:35:22  abbott
// Summary: New iterator for loop syntax (redmine 1436)
//
// Revision 1.60  2022/02/18 14:12:00  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.59  2021/10/04 09:00:49  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.58  2021/10/04 07:55:04  bigatti
// Summary: using new homomorphism call, instead of "apply"  redmine 1598
//
// Revision 1.57  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.56  2018/07/04 13:12:09  bigatti
// -- better err mesg ERR::ExpectedCoeffsInField (ReqCoeffsInField)
//
// Revision 1.55  2018/06/27 12:15:19  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.54  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.53  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.52  2018/05/18 16:38:51  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.51  2018/05/18 12:24:47  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.50  2018/05/17 15:52:18  bigatti
// -- added include SparsePolyIter
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.49  2018/04/04 12:36:01  bigatti
// -- just indentation
//
// Revision 1.48  2018/02/22 16:56:41  bigatti
// -- error message
//
// Revision 1.47  2017/12/01 17:19:08  bigatti
// -- just spaces
//
// Revision 1.46  2017/11/29 17:41:48  bigatti
// -- added GBasisRealSolveCore
//
// Revision 1.45  2017/11/24 17:46:39  bigatti
// -- renamed GBasisSelfSat --> GBasisSelfSatCore
// -- added GBasisSelfSat in cpkg5
//
// Revision 1.44  2017/11/23 12:34:56  bigatti
// -- added GBasisSelfSat
//
// Revision 1.43  2017/05/22 21:11:44  abbott
// Summary: Replaced IsOne by IsInvertible (in saturate)
//
// Revision 1.42  2017/05/18 06:05:34  bigatti
// -- fixed ComputeSaturationByPrincipal with f = 1
//
// Revision 1.41  2017/05/17 15:59:48  bigatti
// -- calling NewPolyRing instead of NewPolyRing_DMPI in various places
//
// Revision 1.40  2017/04/18 09:26:37  bigatti
// -- removed StatLevel argument (now using VerbosityLevel)
//
// Revision 1.39  2016/11/23 14:00:45  bigatti
// -- now GBasis computes REducedBGasis
//
// Revision 1.38  2016/11/10 20:34:34  bigatti
// -- ComputeSaturationByPrincipal using factorization only for K == QQ
//
// Revision 1.37  2016/11/10 16:17:28  bigatti
// -- new ComputeSaturationByPrincipal (using factor)
//
// Revision 1.36  2015/12/04 15:22:24  bigatti
// -- renamed ComputeSSaturation into ComputeSaturation
//
// Revision 1.35  2015/10/01 10:15:26  bigatti
// -- just a newline
//
// Revision 1.34  2015/06/11 16:57:25  bigatti
// -- using new functions monomial(ring, pp) and monomial(ring, expv)
//
// Revision 1.33  2015/05/20 12:58:37  bigatti
// -- renamed CColon --> Colon
// -- some local name changes
//
// Revision 1.32  2015/03/04 10:29:26  bigatti
// -- added myDoGBasisElimFirst
// -- changed argument names: thePL --> inPL
//
// Revision 1.31  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.30  2014/07/30 14:10:18  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.29  2014/07/14 15:07:57  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.28  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.27  2014/07/08 08:39:14  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.26  2014/07/07 13:02:17  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.25  2014/04/30 16:14:38  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.24  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.23  2014/03/21 13:08:09  bigatti
// -- cosmetics
//
// Revision 1.22  2013/10/30 09:55:56  bigatti
// -- dealing with 0 generator in ComputeColonByPrincipal
//
// Revision 1.21  2013/06/12 08:54:07  bigatti
// -- added computation of MinGens (in ComputeGBasis)
// -- changed some "the" in "out"/"in" in argument names
//
// Revision 1.20  2013/02/21 17:16:45  bigatti
// -- changed syntax for ComputeSyz
//
// Revision 1.19  2013/02/12 16:36:49  bigatti
// -- added temporary (ehm..) function IsHomogGrD0 for GradingDim==0
//
// Revision 1.18  2013/01/31 11:43:22  bigatti
// -- added Stats field to ComputeXXXGBasis for returning statistics
//
// Revision 1.17  2012/10/03 12:24:46  bigatti
// -- some (..)ByPrincipal now take a RingElem instead of PolyList
// -- some cleaning
//
// Revision 1.16  2012/10/02 16:44:41  bigatti
// -- just aesthetics
//
// Revision 1.15  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.14  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.13  2011/12/07 15:54:34  bigatti
// -- renamed ambiguous "operator<" and hidden into anonymous namespace
// -- renamed ambiguous "operator==" into AreEqualLPPs (used by Saturation)
//
// Revision 1.12  2011/12/05 16:32:10  bigatti
// -- fixed bug about saturation (by non-principal ideal)
//
// Revision 1.11  2010/07/27 07:37:13  bigatti
// -- new class GBCriteria, simplified GReductor ctor
//
// Revision 1.10  2009/04/27 13:18:21  bigatti
// -- added SaturatingAlg flag
//
// Revision 1.9  2009/01/30 13:41:50  bigatti
// -- enum instead of bool arguments
//
// Revision 1.8  2008/09/19 13:33:42  bigatti
// -- added: Sat algorithm (M.Caboara)
//
// Revision 1.7  2008/09/19 11:34:15  bigatti
// -- new mechanism for passing verbosity level (or StatLevel)
//    [only partially tested]
//
// Revision 1.6  2008/07/09 16:09:47  abbott
// Cosmetic tidying.  Removed bogus comments. Added missing & (for C++ const reference).
//
// Revision 1.5  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.4  2007/11/09 10:45:52  bigatti
// -- [caboara] preparation for self-saturating algorithm
//
// Revision 1.3  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/03/27 13:19:21  bigatti
// -- removed printout
//
// Revision 1.1  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.17  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.16  2007/03/07 17:04:31  cocoa
// -- several changes by M.Caboara: more operations on ideals,
//    exception cleaner, coding conventions, WSugar, dynamic
//
// Revision 1.15  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.14  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.13  2007/01/25 12:05:20  bigatti
// -- fixed: removed coefficients in LT output
//
// Revision 1.12  2006/12/21 13:48:33  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.11  2006/12/04 13:55:54  cocoa
// -- added: sugar for GradingDim > 0  (called wsugar)
//
// Revision 1.10  2006/11/24 17:01:42  cocoa
// -- reorganized includes of header files
//
// Revision 1.9  2006/11/22 14:43:32  cocoa
// -- minor cleaning (indicated by Intel compiler)
//
// Revision 1.8  2006/11/20 14:57:17  cocoa
// -- added: (standard) sugar for modules
// -- fixed: non-homogeneous sysygies
// -- minor fixes     [M.Caboara]
//
// Revision 1.7  2006/10/06 16:46:17  cocoa
// -- syzygies for non-homogenous polynomials (Max)
// -- wip: evolution of Groebner Framework (Max)
//
// Revision 1.6  2006/10/06 10:48:34  cocoa
// Removed #include references to GradedFreeModule.H
//
// Revision 1.5  2006/08/17 09:35:49  cocoa
// -- changed: checks for homogeneity of input
//
// Revision 1.4  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.3  2006/07/18 10:52:24  cocoa
// -- added check for non-homogeneous input in ComputeElim
//
// Revision 1.2  2006/06/12 13:36:35  cocoa
// -- fixed embarrassing bug in WithDenominator1Hom clearing denominators
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.17  2006/05/16 08:59:16  cocoa
// -- added function for interactive Groebner
//
// Revision 1.16  2006/05/12 17:03:16  cocoa
// -- swapped arguments in homogenized
//
// Revision 1.15  2006/05/11 16:00:22  cocoa
// -- fixed spelling of "homogenize"
//
// Revision 1.14  2006/04/27 14:01:11  cocoa
// -- tidied up include files (using GTypes.H)
//
// Revision 1.13  2006/04/26 10:02:12  cocoa
// -- added GReductor::ourDefaultStatLevel variable to allow CoCoAServer
//    to set statistics level
//
// Revision 1.12  2006/04/21 16:47:06  cocoa
// -- new syntax for ComputeGBasis by Max
//
// Revision 1.11  2006/04/11 16:43:16  cocoa
// -- changed call of LPPNoCopy --> LPP (does not make a copy)
//
// Revision 1.10  2006/04/11 16:22:40  cocoa
// -- added: Elim, LT
//
