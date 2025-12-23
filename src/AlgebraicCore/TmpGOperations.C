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
      MakeMonic(outGB);  //  2016-11-22: make monic
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
  

  namespace // anonymous
  { // namespace // anonymous ----------------------------------------------
    
    bool IsEveryWDegLtEq(const PolyList& F, long D)
    {
      for (auto& f:F)  if (wdeg(f)[0] > D) return false;
      return true;
    }
    
  } // namespace // anonymous ----------------------------------------------
  

  void ComputeGBasisTrunc(PolyList& outGB, PolyList& outMinGens, long& TruncDeg, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGB.clear();
      outMinGens.clear();
      return;
    }
    const SparsePolyRing P(owner(inGens));
    CoCoA_ASSERT_ALWAYS(TruncDeg >= 0); // user TruncDeg must be >=0
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (GradingDim(P)!=1)  CoCoA_THROW_ERROR1(ERR::ReqGradingDim1);
    if (!IsHomog(inGens))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    bool IsSatAlg=false;
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    {
      //---------------------------------------------------
      const SparsePolyRing Rx = NewPolyRing(BaseRing(CoeffRing(P)), symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(inGens), IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(inGens, Rx));
      GBR.mySetTruncDeg(TruncDeg); // input value
      //---------------------------------------------------
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList TmpGB;
      PolyList TmpMinGens;
      GBR.myCopyGBasis(TmpGB);
      outGB = WithDenominator1Hom(TmpGB, P);
      MakeMonic(outGB);  //  2016-11-22: make monic
      if (IsEveryWDegLtEq(inGens, TruncDeg)) GBR.myCopyMinGens(TmpMinGens);
      outMinGens = WithDenominator1Hom(TmpMinGens, P);
      TruncDeg = GBR.myTruncDeg(); // update TruncDeg (if changed/complete)
    }
    else
    {
      //---------------------------------------------------
      GRingInfo GRI(P, IsHomogGrD0(inGens), IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, inGens);
      GBR.mySetTruncDeg(TruncDeg); // input value
      //---------------------------------------------------
      GBR.myDoGBasis();// homog input standard alg interred
      GBR.myCopyGBasis(outGB);
      outMinGens.clear();//just to remember to clean this up
      if (IsEveryWDegLtEq(inGens, TruncDeg)) GBR.myCopyMinGens(outMinGens);
      TruncDeg = GBR.myTruncDeg(); // update TruncDeg (if changed/complete)
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
      MakeMonic(outGB);  //  2016-11-22: make monic
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
