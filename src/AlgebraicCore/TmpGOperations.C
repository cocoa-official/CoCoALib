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
  { // namespace // anonymous ----------------------------------------------
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

  } // namespace // anonymous ----------------------------------------------
      

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
    VectorList tmpGB;
    VectorList tmpMinGens;
    GBR.myCopyGBasis(tmpGB);
    if (GradingDim(OldP)>0 && IsHomog(inGens)) GBR.myCopyMinGens(tmpMinGens);
    outGB.clear(); // just to remember to clean this up
    outMinGens.clear(); // just to remember to clean this up
    //    outGB =      DeEmbedPolyList(tmpGB, GRI);
    //    outMinGens = DeEmbedPolyList(tmpMinGens, GRI);
    outGB =      tmpGB;
    outMinGens = tmpMinGens;
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
      PolyList tmpGB;
      PolyList tmpMinGens;
      GBR.myCopyGBasis(tmpGB);
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myCopyMinGens(tmpMinGens);
      outGB.clear();//just to remember to clean this up
      outMinGens.clear();//just to remember to clean this up
      outGB = WithDenominator1Hom(tmpGB, SPR);
      MakeMonic(outGB);  //  2016-11-22: make monic
      outMinGens = WithDenominator1Hom(tmpMinGens, SPR);
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
      PolyList tmpGB;
      PolyList tmpMinGens;
      GBR.myCopyGBasis(tmpGB);
      outGB = WithDenominator1Hom(tmpGB, P);
      MakeMonic(outGB);  //  2016-11-22: make monic
      if (IsEveryWDegLtEq(inGens, TruncDeg)) GBR.myCopyMinGens(tmpMinGens);
      outMinGens = WithDenominator1Hom(tmpMinGens, P);
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
  
  
  void ComputeGBasisSelfSatCore(PolyList& outGB, const PolyList& inGens, const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGB.clear();
      return;
    }
    //    bool IsSatAlg=true;
    SparsePolyRing SPR(owner(inGens));
    PolyList tmpGB;
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  WithoutDenominators by EmbedVectorList
    //                   :  WithDenominator1Hom by DeEmbedPolyList
    if (!IsField(CoeffRing(SPR)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(SPR)))
    {
      const ring R = BaseRing(CoeffRing(SPR));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(SPR)), ordering(PPM(SPR)));
      GRingInfo GRI(Rx, IsHomogGrD0(inGens), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, WithoutDenominators(inGens, Rx),
                    GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      GBR.myCopyGBasis(tmpGB);
      tmpGB = WithDenominator1Hom(tmpGB, SPR);
    }
    else
    {
      GRingInfo GRI(SPR, IsHomogGrD0(inGens), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, inGens, GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      GBR.myCopyGBasis(tmpGB);
    }
    swap(outGB, tmpGB);
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
      PolyList tmpGB;
      PolyList tmpMinGens;
      GBR.myCopyGBasis(tmpGB);
      if (GradingDim(SPR)>0 && IsHomog(inGens)) GBR.myCopyMinGens(tmpMinGens);
      outGB.clear();//just to remember to clean this up
      outGB = WithDenominator1Hom(tmpGB, SPR);
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


  void ComputeElim(PolyList& outGens, const PolyList& inGens, ConstRefPPMonoidElem inds,
                   const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeElim");
    VERBOSE(99) << "-- called --" << std::endl;
    if (inGens.empty())
    {
      outGens.clear();
      return;
    }
    PolyList tmpElimGens;
    const SparsePolyRing OldSPR=owner(inGens);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    bool IsHomogGrD0PL = false;
    if  (GradingDim(OldSPR)!=0 && IsHomogGrD0(inGens)) IsHomogGrD0PL = true;
    SparsePolyRing NewSPR=MakeElimRingFromOld(OldSPR,IndexList, IsHomogGrD0PL);
    RingHom OldToNew = PolyAlgebraHom(OldSPR, NewSPR, indets(NewSPR));
    RingHom NewToOld = PolyAlgebraHom(NewSPR, OldSPR, indets(OldSPR));
    PolyList NewPL;
    for(const RingElem& g: inGens)
      NewPL.push_back(OldToNew(g));
    PPMonoidElem ElimIndsProd(LPP(OldToNew(monomial(OldSPR,inds))));
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis(GB, MinGens, NewPL, CheckForTimeout);
    for (const RingElem& g: GB)
    {
      if ( IsConstant(g) ) // redmine #1647
      {
        tmpElimGens.clear();
        VERBOSE(99) << "g is constant" << std::endl; /////////////////
        tmpElimGens.push_back(OldSPR->myOne());
        break;
      }
      if (IsCoprime(LPP(g), ElimIndsProd))
        tmpElimGens.push_back(NewToOld(g));
    }
    swap(outGens,tmpElimGens);
    return;
  }//ComputeElim


  RingElem ComputeElimFirst(const PolyList& inGens, ConstRefPPMonoidElem inds,
                            const CpuTimeLimit& CheckForTimeout)
  {
    const SparsePolyRing OldSPR=owner(inGens);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    bool IsHomogGrD0PL = (GradingDim(OldSPR)==0 ? false : IsHomogGrD0(inGens));
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
      GRingInfo GRI(NewSPR,IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, OldToNew(inGens));
      return NewToOld(GBR.myDoGBasisElimFirst(ElimIndsProd));
    }
    return zero(OldSPR); // just to keep the compiler quiet
  }


  void ComputeElim(VectorList& /*outGens*/,
                   const VectorList& /*theVL*/,
                   ConstRefPPMonoidElem /*inds*/,
                   const CpuTimeLimit& /*CheckForTimeout*/)
  {
    // Remember to check if ordering is TOPOS ordering here
    // if not, use hom
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeElim


  void ComputeSyz(VectorList& theSyzResult, const FreeModule& SyzFM, const VectorList& inGens,
            const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsHomogGrD0(inGens))
      CoCoA_THROW_ERROR2(ERR::NYI, "Not Yet Tested for non-homog input");
    if (inGens.empty())
    {
      theSyzResult.clear();
      return;
    }
    bool IsSatAlg=false;
    VectorList SyzResult;
    ModOrdTypeForcing MOType;
    if (IsHomogGrD0(inGens))
      MOType=WDegPosTO;
    else
      MOType=PosWDegTO;
    const FreeModule FM=owner(inGens);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM,MOType));
    const SparsePolyRing OldP(RingOf(FM));
    // Note: the GRI should build itself SyzFM and NewP from the data and deduce FM and OldP.
    //       All the embedding/deembedding functions should be members of GRI.
    GRingInfo GRI(NewP,OldP,FM,SyzFM,IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=SyzEmbedVectorList(inGens,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    SyzResult=DeEmbedPolyList(OutputPolyList,GRI,NumCompts(FM));
    swap(theSyzResult,SyzResult);
    return;
  }//ComputeSyz


  void ComputeSyz(VectorList& theSyzResult, const FreeModule& SyzFM, const PolyList& inGens,
            const CpuTimeLimit& CheckForTimeout)
  {
    //if (!IsHomogGrD0(inGens))
    //  CoCoA_THROW_ERROR2("Not Yet Tested for non-homogeneous input", "ComputeSyz(const VectorList&)");
    if (inGens.empty())
    {
      theSyzResult.clear();
      return;
    }
    bool IsSatAlg=false;
    VectorList SyzResult;
    ModOrdTypeForcing MOType;
    if (IsHomogGrD0(inGens))
      MOType=WDegPosTO;
    else
      MOType=PosWDegTO;
    const SparsePolyRing OldP(owner(inGens));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbedding(OldP,MOType));
    GRingInfo GRI(NewP,OldP,SyzFM,IsHomogGrD0(inGens),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=SyzEmbedPolyList(inGens,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    SyzResult=DeEmbedPolyList(OutputPolyList,GRI,1);
    swap(theSyzResult,SyzResult);
    return;
  }//ComputeSyz


  void ComputeIntersection(VectorList& theIntersectionResult,
                           const VectorList& inGens1,
                           const VectorList& inGens2,
                           const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    VectorList IntersectionResult;
    const FreeModule FM=owner(inGens1);
    //const SparsePolyRing NewP(MakeNewPRingFromModule(FM,WDegPosTO));
    const SparsePolyRing NewP(MakeNewPRingFromModulePosFirst(FM,IsHomogGrD0(inGens1)&&IsHomogGrD0(inGens2)));
    const SparsePolyRing OldP(RingOf(FM));
    GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(inGens1)&&IsHomogGrD0(inGens2),
                  IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=IntEmbedVectorLists(inGens1,inGens2,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    IntersectionResult=DeEmbedPolyList(OutputPolyList,GRI,NumCompts(FM));
    swap(IntersectionResult,theIntersectionResult);
    return;
  }//ComputeIntersection


  void ComputeIntersection(PolyList& theIntersectionResult,
                           const PolyList& inGens1,
                           const PolyList& inGens2,
                           const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    PolyList IntersectionResult;
    const SparsePolyRing OldP(owner(inGens1));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbeddingPosFirst(OldP,IsHomogGrD0(inGens1)&&IsHomogGrD0(inGens2)));
    GRingInfo GRI(NewP,OldP,IsHomogGrD0(inGens1)&&IsHomogGrD0(inGens2),IsSatAlg,
                  NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=IntEmbedPolyLists(inGens1,inGens2,GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    IntersectionResult=DeEmbedPolyListToPL(OutputPolyList,GRI,1);
    swap(theIntersectionResult,IntersectionResult);
    return;
  }//ComputeIntersection


// Colon by a single vector
  void ComputeColonByPrincipal(PolyList& theColonResult, const VectorList& inGens1, const VectorList& inGens2,
            const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    //    for (VectorList::const_iterator it=inGens1.begin();it!=inGens1.end();++it)
    //      if (!IsHomogGrD0(*it))
      if (!IsHomogGrD0(inGens1))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
      //    for (VectorList::const_iterator it=inGens2.begin();it!=inGens2.end();++it)
      //      if (!IsHomogGrD0(*it))
      if (!IsHomogGrD0(inGens2))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
    PolyList ColonResult;
    const FreeModule FM=owner(inGens1);
    const SparsePolyRing NewP(MakeNewPRingFromModule(FM,WDegPosTO));
    const SparsePolyRing OldP(RingOf(FM));
    GRingInfo GRI(NewP,OldP,FM,FM,IsHomogGrD0(inGens1)&&IsHomogGrD0(inGens2),
                  IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=ColonEmbedVectorLists(inGens1,inGens2,GRI);
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
                    const VectorList& inGens1,
                    const VectorList& inGens2,
                    const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens1.empty() && inGens2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    PolyList ColonResult;
    if (inGens2.empty())
    {
      ColonResult.push_back(one(RingOf(owner(inGens1))));
      swap(ColonResult, theResult);
      return;
    }
    if (inGens1.empty())
    {
      swap(ColonResult,theResult);
      return;
    }
    PolyList tmp;
    VectorList::const_iterator it=inGens2.begin();
    ComputeColonByPrincipal(ColonResult,inGens1,MakeVectorList(*it), CheckForTimeout);
    ++it;
    for (;it!=inGens2.end();++it)
    {
      ComputeColonByPrincipal(tmp,inGens1,MakeVectorList(*it), CheckForTimeout);
      ComputeIntersection(ColonResult,ColonResult,tmp, CheckForTimeout);
    }
    swap(ColonResult,theResult);
    return;
  }//ComputeColon


  void ComputeColonByPrincipal(PolyList& theResult,
                               const PolyList& inGens1,
                               ConstRefRingElem f,
                               const CpuTimeLimit& CheckForTimeout)
  {
    PolyList tmpColonResult;
    if (IsZero(f))
    {
      tmpColonResult.push_back(one(owner(f)));
      swap(theResult, tmpColonResult);
      return;
    }
    if (inGens1.empty())
    {
      theResult.clear();
      return;
    }
    bool IsSatAlg=false;
    const SparsePolyRing OldP(owner(inGens1));
    const SparsePolyRing NewP(MakeNewPRingForSimpleEmbeddingPosFirst(OldP,IsHomogGrD0(inGens1)&&IsHomogGrD0(f)));
    GRingInfo GRI(NewP,OldP,IsHomogGrD0(inGens1)&&IsHomogGrD0(f),
                  IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=ColonEmbedPolyLists(inGens1, std::vector<RingElem>(1,f), GRI);
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    PolyList OutputPolyList;
    GBR.myCopyGBasis(OutputPolyList);
    tmpColonResult = DeEmbedPolyListToPL(OutputPolyList, GRI, 1);
    swap(theResult, tmpColonResult);
    return;
  }//ComputeColonByPrincipal


// Colon ideal:ideal
  void ComputeColon(PolyList& theResult,
                    const PolyList& inGens1,
                    const PolyList& inGens2,
                    const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens1.empty() && inGens2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    PolyList tmpColonResult;
    if (inGens2.empty())
    {
      tmpColonResult.push_back(one(owner(inGens1)));
      swap(theResult, tmpColonResult);
      return;
    }
    if (inGens1.empty()) // ComputeColonByPrincipal checks this again
    {
      theResult.clear();
      return;
    }
    PolyList tmp;
    PolyList::const_iterator it=inGens2.begin();
    ComputeColonByPrincipal(tmpColonResult, inGens1, *it, CheckForTimeout);
    for (++it; it!=inGens2.end(); ++it)
    {
      ComputeColonByPrincipal(tmp, inGens1, *it, CheckForTimeout);
      ComputeIntersection(tmpColonResult, tmpColonResult, tmp, CheckForTimeout);
      if (tmpColonResult.empty()) break;
    }
    swap(theResult, tmpColonResult);
    return;
  }//ComputeColon


  void ComputeColonByPrincipal(VectorList& /*theResult*/,
                               const VectorList& /*inGens*/,
                               const PolyList& /*inGens*/,
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
    
    
    bool AreEqualLPPsIter(std::vector<PolyList::const_iterator>& theV1,
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
    
    
    // Useful when you have I\subset J and you want to check I==J
    bool AreEqualLPPs(const PolyList& inGens1, const PolyList& inGens2)
    {
      std::vector<PolyList::const_iterator> V1,V2;
      for (PolyList::const_iterator it=inGens1.begin(); it!=inGens1.end(); ++it)
        V1.push_back(it);
      for (PolyList::const_iterator it=inGens2.begin(); it!=inGens2.end(); ++it)
        V2.push_back(it);
      return AreEqualLPPsIter(V1,V2);
    }    
    
  } // anonymous namespace
  
   
//---------------------------------------------------------


  void ComputeSaturationByIrred(/*const*/ PolyList& tmpPL,
                                ConstRefRingElem f,
                                const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturationByIrred");
    VERBOSE(99) << "-- called --" << std::endl;
    const RingElem h = (indets(owner(f))).back();
    tmpPL.push_back(h*f-1);
    ComputeElim(tmpPL, tmpPL, LPP(h), CheckForTimeout);
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

  } // namespace anonymous
  

  void ComputeSaturationHomogByPP(PolyList& outPL,
                                  const PolyList& inGens,
                                  ConstRefPPMonoidElem t,
                                  const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturationHomogByPP");
    VERBOSE(99) << "-- called --" << std::endl;
//    if (IsZero(I))  return ideal(zero(Ph));
    const ring& P = owner(inGens);
    const ring& K = CoeffRing(P); // CoCoA_ASSERT(IsField(k));
    const long GrDim = GradingDim(P);
    const long NumX = NumIndets(P);
    std::vector<symbol> names;  for (long i=0; i<NumX; ++i) { names.push_back(IndetSymbol(P,i)); }
    VERBOSE(90) << "GrDim = " << GrDim << "   W = " << GradingMat(P) << std::endl;
    matrix WRL0 = NewDenseMat(ConcatVer(GradingMat(P), ZeroMat(RingZZ(), 1, NumIndets(P))));
    PolyList tmpPL = inGens;
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
                                    const PolyList& inGens,
                                    ConstRefRingElem f,
                                    const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturationByPrincipal");
    VERBOSE(99) << "-- called --" << std::endl;
    // non-homogeneous
    if (IsInvertible(f))
    {
      PolyList tmpPL(inGens);
      swap(outPL, tmpPL);
      return;
    }
    if (IsMonomial(f) /*IsPP*/ && IsHomog(inGens))
    {
      ComputeSaturationHomogByPP(outPL, inGens, LPP(f), CheckForTimeout);
      return;
    }
    const SparsePolyRing P = owner(inGens);
    const SparsePolyRing NewP = NewPolyRing(CoeffRing(P),
                                            NewSymbols(NumIndets(P)+1));
    const std::vector<RingElem>& x_new = indets(NewP);
    RingHom PToNewP = PolyAlgebraHom(P, NewP, std::vector<RingElem>(x_new.begin(), x_new.end()-1));
    std::vector<RingElem> x = indets(P);
    x.push_back(one(P));
    RingHom NewPToP = PolyAlgebraHom(NewP, P, x);
    PolyList tmpPL = PToNewP(inGens);
    // 2023-05-02: JAA: would be better to detect & handle specially when f is monomial!!
    // If f is not monomial & we can factorize, do a succession of saturations...
    //    if (!IsMonomial(f) && IsQQ(CoeffRing(P)))
    if (IsQQ(CoeffRing(P)))
    {
      const std::vector<RingElem> F = factor(f).myFactors();//const factorization<RingElem> F=factor(f);
      ComputeSaturationByIrred(tmpPL, PToNewP(F[0]), CheckForTimeout);
      for (long i=1; i<len(F); ++i)
        ComputeSaturationByIrred(tmpPL, PToNewP(F[i]), CheckForTimeout);
    }
    else // 2023-05-02: JAA: would be better to use radical(f) below.
      ComputeSaturationByIrred(tmpPL, PToNewP(f), CheckForTimeout);
    tmpPL = NewPToP(tmpPL);
    swap(outPL, tmpPL);
  } // ComputeSaturationByPrincipal


  void ComputeSaturation(PolyList& theSaturationResult,
                         const PolyList& inGens1,
                         const PolyList& inGens2,
                         const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturation");
    VERBOSE(99) << "-- called --" << std::endl;
    if (inGens1.empty() && inGens2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    if (inGens2.empty())
    {
      theSaturationResult.clear();// this or swap? this look better
      theSaturationResult.push_back(one(owner(inGens1)));
      return;
    }
    if (inGens1.empty())
    {
      theSaturationResult.clear();
      return;
    }

    PolyList tmpPL2;
    if (len(inGens2)==1)
      ComputeSaturationByPrincipal(tmpPL2, inGens1, inGens2.front(), CheckForTimeout);
    else
    {
      PolyList tmpPL1;
      ComputeColon(tmpPL1, inGens1, inGens2, CheckForTimeout);
      ComputeColon(tmpPL2, tmpPL1, inGens2, CheckForTimeout);
      while (!AreEqualLPPs(tmpPL1, tmpPL2))
      {
        swap(tmpPL1,tmpPL2);
        ComputeColon(tmpPL2, tmpPL1, inGens2, CheckForTimeout);
      }
    }
    swap(theSaturationResult,tmpPL2);
  }//ComputeSaturation


  void ComputeSaturationByPrincipal(VectorList& /*theSaturation*/,
                                    const VectorList& /*inGens*/,
                                    const PolyList& /*inGens*/,
                                    const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeSaturationByPrincipal


  void ComputeSaturation(VectorList& /*theSaturation*/,
                         const VectorList& /*inGens1*/,
                         const PolyList& /*inGens2*/,
                         const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeSaturation


  void ComputeHomogenization(VectorList& /*theHomogResult*/,
                             const VectorList& /*inGens1*/,
                             const PolyList& /*theIndets*/,
                             const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeHomogenization


  void ComputeHomogenization(PolyList& outGens,
                             const PolyList& inGens,
                             const PolyList& HomogIndets,
                             const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outGens.clear();
      return;
    }
    PolyList tmpHGens;
    const SparsePolyRing P=owner(inGens);
    RingElem IndetProd(product(HomogIndets));
    for (const RingElem& g: inGens)
      tmpHGens.push_back(homog(g, HomogIndets));
    ComputeSaturationByPrincipal(tmpHGens,tmpHGens,IndetProd,CheckForTimeout);
    swap(outGens, tmpHGens);
    return;
  }//ComputeHomogenization


// WARN: it supposes ComputeSaturationByPrincipal returns a GB
  bool RadicalMembership(const PolyList& inGens, ConstRefRingElem f,
                         const CpuTimeLimit& CheckForTimeout)
  {
    PolyList tmpGens;
    ComputeSaturationByPrincipal(tmpGens, inGens, f, CheckForTimeout);
    if (len(tmpGens) != 1) return false;
    return IsInvertible(tmpGens.front());
  }//RadicalMembership

  
  void ComputeLT(VectorList& /*theLTResult*/,
                 const VectorList& /*inGens*/,
                 const CpuTimeLimit& /*CheckForTimeout*/)
  {
    // Waiting for LT of a ModuleElem
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return;
  }//ComputeLT


  void ComputeLT(PolyList& outLTs, const PolyList& inGens,
                 const CpuTimeLimit& CheckForTimeout)
  {
    if (inGens.empty())
    {
      outLTs.clear();
      return;
    }
    PolyList tmpGB;
    PolyList tmpMinGens; // useless
    ComputeGBasis(tmpGB, tmpMinGens, inGens, CheckForTimeout);
    PolyList tmpLTs;
    SparsePolyRing P(owner(*(tmpGB.begin())));
    for (const RingElem& g: tmpGB)
      tmpLTs.push_back(monomial(P,LPP(g)));
    swap(outLTs, tmpLTs);
  }//ComputeLT

}// end namespace cocoa
