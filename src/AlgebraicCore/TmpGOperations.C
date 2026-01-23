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

    bool IsHomogGrD0(const ModuleElem v)
    {
      if (GradingDim(RingOf(owner(v)))==0) return true;
      return IsHomog(v);
    }


    PolyList KxToRx(const PolyList& F, SparsePolyRing Rx)
    {
      if (F.empty()) return F;
      CoCoA_ASSERT(HasUniqueOwner(F));
      const SparsePolyRing Kx(owner(F));
      CoCoA_ASSERT(IsFractionFieldOfGCDDomain(CoeffRing(Kx)));
      CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
      CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
      PolyList F_Rx;
      RingElem f_Rx(Rx);
      std::vector<long> expv(NumIndets(Rx));
      for (const RingElem& f: F)
      {
        const RingElem d(CommonDenom(f));
        f_Rx = 0;
        for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
          PushBack(f_Rx, num(d*coeff(it)), exponents(expv,PP(it))); // same ord
        F_Rx.push_back(f_Rx);
      }
      return F_Rx;
    }


    PolyList RxToKx(const PolyList& F, SparsePolyRing Kx)
    {
      CoCoA_ASSERT(IsFractionField(CoeffRing(Kx)));
      if (F.empty()) return F;
      CoCoA_ASSERT(HasUniqueOwner(F));
      const SparsePolyRing Rx(owner(F));
      CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
      CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
      PolyList outF;
      outF.reserve(len(F));
      const RingHom EmbedIntoFrF = EmbeddingHom(CoeffRing(Kx));
      const RingHom EmbedCoeff   = CoeffEmbeddingHom(Kx);
      const RingHom phi = PolyRingHom(Rx, Kx, EmbedCoeff(EmbedIntoFrF), indets(Kx));
      //    transform(F.begin(), F.end(), back_inserter(outF), phi);
      return phi(F);
    }

  } // namespace // anonymous ----------------------------------------------
      

// GBasis ==================================================================

  void ComputeGBasis2(VectorList& GB_out, VectorList& MinGens_out,
                      const VectorList& G_in, const CpuTimeLimit& CheckForTimeout)
  {
    if (G_in.empty())
    {
      GB_out.clear();
      MinGens_out.clear();
      return;
    }
    const FreeModule FM=owner(G_in);
    const SparsePolyRing P_new(MakeNewPRingFromModule(FM));
    const SparsePolyRing P(RingOf(FM));
    if (!IsField(CoeffRing(P)))
      CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool IsSatAlg=false;
    const GRingInfo GRI(P_new,P,FM,FM,IsHomogGrD0(G_in),
                        IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys=EmbedVectorList(G_in,GRI);
    if (EmbeddedPolys.empty())
    {
      GB_out.clear();
      MinGens_out.clear();
      return;
    }
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    VectorList GB_tmp;
    VectorList MinGens_tmp;
    GBR.myCopyGBasis_module(GB_tmp);
    if (GradingDim(P)>0 && IsHomog(G_in)) GBR.myCopyMinGens_module(MinGens_tmp);
    GB_out.clear(); // just to remember to clean this up
    MinGens_out.clear(); // just to remember to clean this up
    //    GB_out =      DeEmbedPolyList(GB_tmp, GRI);
    //    MinGens_out = DeEmbedPolyList(MinGens_tmp, GRI);
    GB_out =      GB_tmp;
    MinGens_out = MinGens_tmp;
  }//ComputeGBasis


  void ComputeGBasis2(PolyList& GB_out, PolyList& MinGens_out,
                      const PolyList& G_in, const CpuTimeLimit& CheckForTimeout)
  {
    if (G_in.empty())
    {
      GB_out.clear();
      MinGens_out.clear();
      return;
    }
    SparsePolyRing P(owner(G_in));
    // MAX: Adding to ComputeGBasis2 for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  KxToRx by EmbedVectorList
    //                   :  RxToKx by DeEmbedPolyList
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool IsSatAlg=false;
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    {
      const ring R = BaseRing(CoeffRing(P));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G_in),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G_in, Rx));
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList GB_tmp = RxToKx(GBR.myExportGBasis(), P);
      PolyList MinGens_tmp;
      if (GradingDim(P)>0 && IsHomog(G_in)) MinGens_tmp = RxToKx(GBR.myExportMinGens(), P);
      MakeMonic(GB_tmp);  //  2016-11-22: make monic
      swap(GB_out, GB_tmp);
      swap(MinGens_out, MinGens_tmp);
    }
    else
    {
      GRingInfo GRI(P,IsHomogGrD0(G_in),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G_in);
      GBR.myDoGBasis();// homog input standard alg interred
      GB_out = GBR.myExportGBasis();
      if (GradingDim(P)>0 && IsHomog(G_in)) MinGens_out = GBR.myExportMinGens();
    }          
  }//ComputeGBasis2
  

  namespace
  { // namespace // anonymous ----------------------------------------------
    
    bool IsEveryWDegLtEq(const PolyList& F, long D)
    {
      for (auto& f:F)  if (wdeg(f)[0] > D) return false;
      return true;
    }
    
  } // namespace // anonymous ----------------------------------------------
  

  void ComputeGBasisTrunc2(PolyList& GB_out, PolyList& MinGens_out,
                           long& TruncDeg, const PolyList& G_in, const CpuTimeLimit& CheckForTimeout)
  {
    if (G_in.empty())
    {
      GB_out.clear();
      MinGens_out.clear();
      return;
    }
    const SparsePolyRing P(owner(G_in));
    CoCoA_ASSERT_ALWAYS(TruncDeg >= 0); // user TruncDeg must be >=0
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (GradingDim(P)!=1)  CoCoA_THROW_ERROR1(ERR::ReqGradingDim1);
    if (!IsHomog(G_in))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    bool IsSatAlg=false;
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    { //---------------------------------------------------
      const SparsePolyRing Rx = NewPolyRing(BaseRing(CoeffRing(P)), symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G_in), IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G_in, Rx));
      GBR.mySetTruncDeg(TruncDeg); // input value
      //---------------------------------------------------
      GBR.myDoGBasis();// homog input standard alg interred
      PolyList GB_tmp = RxToKx(GBR.myExportGBasis(), P);
      PolyList MinGens_tmp;
      MakeMonic(GB_tmp);  //  2016-11-22: make monic
      if (IsEveryWDegLtEq(G_in, TruncDeg)) MinGens_tmp = RxToKx(GBR.myExportMinGens(), P);
      swap(GB_out, GB_tmp);
      swap(MinGens_out, MinGens_tmp);
      TruncDeg = GBR.myTruncDeg(); // update TruncDeg (if changed/complete)
    }
    else
    { //---------------------------------------------------
      GRingInfo GRI(P, IsHomogGrD0(G_in), IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G_in);
      GBR.mySetTruncDeg(TruncDeg); // input value
      //---------------------------------------------------

      GBR.myDoGBasis();// homog input standard alg interred
      GB_out = GBR.myExportGBasis();
      if (IsEveryWDegLtEq(G_in, TruncDeg)) MinGens_out = GBR.myExportMinGens();
      TruncDeg = GBR.myTruncDeg(); // update TruncDeg (if changed/complete)
    }
  }//ComputeGBasisTrunc
  
  
  PolyList ComputeGBasisSelfSatCore(const PolyList& G, const CpuTimeLimit& CheckForTimeout)
  {
    if (G.empty()) return G;
    //    bool IsSatAlg=true;
    SparsePolyRing P(owner(G));
    // MAX: Adding to ComputeGBasis for modules should suffice.
    // The new ring should be created by MakeNewPRingFromModule
    // and the embeddings:  KxToRx by EmbedVectorList
    //                   :  RxToKx by DeEmbedPolyList
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    {
      const ring R = BaseRing(CoeffRing(P));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G, Rx), GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      return RxToKx(GBR.myExportGBasis(), P);
    }
    else
    {
      GRingInfo GRI(P, IsHomogGrD0(G), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G, GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      return GBR.myExportGBasis();
    }
  }


  PolyList ComputeGBasisRealSolve(const PolyList& G, const CpuTimeLimit& CheckForTimeout)
  {
    if (G.empty()) return G;
    SparsePolyRing P(owner(G));
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (characteristic(P)!=0)  CoCoA_THROW_ERROR1(ERR::ReqChar0);
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    {
      const ring R = BaseRing(CoeffRing(P));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G), false/*IsSatAlg*/,NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G, Rx));
      GBR.myDoGBasisRealSolve();// homog input standard alg interred
      PolyList GB_tmp = RxToKx(GBR.myExportGBasis(), P);
      MakeMonic(GB_tmp);  //  2016-11-22: make monic
      return GB_tmp;
    }
    else
    {
      GRingInfo GRI(P,IsHomogGrD0(G),false/*IsSatAlg*/,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G);
      GBR.myDoGBasisRealSolve();
      return GBR.myExportGBasis();
    }          
  }


  PolyList ComputeElim(const PolyList& G, ConstRefPPMonoidElem inds,
                       const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeElim");
    VERBOSE(99) << "-- called --" << std::endl;
    if (G.empty())  return G;
    const SparsePolyRing P = owner(G);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    bool HomogInput = false;
    if  (GradingDim(P)!=0 && IsHomog(G)) HomogInput = true;
    SparsePolyRing P_new=MakeElimRing(P, IndexList, HomogInput);
    RingHom PToNew = PolyAlgebraHom(P, P_new, indets(P_new));
    RingHom NewToP = PolyAlgebraHom(P_new, P, indets(P));
    PolyList NewGens;
    for(const RingElem& g: G) NewGens.push_back(PToNew(g));
    PPMonoidElem ElimIndsProd(LPP(PToNew(monomial(P,inds))));
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis2(GB, MinGens, NewGens, CheckForTimeout);
    PolyList ElimGens;
    for (const RingElem& g: GB)
    {
      if ( IsConstant(g) ) // redmine #1647
      {
        ElimGens.clear();
        VERBOSE(99) << "g is constant" << std::endl; /////////////////
        ElimGens.push_back(one(P));
        break;
      }
      if (IsCoprime(LPP(g), ElimIndsProd))
        ElimGens.push_back(NewToP(g));
    }
    return ElimGens;
  }//ComputeElim


  RingElem ComputeElimFirst(const PolyList& G, ConstRefPPMonoidElem inds,
                            const CpuTimeLimit& CheckForTimeout)
  {
    const SparsePolyRing P=owner(G);
    std::vector<long> IndexList=PPMonoidElem2IndexList(inds);
    bool IsHomogGrD0PL = (GradingDim(P)==0 ? false : IsHomogGrD0(G));
    if (!IsHomogGrD0PL)  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    bool IsSatAlg = false;
    SparsePolyRing P_new = MakeElimRing(P,IndexList, IsHomogGrD0PL);
    RingHom OldToNew = PolyAlgebraHom(P, P_new, indets(P_new));
    RingHom NewToOld = PolyAlgebraHom(P_new, P, indets(P));
    PPMonoidElem ElimIndsProd = LPP(OldToNew(monomial(P,inds)));
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
      CoCoA_THROW_ERROR2(ERR::NYI, "Only for FFp");
    else
    {
      GRingInfo GRI(P_new,IsHomogGrD0(G),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, OldToNew(G));
      return NewToOld(GBR.myDoGBasisElimFirst(ElimIndsProd));
    }
    return zero(P); // just to keep the compiler quiet
  }


  VectorList ComputeElim(const VectorList& G, ConstRefPPMonoidElem /*inds*/,
                         const CpuTimeLimit& /*CheckForTimeout*/)
  {
    // Remember to check if ordering is TOPOS ordering here
    // if not, use hom
    CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
    return G; // just to keep the compiler quiet
  }//ComputeElim


  VectorList ComputeSyz(const FreeModule& SyzFM, const VectorList& G,
                        const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsHomogGrD0(G))
      CoCoA_THROW_ERROR2(ERR::NYI, "Not Yet Tested for non-homog input");
    if (G.empty())  return G;
    bool IsSatAlg=false;
    ModOrdTypeForcing MOType = (IsHomogGrD0(G) ? WDegPosTO : PosWDegTO);
    const FreeModule FM=owner(G);
    const SparsePolyRing P(RingOf(FM));
    const SparsePolyRing P_new(MakeNewPRingFromModule(FM,MOType));
    // Note: the GRI should build itself SyzFM and P_new from the data and deduce FM and P.
    //       All the embedding/deembedding functions should be members of GRI.
    GRingInfo GRI(P_new,P,FM,SyzFM,IsHomogGrD0(G),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, SyzEmbedVectorList(G,GRI));
    GBR.myDoGBasis();
    return DeEmbedPolyList(GBR.myExportGBasis(), GRI,NumCompts(FM));
  }//ComputeSyz


  VectorList ComputeSyz(const FreeModule& SyzFM, const PolyList& G,
                        const CpuTimeLimit& CheckForTimeout)
  {
    //if (!IsHomogGrD0(G))
    //  CoCoA_THROW_ERROR2("Not Yet Tested for non-homogeneous input", "ComputeSyz(const VectorList&)");
    if (G.empty()) return VectorList{zero(SyzFM)};
    bool IsSatAlg=false;
    ModOrdTypeForcing MOType = (IsHomogGrD0(G) ? WDegPosTO : PosWDegTO);
    const SparsePolyRing P(owner(G));
    const SparsePolyRing P_new(MakeNewPRingForSimpleEmbedding(P,MOType));
    GRingInfo GRI(P_new,P,SyzFM,IsHomogGrD0(G),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, SyzEmbedPolyList(G,GRI));
    GBR.myDoGBasis();
    return DeEmbedPolyList(GBR.myExportGBasis(), GRI,1);
  }//ComputeSyz


  VectorList ComputeIntersection(const VectorList& G1, const VectorList& G2,
                                 const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    const FreeModule FM=owner(G1);
    //const SparsePolyRing P_new(MakeNewPRingFromModule(FM,WDegPosTO));
    bool HomogInput = IsHomogGrD0(G1) && IsHomogGrD0(G2);
    const SparsePolyRing P_new(MakeNewPRingFromModulePosFirst(FM,HomogInput));
    const SparsePolyRing P(RingOf(FM));
    GRingInfo GRI(P_new, P, FM, FM, HomogInput,
                  IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, IntEmbedVectorLists(G1, G2, GRI));
    GBR.myDoGBasis();// homog input standard alg interred
    return DeEmbedPolyList(GBR.myExportGBasis(), GRI,NumCompts(FM));
  }//ComputeIntersection


  PolyList ComputeIntersection(const PolyList& G1, const PolyList& G2,
                               const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    const SparsePolyRing P(owner(G1));
    const bool HomogInput = IsHomogGrD0(G1) && IsHomogGrD0(G2);
    const SparsePolyRing P_new(MakeNewPRingForSimpleEmbeddingPosFirst(P, HomogInput));
    GRingInfo GRI(P_new, P, HomogInput,
                  IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, IntEmbedPolyLists(G1, G2, GRI));
    GBR.myDoGBasis();// homog input standard alg interred
    return DeEmbedPolyListToPL(GBR.myExportGBasis(), GRI, 1);
  }//ComputeIntersection


  namespace
  { // namespace // anonymous ----------------------------------------------

    PolyList ComputeColonByPrincipal(const VectorList& G1, const ModuleElem& v,
                                     const CpuTimeLimit& CheckForTimeout)
    {
      CoCoA_ASSERT_ALWAYS(!G1.empty());  // Anna togli always
      bool IsSatAlg=false;
      if (!IsHomogGrD0(G1))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
      if (!IsHomogGrD0(v))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
      const FreeModule FM = owner(G1);
      const SparsePolyRing P_new(MakeNewPRingFromModule(FM,WDegPosTO));
      const SparsePolyRing P(RingOf(FM));
      GRingInfo GRI(P_new,P,FM,FM,IsHomogGrD0(G1)&&IsHomogGrD0(v),
                    IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, ColonEmbedVectorLists(G1, VectorList(1,v), GRI));
      GBR.myDoGBasis();// homog input standard alg interred
      return DeEmbedPolyListToPL(GBR.myExportGBasis(), GRI, NumCompts(FM));
    }
    

    PolyList ComputeColonByPrincipal(const PolyList& G, ConstRefRingElem f,
                                     const CpuTimeLimit& CheckForTimeout)
    {
      CoCoA_ASSERT_ALWAYS(!G.empty());  //-- Anna  togli always
      PolyList tmpColonResult;
      if (IsZero(f))  return std::vector<RingElem>{one(owner(f))};
      bool IsSatAlg=false;
      const SparsePolyRing P(owner(G));
      const SparsePolyRing P_new(MakeNewPRingForSimpleEmbeddingPosFirst(P,IsHomogGrD0(G)&&IsHomogGrD0(f)));
      GRingInfo GRI(P_new, P, IsHomogGrD0(G) && IsHomogGrD0(f),
                    IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GPolyList EmbeddedPolys=ColonEmbedPolyLists(G, std::vector<RingElem>{f}, GRI);
      GReductor GBR(GRI, EmbeddedPolys);
      GBR.myDoGBasis();// homog input standard alg interred
      return DeEmbedPolyListToPL(GBR.myExportGBasis(), GRI, 1);
    }

    VectorList ComputeColonByPrincipal(const VectorList& G, const PolyList& /*G_in*/,
                                       const CpuTimeLimit& /*CheckForTimeout*/)
    {
      CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
      return G; // just to keep the compiler quiet
    }
  

  } // namespace // anonymous ----------------------------------------------


  PolyList ComputeColon(const VectorList& G1, const VectorList& G2,
                        const CpuTimeLimit& CheckForTimeout)
  {
    if (G1.empty() && G2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    if (G2.empty())  return std::vector<RingElem>{one(RingOf(owner(G1)))};
    if (G1.empty())  return std::vector<RingElem>{};
    PolyList aux;
    VectorList::const_iterator it=G2.begin();
    PolyList tmpColon = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
    for (++it; it!=G2.end(); ++it)
    {
      aux = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
      tmpColon = ComputeIntersection(tmpColon, aux, CheckForTimeout);
      if (tmpColon.empty()) break;  /// anna non puo' essere empty, giusto?
    }
    return tmpColon;
  }


  PolyList ComputeColon(const PolyList& G1, const PolyList& G2,
                        const CpuTimeLimit& CheckForTimeout)
  {
    if (G1.empty() && G2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    if (G2.empty()) return std::vector<RingElem>{one(owner(G1))};
    if (G1.empty()) return G1;
    PolyList aux;
    PolyList::const_iterator it=G2.begin();
    PolyList tmpColon = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
    for (++it; it!=G2.end(); ++it)
      if (!IsZero(*it))
      {
        aux = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
        tmpColon = ComputeIntersection(tmpColon, aux, CheckForTimeout);
        if (tmpColon.empty()) break;
      }
    return tmpColon;
  }


  namespace
  { // namespace // anonymous ----------------------------------------------
    // These procedures are for EqualLPPs used by Saturation
  
    // PolyList::const_iterators are ordered according to LPP of their polys
    bool ByLPP(const PolyList::const_iterator& it,
               const PolyList::const_iterator& it2)
    {
      return (SparsePolyRingPtr(owner(*it))->myCmpLPP(raw(*it),raw(*it2)) == -1);
    }
    
    
    bool AreEqualLPPsIter(std::vector<PolyList::const_iterator>& v1,
                          std::vector<PolyList::const_iterator>& v2)
    {  // Auxiliary: v1 and v2 have same length
      stable_sort(v1.begin(), v1.end(), ByLPP);
      stable_sort(v2.begin(), v2.end(), ByLPP);
      const long len1 = len(v1);
      for (int i=0; i!=len1; ++i)
        if (LPP(*(v1[i])) != LPP(*(v2[i])))  return false;
      return true;
    }
    
    
    // Useful when you have I\subset J and you want to check I==J
    bool AreEqualLPPs(const PolyList& G1, const PolyList& G2)
    {
      if (len(G1) != len(G2))  return false;
      std::vector<PolyList::const_iterator> v1,v2;
      for (PolyList::const_iterator it=G1.begin(); it!=G1.end(); ++it)
        v1.push_back(it);
      for (PolyList::const_iterator it=G2.begin(); it!=G2.end(); ++it)
        v2.push_back(it);
      return AreEqualLPPsIter(v1,v2);
    }

    
    PolyList ComputeSaturationByIrred(const PolyList& G, ConstRefRingElem f,
                                      const CpuTimeLimit& CheckForTimeout)
    {
      VerboseLog VERBOSE("ComputeSaturationByIrred");
      VERBOSE(99) << "-- called --" << std::endl;
      const RingElem h = (indets(owner(f))).back();
      PolyList tmpPL = G;
      tmpPL.push_back(h*f-1);
      return ComputeElim(tmpPL, LPP(h), CheckForTimeout);
    } // ComputeSaturationByIrred

    
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

    
    PolyList ComputeSaturationHomogByPP(const PolyList& G,
                                        ConstRefPPMonoidElem t,
                                        const CpuTimeLimit& CheckForTimeout)
    {
      VerboseLog VERBOSE("ComputeSaturationHomogByPP");
      VERBOSE(99) << "-- called --" << std::endl;
      //    if (IsZero(I))  return ideal(zero(Ph));
      const ring& P = owner(G);
      const ring& K = CoeffRing(P); // CoCoA_ASSERT(IsField(k));
      const long GrDim = GradingDim(P);
      const long NumX = NumIndets(P);
      std::vector<symbol> names;  for (long i=0; i<NumX; ++i) { names.push_back(IndetSymbol(P,i)); }
      VERBOSE(90) << "GrDim = " << GrDim << "   W = " << GradingMat(P) << std::endl;
      matrix WRL0 = NewDenseMat(ConcatVer(GradingMat(P), ZeroMat(RingZZ(), 1, NumIndets(P))));
      PolyList tmpGens = G;
      PolyList MinGens; // unused
      std::vector<long> expv = exponents(t);
      for (long i=0; i<NumIndets(P); ++i)
        if (expv[i] != 0)
        {
          VERBOSE(95) << "doing indet " << i << " = " << indet(P,i) << std::endl;
          SetEntry(WRL0, GrDim,i, -1);
          ring PzDegRevLex = NewPolyRing(K, names, NewMatrixOrdering(MakeTermOrdMat(WRL0), GrDim));
          SetEntry(WRL0, GrDim,i, 0);
          RingHom phi = PolyAlgebraHom(owner(tmpGens[0]), PzDegRevLex, indets(PzDegRevLex));
          ComputeGBasis2(tmpGens, MinGens, phi(tmpGens), CheckForTimeout);
          for (RingElem& g: tmpGens)  g = SatByIndet(g,i);
        }
      for (auto& g:tmpGens)
      {
        if ( IsZero(g) )  CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
        if ( IsConstant(g) ) // redmine #1647
        {
          VERBOSE(99) << "g is constant" << std::endl; /////////////////
          return std::vector<RingElem>(1, one(P));
        }
      }
      RingHom psi = PolyAlgebraHom(owner(tmpGens[0]), P, indets(P));
      return psi(tmpGens);
    } // ComputeSaturationHomogByPP

  
    PolyList ComputeSaturationByPrincipal(const PolyList& G, ConstRefRingElem f,
                                          const CpuTimeLimit& CheckForTimeout)
    {
      VerboseLog VERBOSE("ComputeSaturationByPrincipal");
      VERBOSE(99) << "-- called --" << std::endl;
      // non-homogeneous
      if (IsInvertible(f))  return G;
      if (IsMonomial(f) /*IsPP*/ && IsHomog(G))
        return ComputeSaturationHomogByPP(G, LPP(f), CheckForTimeout);
      const SparsePolyRing P = owner(G);
      const SparsePolyRing P_new = NewPolyRing(CoeffRing(P),
                                               NewSymbols(NumIndets(P)+1));
      const std::vector<RingElem>& x_new = indets(P_new);
      RingHom PToP_new = PolyAlgebraHom(P, P_new, std::vector<RingElem>(x_new.begin(), x_new.end()-1));
      std::vector<RingElem> x = indets(P);
      x.push_back(one(P));
      RingHom P_newToP = PolyAlgebraHom(P_new, P, x);
      PolyList tmpGens = PToP_new(G);
      // 2023-05-02: JAA: would be better to detect & handle specially when f is monomial!!
      // If f is not monomial & we can factorize, do a succession of saturations...
      //    if (!IsMonomial(f) && IsQQ(CoeffRing(P)))
      if (IsQQ(CoeffRing(P)))
      {
        const std::vector<RingElem> F = factor(f).myFactors();//const factorization<RingElem> F=factor(f);
        tmpGens = ComputeSaturationByIrred(tmpGens, PToP_new(F[0]), CheckForTimeout);
        for (long i=1; i<len(F); ++i)
          tmpGens = ComputeSaturationByIrred(tmpGens, PToP_new(F[i]), CheckForTimeout);
      }
      else // 2023-05-02: JAA: would be better to use radical(f) below.
        tmpGens = ComputeSaturationByIrred(tmpGens, PToP_new(f), CheckForTimeout);
      return P_newToP(tmpGens);
    } // ComputeSaturationByPrincipal


    VectorList ComputeSaturationByPrincipal(const VectorList& GensM, const PolyList& /*G_in*/,
                                            const CpuTimeLimit& /*CheckForTimeout*/)
    {
      CoCoA_THROW_ERROR2(ERR::NYI, "VectorList");
      return GensM; // just to keep the compiler quiet
    }


  } // namespace // anonymous ----------------------------------------------


  PolyList ComputeSaturation(const PolyList& G1, const PolyList& G2,
                             const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeSaturation");
    VERBOSE(99) << "-- called --" << std::endl;
    if (G1.empty() && G2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both G1 and G2 are empty");
    if (G2.empty())  return std::vector<RingElem>{one(owner(G1))};
    if (G1.empty())  return G1;
    if (len(G2)==1)
      return ComputeSaturationByPrincipal(G1, G2.front(), CheckForTimeout);
    PolyList aux = ComputeColon(G1, G2, CheckForTimeout);
    PolyList tmpGens = ComputeColon(aux, G2, CheckForTimeout);
    while (!AreEqualLPPs(aux, tmpGens))
    {
      swap(aux, tmpGens);
      tmpGens = ComputeColon(aux, G2, CheckForTimeout);
    }
    return tmpGens;
  }//ComputeSaturation


  VectorList ComputeSaturation(const VectorList& G1, const PolyList& /*G2*/,
                               const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "for modules");
    return G1;  // just to keep the compiler quiet
  }


  VectorList ComputeHomogenization(const VectorList& G, const PolyList& /*theIndets*/,
                             const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "for modules");
    return G;  // just to keep the compiler quiet
  }


  PolyList ComputeHomogenization(const PolyList& G, const PolyList& HomogIndets,
                                 const CpuTimeLimit& CheckForTimeout)
  {
    if (G.empty()) return G;
    RingElem IndetProd(product(HomogIndets));
    PolyList tmpHGens;
    for (const RingElem& g: G) tmpHGens.push_back(homog(g, HomogIndets));
    return ComputeSaturationByPrincipal(tmpHGens,IndetProd,CheckForTimeout);
  }


// WARN: it supposes ComputeSaturationByPrincipal returns a GB
  bool RadicalMembership(const PolyList& G, ConstRefRingElem f,
                         const CpuTimeLimit& CheckForTimeout)
  {
    PolyList tmpGens;
    tmpGens = ComputeSaturationByPrincipal(G, f, CheckForTimeout);
    if (len(tmpGens) != 1) return false;
    return IsInvertible(tmpGens.front());
  }

  
}// end namespace cocoa
