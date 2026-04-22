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
#include "CoCoA/ModuleOrdering.H" // for WDegPosnOrd
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

namespace CoCoA
{


  namespace // anonymous
  { // namespace // anonymous ----------------------------------------------
    bool IsHomogGrD0(const std::vector<RingElem>& L)
    {
      if (L.empty()) return true;
      if (GradingDim(owner(L[0]))==0) return true;
      return IsHomog(L);
    }

    bool IsHomogGrD0(const std::vector<ModuleElem>& L)
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


    std::vector<RingElem> KxToRx(const std::vector<RingElem>& F, SparsePolyRing Rx)
    {
      if (F.empty()) return F;
      CoCoA_ASSERT(HasUniqueOwner(F));
      const SparsePolyRing Kx(owner(F[0]));
      CoCoA_ASSERT(IsFractionFieldOfGCDDomain(CoeffRing(Kx)));
      CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
      CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
      std::vector<RingElem> F_Rx;  F_Rx.reserve(len(F));
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


    std::vector<RingElem> RxToKx(const std::vector<RingElem>& F, SparsePolyRing Kx)
    {
      CoCoA_ASSERT(IsFractionField(CoeffRing(Kx)));
      if (F.empty()) return F;
      CoCoA_ASSERT(HasUniqueOwner(F));
      const SparsePolyRing Rx(owner(F[0]));
      CoCoA_ASSERT(BaseRing(CoeffRing(Kx)) == CoeffRing(Rx));
      CoCoA_ASSERT(ordering(PPM(Kx)) == ordering(PPM(Rx)));
      const RingHom RToK = EmbeddingHom(CoeffRing(Kx));
      const RingHom phi = PolyRingHom(Rx, Kx, RToK, indets(Kx));
      return phi(F);
    }

    ////////////////////////////////////
    // from TmpGReductor 2026-04-10
    ////////////////////////////////////
    
    enum ModOrdTypeForcing {NoForcing, PosWDegTO, WDegTOPos, WDegPosTO};


    ModOrdTypeForcing ModuleOrderType(const FreeModule& M)
    {
      if (IsOrdPosn(ordering(M))) return WDegTOPos;
      if (IsWDegPosnOrd(ordering(M))) return WDegPosTO;
      return PosWDegTO;
    } // ModOrdType

  
    // This is OK for the non-homogeneous case
    // For the homogenous case, PosTo this is inefficient, since
    // the Deg rows in the To part are useless.
    SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM,
                                          ModOrdTypeForcing MOType)
    {
      const ModuleOrdering MTO = ordering(FM);
      const SparsePolyRing OldP=RingOf(FM);
      const long NumOldInds=NumIndets(OldP);
      long GrDim;
      if (MOType==PosWDegTO)
        GrDim=0;// Set simple sugar on
      else
        GrDim=GradingDim(OldP);
      const long NumNewInds=NumOldInds+GrDim+1;

      ConstMatrixView OldOrdOMat = OrdMat(OldP);

      matrix NewOrdMat(NewDenseMat(RingZZ(), NumNewInds, NumNewInds));
      ////std::clog<<"NewOrdMat starts as "<<NewOrdMat<<std::endl;
      if (MOType == NoForcing)  MOType = ModuleOrderType(FM);

      switch (MOType)
      {
      case PosWDegTO:
        // Setting the module component ordering
        SetEntry(NewOrdMat, 0, NumNewInds-1, 1); 	
        // Part common to IsWDegPosnOrd and IsOrdPosn
        // Setting the Grading: the OldGrading		
        for (long i=0; i < GrDim+1; ++i)			
          for (long j=0; j < NumOldInds; ++j)		
            SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
        // Setting the Grading: the NewGrading		
        for (long i=1; i < GrDim+1; ++i)			
          SetEntry(NewOrdMat, i, i+NumOldInds-1, 1); 	
        // Setting the TO ordering
        for (long i=GrDim; i < NumOldInds; ++i)
          for (long j=0; j < NumOldInds; ++j)
            SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
        break;

      case WDegTOPos:
        // Part common to IsWDegPosnOrd and IsOrdPosn
        // Setting the Grading: the OldGrading		
        for (long i=0; i < GrDim; ++i)			
          for (long j=0; j < NumOldInds; ++j)		
            SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
        // Setting the Grading: the NewGrading		
        for (long i=0; i < GrDim; ++i)			
          SetEntry(NewOrdMat, i, i+NumOldInds, 1);
        // Setting the TO	
        for (long i=GrDim; i < NumOldInds; ++i)
          for (long j=0; j < NumOldInds; ++j)
            SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
        // Setting the module component ordering
        SetEntry(NewOrdMat, NumNewInds-1-GrDim, NumNewInds-1, 1); 	
        break;

      case WDegPosTO:; // This is the default
      default:
        // Part common to IsWDegPosnOrd and IsOrdPosn
        // Setting the Grading: the OldGrading		
        for (long i=0; i < GrDim; ++i)			
          for (long j=0; j < NumOldInds; ++j)		
            SetEntry(NewOrdMat, i, j, OldOrdOMat(i,j));
      // Setting the Grading: the NewGrading		
      for (long i=0; i < GrDim; ++i)			
        SetEntry(NewOrdMat, i, i+NumOldInds, 1); 	
      // Setting the module component ordering
      SetEntry(NewOrdMat, GrDim, NumNewInds-1, 1);
      // Setting the TO ordering
      for (long i=GrDim; i < NumOldInds; ++i)
        for (long j=0; j < NumOldInds; ++j)
          SetEntry(NewOrdMat, i+1, j, OldOrdOMat(i,j));
      break;
      }
      // Filling the matrix
      for (long i=0; i < GrDim; ++i)
        for (long j=0; j < NumOldInds; ++j)
          SetEntry(NewOrdMat, NumOldInds+i+1, j, OldOrdOMat(i,j));

      const PPOrdering MatNewOrd = NewMatrixOrdering(NewOrdMat, GrDim);

      const std::vector<symbol> IndetNames = NewSymbols(NumOldInds + GrDim + 1);
      //---> for DEBUGGING choose these IndetNames:
      //   std::vector<symbol> IndetNames = SymbolRange("x", 0, NumOldInds-1);
      //   if ( GrDim==1 ) IndetNames.push_back(symbol("s"));  // indet representing shift
      //   else
      //     for ( long i=0 ; i<GrDim ; ++i )
      //       IndetNames.push_back(symbol("s",i));  // indet representing shift
      //   IndetNames.push_back(symbol("e"));  // indet representing module component
      //---> for DEBUGGING
      SparsePolyRing NewP(NewPolyRing(CoeffRing(OldP),IndetNames,MatNewOrd));
      return NewP;
    } // MakeNewPRingFromModule


    // Called by syz, indirectly by intersection, colonbyprincipal (ideal)
    SparsePolyRing MakeNewPRingForP2(const SparsePolyRing& OldP,
                                     ModOrdTypeForcing MOType)
    {
      std::vector<degree> InputShifts;
      degree tmp(GradingDim(OldP));
      InputShifts.push_back(tmp);
      InputShifts.push_back(tmp);
      const FreeModule FM=NewFreeModule(OldP, InputShifts, WDegPosnOrd);
      return MakeNewPRingFromModule(FM, MOType);
    }


    // returns the poly ring equivalent with OldP^2, same grading
    // The ordering is WDegPosnOrd if MOType==NoForcing or MOType
    // Called by intersection, colonbyprincipal
    SparsePolyRing MakeNewPRingForP2_PosTO(const SparsePolyRing& OldP,
                                           bool HomogInput)
    {
      if (HomogInput) return MakeNewPRingForP2(OldP, WDegPosTO);
      else return MakeNewPRingForP2(OldP, PosWDegTO);
    }

    // Called by ComputeGBasis2 (module), ComputeSyz (module), ComputeColonByPrincipal (module)
    SparsePolyRing MakeNewPRingFromModule(const FreeModule& FM)
    { return MakeNewPRingFromModule(FM, NoForcing); }


    // Called by ComputeIntersection (module)
    SparsePolyRing MakeNewPRingFromModule_PosTO(const FreeModule& FM,
                                                bool HomogInput)
    {
      if (HomogInput) return MakeNewPRingFromModule(FM, WDegPosTO);
      else            return MakeNewPRingFromModule(FM, PosWDegTO);
    }


    ////////////////////////////////////////////
    // Embed....
    // from TmpGReductor 2026-04-10
    ////////////////////////////////////


    //GPoly EmbedPoly(ConstRefRingElem the_p, const GRingInfo& theGRI,
    //                const degree& the_d, const long CompIndex)
    //{
    //const RingHom& phi=theGRI.myOldP2NewP();
    //return GPoly(phi(the_p)*theGRI.myE(CompIndex)*theGRI.myY(the_d), theGRI);
    //}


    GPoly EYd(const GRingInfo& theGRI, const long CompIndex, const degree& d)
    {  return GPoly(theGRI.myE(CompIndex)*theGRI.myY(d), theGRI);  }


    GPoly EmbedPoly(ConstRefRingElem p, const GRingInfo& theGRI,
                    const long CompIndex)
    {
      const RingHom& phi=theGRI.myOldP2NewP();
      return GPoly(phi(p)*theGRI.myE(CompIndex), theGRI);
    }


    GPolyList EmbedPolyList(const std::vector<RingElem>& F,
                            const GRingInfo& GRI,
                            const long CompIndex)
    {
      GPolyList result;
      if (F.empty())  return result;
      for (const RingElem& f: F) result.push_back(EmbedPoly(f, GRI, CompIndex));
      return result;
    }


    GPolyList EmbedPolyListNo0(const std::vector<RingElem>& F,
                               const GRingInfo& GRI,
                               const long CompIndex)
    {
      GPolyList result;
      if (F.empty())  return result;
      for (const RingElem& f: F)
        if (!IsZero(f)) result.push_back(EmbedPoly(f, GRI, CompIndex));
      return result;
    }

    /*
      This realizes the embedding FM-->NewP of a vector v.
      e_i->EY[i] and OldP2NewP gives the RingHom between BaseRing(FM) and NewP
    */
    GPoly EmbedVector(const ModuleElem& v,
                      const GRingInfo& theGRI,
                      const long FromCompt)			
    {
      RingElem p(theGRI.myNewSPR());
      const RingHom& phi=theGRI.myOldP2NewP();
      for (long i=0; i<NumCompts(owner(v)); ++i)
        p += phi(v[i]) * theGRI.myEY(i+FromCompt);
      return GPoly(p, theGRI);
    } // EmbedVector


    // GPoly EmbedVector(const ModuleElem& v,
    //                   const GRingInfo& theGRI)			
    // { return EmbedVector(v, theGRI, 0); }


    GPolyList EmbedVectorList(const std::vector<ModuleElem>& VL,
                              const GRingInfo& theGRI,
                              const long FromCompt)
    {
      GPolyList outPL;
      if (VL.empty())  return outPL;
      for (const auto& v: VL)
        if (!IsZero(v))
          outPL.push_back(EmbedVector(v, theGRI, FromCompt));
      return outPL;
    }


    // GPolyList EmbedVectorList(const std::vector<ModuleElem>& VL, const GRingInfo& GRI)
    // { return EmbedVectorList(VL, GRI, 0); }


    GPolyList SyzEmbedVectorList(const std::vector<ModuleElem>& InputVectorList,
                                 const GRingInfo& GRI)
    {
      GPolyList outPL;
      if (InputVectorList.empty())  return outPL;
      const SparsePolyRing NewP=GRI.myNewSPR();
      outPL = EmbedVectorList(InputVectorList, GRI, 0);
      long k=NumCompts(GRI.myFreeModule());
      if (GRI.myInputAndGrading()==NONHOMOG_GRADING)
        for (GPoly& p: outPL)
        { // Added by JAA 2012/10/11
          RingElem Ek = GRI.myE(k++); // previous k
          p.myAppendClear(Ek);
        }
      else
        for (GPoly& p: outPL)
        { // Added by JAA 2012/10/11
          RingElem EkY = GRI.myE(k++) * GRI.myY(wdeg(p)); // previous k
          p.myAppendClear(EkY);
        }
      return outPL;
    } // SyzEmbedVectorList


    // Rewritten by AMB 2026/04/14
    // GPolyList SyzEmbedPolyList(const std::vector<RingElem>& F,
    //                            const GRingInfo& theGRI)
    // {
    //   GPolyList outPL;
    //   if (F.empty())  return outPL;
    //   const SparsePolyRing NewP=theGRI.myNewSPR();
    //   outPL = EmbedPolyList(F, theGRI, 0);
    //   RingElem SyzPP(NewP); // Gives the right degree to p+E^i
    //   degree d(GradingDim(NewP));
    //   long k=1;
    //   for (GPolyList::iterator it=outPL.begin();it!=outPL.end();++it,++k) // for (auto& g: outPL) BUT MUST ALSO INCREMENT k !!
    //   {
    //     SyzPP=one(NewP);
    //     d=wdeg(*it);
    //     for (long j=0; j < GradingDim(NewP); ++j)
    //       SyzPP*=power(theGRI.myY(j),d[j]);
    //     if (theGRI.myInputAndGrading()==NONHOMOG_GRADING)
    //     { // Added by JAA 2012/10/11
    //       RingElem Ek = theGRI.myE(k);
    //       (*it).myAppendClear(Ek);
    //     }
    //     else
    //     { // Added by JAA 2012/10/11
    //       RingElem EkSyzPP = theGRI.myE(k)*SyzPP;
    //       (*it).myAppendClear(EkSyzPP);
    //     }
    //   }
    //   return outPL;
    // } // SyzEmbedPolyList


    GPolyList SyzEmbedPolyList(const std::vector<RingElem>& F,
                               const GRingInfo& theGRI)
    {
      GPolyList outPL;
      if (F.empty())  return outPL;
      const SparsePolyRing NewP=theGRI.myNewSPR();
      outPL = EmbedPolyList(F, theGRI, 0);
      RingElem Ek;
      long k=1;
      for (GPolyList::iterator it=outPL.begin();it!=outPL.end();++it,++k) // for (auto& g: outPL) BUT MUST ALSO INCREMENT k !!
      {
        Ek = theGRI.myE(k);
        if (theGRI.myInputAndGrading() != NONHOMOG_GRADING)
          Ek *= theGRI.myY(wdeg(*it)); // Gives the right degree to p+E^k
        (*it).myAppendClear(Ek);
      }
      return outPL;
    }


    GPolyList IntEmbedPolyLists(const std::vector<RingElem>& F1,
                                const std::vector<RingElem>& F2,
                                const GRingInfo& GRI)
    {
      GPolyList F1_gp = EmbedPolyListNo0(F1, GRI, 0);
      GPolyList F2_gp = EmbedPolyListNo0(F2, GRI, 0);
      GPolyList F2_gp1 = EmbedPolyListNo0(F2, GRI, 1);
      for (auto it=F2_gp.begin(), it1=F2_gp1.begin(); it!=F2_gp.end(); ++it,++it1)
        (*it).myAppendClear(*it1);
      F1_gp.splice(F1_gp.end(), F2_gp);
      return F1_gp;
    }


    GPolyList IntEmbedVectorLists(const std::vector<ModuleElem>& G1,
                                  const std::vector<ModuleElem>& G2,
                                  const GRingInfo& theGRI)
    {
      const long NC = NumCompts(owner(G1[0]));
      GPolyList G1_gp = EmbedVectorList(G1, theGRI, 0);
      GPolyList G2_gp = EmbedVectorList(G2, theGRI, 0);
      GPolyList G2_gpNC  = EmbedVectorList(G2, theGRI, NC);
      for (auto it=G2_gp.begin(), itNC=G2_gpNC.begin(); it!=G2_gp.end(); ++it,++itNC)
        (*it).myAppendClear(*itNC);
      G1_gp.splice(G1_gp.end(), G2_gp);
      return G1_gp;
    }


    GPolyList ColonEmbedVectorLists(const std::vector<ModuleElem>& VL,
                                    const ModuleElem& v,
                                    const GRingInfo& theGRI)
    {
      GPolyList FirstPart;
      if (VL.empty())  return FirstPart;
      FirstPart = EmbedVectorList(VL, theGRI, 0);
      GPoly GP1 = EmbedVector(v, theGRI, 0);
      GPoly GP2 = EYd(theGRI, NumCompts(owner(v)), wdeg(v));
      GP1.myAppendClear(GP2);
      FirstPart.push_back(GP1);
      return FirstPart;
    }


    GPolyList ColonEmbedPolyLists(const std::vector<RingElem>& G,
                                  const RingElem& f,
                                  const GRingInfo& theGRI)
    {
      GPolyList FirstPart;
      if (G.empty())  return FirstPart;
      FirstPart = EmbedPolyList(G, theGRI, 0);
      GPoly GP1 = EmbedPoly(f, theGRI, 0);
      GPoly GP2 = EYd(theGRI, 1, wdeg(f));
      GP1.myAppendClear(GP2);
      FirstPart.push_back(GP1);
      return FirstPart;
    }


    // void ColonEmbedGPolyList(GPolyList& theGPL, GPoly& the_gp)
    // {
    //   const GRingInfo& GRI(the_gp.myGRingInfo());
    //   CoCoA_ASSERT(theGPL.begin()->myGRingInfo()==GRI);
    //   const long NC = NumCompts(GRI.myFreeModule());
    //   RingElem tmp =  GRI.myE(NC)*GRI.myY(wdeg(the_gp)); // JAA 2012-10-11
    //   the_gp.myAppendClear(tmp);                         // JAA 2012-10-11
    //   theGPL.push_back(GPoly(GRI));
    //   theGPL.back().AssignClear(the_gp);
    // } // ColonEmbedGPolyList


    ////////////////////////////////////////////
    // DeEmbed....
    // from TmpGReductor 2026-04-10
    ////////////////////////////////////

    // identical copy in TmpGReductor:
    ModuleElem DeEmbedPoly(ConstRefRingElem g,
                           const GRingInfo& theGRI,
                           const long FromCompt)
    {
      const SparsePolyRing OldP=theGRI.myOldSPR();
      const SparsePolyRing NewP=theGRI.myNewSPR();
      const FreeModule FM=theGRI.myOutputFreeModule();
      ModuleElem v(FM);
      
      const std::vector<ModuleElem>& e = gens(FM);
      
      RingElem tmp(OldP);
      for (SparsePolyIter i=BeginIter(g); !IsEnded(i); ++i)
      {
        tmp = theGRI.myNewP2OldP()(monomial(NewP,coeff(i),PP(i)));
        CoCoA_ASSERT(theGRI.myCompt_orig(PP(i)) - FromCompt >= 0);
        CoCoA_ASSERT(theGRI.myCompt_orig(PP(i)) - FromCompt < NumCompts(FM));
        v += tmp * e[theGRI.myCompt_orig(PP(i)) - FromCompt];
      }
      return v;
    }


    // Polys whose LPP has compt_work > max-FromCompt,
    // i.e. compt_orig < FromCompt,  are ignored
    std::vector<ModuleElem> DeEmbedPolyList(const std::vector<RingElem>& G,
                                            const GRingInfo& theGRI,
                                            const long FromCompt)
    {
      std::vector<ModuleElem> G_out;
      if (G.empty())  return G_out;
      G_out.reserve(len(G));
      for (const RingElem& g: G) // vvvvvvvvvv  only for PosTO 
        if (theGRI.myCompt_work(LPP(g)) <= theGRI.myCompt_OrigToWork(FromCompt))
          G_out.push_back(DeEmbedPoly(g, theGRI, FromCompt));
      return G_out;
    }


    std::vector<RingElem> DeEmbedPolyListToPL(const std::vector<RingElem>& G,
                                              const GRingInfo& theGRI,
                                              const long FromCompt)
    {
      VerboseLog VERBOSE("DeEmbedPolyListToPL");
      std::vector<RingElem> G_out;
      if (G.empty())  return G_out;
      G_out.reserve(len(G));
      for (const RingElem& g: G)
        if (theGRI.myCompt_work(LPP(g)) <= theGRI.myCompt_OrigToWork(FromCompt))
        {
          VERBOSE(100) << "LPP(g) = " << LPP(g) << std::endl; /////////////////
          G_out.push_back(theGRI.myNewP2OldP()(g));
          // if ( IsConstant(outPL.last()) ) // redmine #1647 ////' doesn't happen
          // {
          //   outPL.clear();
          //   outPL.push_back(theGRI.myOldSPR()->myOne());
          //   break;
          // }
        }
      return G_out;
    }

    ////////////////////////////////////////////
    // homog....
    // from TmpGReductor 2026-04-10
    ////////////////////////////////////

    // The ordering is supposed to be Deg compatible
    RingElem homog(ConstRefRingElem the_p, const std::vector<RingElem>& theY)
    {
      const SparsePolyRing SPR=owner(the_p);
      if (GradingDim(SPR) != len(theY))
        CoCoA_THROW_ERROR2(ERR::BadArg, "incompatible GradingDim");
      RingElem the_hp(SPR);
      RingElem tmp(SPR);
      degree MaxDeg(GradingDim(SPR));
      degree TmpDeg(GradingDim(SPR));
      // Compute the maximum degree
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it);++it)
        MaxDeg=top(MaxDeg,wdeg(PP(it)));
      // Homogenizing
      for (SparsePolyIter it=BeginIter(the_p);!IsEnded(it);++it)
      {
        //std::clog<<"Homogenized:tmp"<<tmp<<endl;
        tmp = monomial(SPR,coeff(it),PP(it));
        TmpDeg=wdeg(PP(it));
        for (long i=0; i != len(theY); ++i)
          tmp*=power(theY[i],(MaxDeg[i]-TmpDeg[i]));
        SPR->myAddClear(raw(the_hp), raw(tmp));
      }
      return the_hp;
    } // homog (multihomog)


    void homogenized(ModuleElem& /*the_hv*/,
                     const ModuleElem& /*the_v*/,
                     const GRingInfo& /*theGRI*/)
    {
      CoCoA_THROW_ERROR1(ERR::NYI);
    } // homogenized


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
    const FreeModule FM=owner(G_in[0]);
    const SparsePolyRing P_work(MakeNewPRingFromModule(FM));
    const SparsePolyRing P(RingOf(FM));
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool IsSatAlg=false;
    const GRingInfo GRI(P_work,P,FM,FM,IsHomogGrD0(G_in),
                        IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GPolyList EmbeddedPolys = EmbedVectorList(G_in, GRI, 0); // removes 0s
    if (EmbeddedPolys.empty())
    {
      GB_out.clear();
      MinGens_out.clear();
      return;
    }
    GReductor GBR(GRI, EmbeddedPolys);
    GBR.myDoGBasis();// homog input standard alg interred
    VectorList GB_tmp = GBR.myExportGBasis_module();
    VectorList MinGens_tmp;
    if (GradingDim(P)>0 && IsHomog(G_in))
      MinGens_tmp = GBR.myExportMinGens_module();
    std::swap(GB_out, GB_tmp);
    std::swap(MinGens_out, MinGens_tmp);
  } // ComputeGBasis2


  void ComputeGBasis2(std::vector<RingElem>& GB_out, std::vector<RingElem>& MinGens_out,
                      const std::vector<RingElem>& G_in, const CpuTimeLimit& CheckForTimeout)
  {
    if (G_in.empty())
    {
      GB_out.clear();
      MinGens_out.clear();
      return;
    }
    SparsePolyRing P(owner(G_in[0]));
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool IsSatAlg=false;
    std::vector<RingElem> GB_tmp;
    std::vector<RingElem> MinGens_tmp;
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    {
      const ring R = BaseRing(CoeffRing(P));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G_in),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G_in, Rx));
      GBR.myDoGBasis();// homog input standard alg interred
      GB_tmp = monic(RxToKx(GBR.myExportGBasis(), P)); // 2016-11-22: make monic
      if (GradingDim(P)>0 && IsHomog(G_in))
        MinGens_tmp = RxToKx(GBR.myExportMinGens(), P);
    }
    else
    {
      GRingInfo GRI(P,IsHomogGrD0(G_in),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G_in);
      GBR.myDoGBasis();// homog input standard alg interred
      GB_tmp = GBR.myExportGBasis();
      if (GradingDim(P)>0 && IsHomog(G_in))
        MinGens_tmp = GBR.myExportMinGens();
    }
    swap(GB_out, GB_tmp);
    swap(MinGens_out, MinGens_tmp);
  } // ComputeGBasis2
  

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
    const SparsePolyRing P(owner(G_in[0]));
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (GradingDim(P)!=1)  CoCoA_THROW_ERROR1(ERR::ReqGradingDim1);
    if (!IsHomog(G_in))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    if (TruncDeg < 0)  CoCoA_THROW_ERROR2(ERR::ReqNonNegative, "user TruncDeg");
    bool IsSatAlg=false;
    PolyList GB_tmp;
    PolyList MinGens_tmp;
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    { //---------------------------------------------------
      const SparsePolyRing Rx = NewPolyRing(BaseRing(CoeffRing(P)), symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G_in), IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G_in, Rx));
      GBR.mySetTruncDeg(TruncDeg); // input value
      //---------------------------------------------------
      GBR.myDoGBasis();// homog input standard alg interred
      GB_tmp = monic(RxToKx(GBR.myExportGBasis(), P)); // 2016-11-22: make monic
      if (IsEveryWDegLtEq(G_in, TruncDeg))
        MinGens_tmp = RxToKx(GBR.myExportMinGens(), P);
      TruncDeg = GBR.myTruncDeg(); // update TruncDeg (if changed/complete)
    }
    else
    { //---------------------------------------------------
      GRingInfo GRI(P, IsHomogGrD0(G_in), IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G_in);
      GBR.mySetTruncDeg(TruncDeg); // input value
      //---------------------------------------------------
      GBR.myDoGBasis();// homog input standard alg interred
      GB_tmp = GBR.myExportGBasis();
      if (IsEveryWDegLtEq(G_in, TruncDeg))
        MinGens_tmp = GBR.myExportMinGens();
      TruncDeg = GBR.myTruncDeg(); // update TruncDeg (if changed/complete)
    }
    swap(GB_out, GB_tmp);
    swap(MinGens_out, MinGens_tmp);
  } // ComputeGBasisTrunc2


  PolyList ComputeGBasis(const PolyList& G, const CpuTimeLimit& CheckForTimeout)
  {
    PolyList GB;
    PolyList MinGens; // useless
    ComputeGBasis2(GB, MinGens, G, CheckForTimeout);
    return GB;
  }


  VectorList ComputeGBasis(const VectorList& G, const CpuTimeLimit& CheckForTimeout)
  {
    VectorList GB;
    VectorList MinGens; // useless
    ComputeGBasis2(GB, MinGens, G, CheckForTimeout);
    return GB;
  }


  PolyList ComputeGBasisSelfSatCore(const PolyList& G, const CpuTimeLimit& CheckForTimeout)
  {
    if (G.empty()) return G;
    //    bool IsSatAlg=true;
    SparsePolyRing P(owner(G[0]));
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
    {
      const ring R = BaseRing(CoeffRing(P));
      SparsePolyRing Rx = NewPolyRing(R, symbols(PPM(P)), ordering(PPM(P)));
      GRingInfo GRI(Rx, IsHomogGrD0(G), true/*IsSatAlg*/, NewDivMaskEvenPowers(), CheckForTimeout);
      GRI.mySetCoeffRingType(CoeffEncoding::FrFldOfGCDDomain);
      GReductor GBR(GRI, KxToRx(G, Rx), GReductor::SaturatingAlg);
      GBR.myDoGBasisSelfSatCore();// homog input standard algorithm interred
      return monic(RxToKx(GBR.myExportGBasis(), P));
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
    SparsePolyRing P(owner(G[0]));
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
      return monic(RxToKx(GBR.myExportGBasis(), P)); // 2016-11-22: make monic
    }
    else
    {
      GRingInfo GRI(P,IsHomogGrD0(G),false/*IsSatAlg*/,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, G);
      GBR.myDoGBasisRealSolve();
      return GBR.myExportGBasis();
    }          
  }


  namespace
  { // namespace // anonymous ----------------------------------------------

    std::vector<long> PPIndices(ConstRefPPMonoidElem t)
    {
      std::vector<long> indices;
      std::vector<long> tmp = exponents(t);
      indices.reserve(len(tmp));
      for (long i=0; i < len(tmp); ++i)
        if (tmp[i]!=0)  indices.push_back(i);
      return indices;
    }


    SparsePolyRing MakeElimRing(const SparsePolyRing& P,
                                const std::vector<long>& IndetsToElim,
                                const bool IsHomog)
    {
      ConstMatrixView weights = GradingMat(P);
      matrix NewOrdMat = ElimMat(IndetsToElim, weights);
      long NewGrDim = 0;
      if (GradingDim(P) != 0 && IsHomog)
      {
        NewOrdMat = ElimHomogMat(IndetsToElim, weights);
        NewGrDim = GradingDim(P);
      }
      const PPOrdering& NewOrd = NewMatrixOrdering(NewOrdMat, NewGrDim);
      std::vector<symbol> IndetNames = NewSymbols(NumIndets(P));
      return NewPolyRing(CoeffRing(P), IndetNames, NewOrd);
    }

  } // namespace // anonymous ----------------------------------------------


  PolyList ComputeElim(const PolyList& G, ConstRefPPMonoidElem inds,
                       const CpuTimeLimit& CheckForTimeout)
  {
    VerboseLog VERBOSE("ComputeElim");
    VERBOSE(99) << "-- called --" << std::endl;
    if (G.empty())  return G;
    const SparsePolyRing P = owner(G[0]);
    bool HomogInput = (GradingDim(P)!=0 && IsHomog(G));
    SparsePolyRing P_work=MakeElimRing(P, PPIndices(inds), HomogInput);
    RingHom PToNew = PolyAlgebraHom(P, P_work, indets(P_work));
    RingHom NewToP = PolyAlgebraHom(P_work, P, indets(P));
    std::vector<RingElem> NewGens;  NewGens.reserve(len(G));
    for (const RingElem& g: G) NewGens.push_back(PToNew(g));
    PPMonoidElem ElimIndsProd(LPP(PToNew(monomial(P,inds))));
    std::vector<RingElem> GB = ComputeGBasis(NewGens, CheckForTimeout);
    std::vector<RingElem> ElimGens;  ElimGens.reserve(len(GB));
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
  } // ComputeElim


  RingElem ComputeElimFirst(const PolyList& G, ConstRefPPMonoidElem inds,
                            const CpuTimeLimit& CheckForTimeout)
  {
    const SparsePolyRing P = owner(G[0]);
    if (!IsField(CoeffRing(P)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    bool HomogInput = (GradingDim(P)!=0 && IsHomogGrD0(G));
    if (!HomogInput)  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    bool IsSatAlg = false;
    SparsePolyRing P_work = MakeElimRing(P, PPIndices(inds), HomogInput);
    RingHom OldToNew = PolyAlgebraHom(P, P_work, indets(P_work));
    RingHom NewToOld = PolyAlgebraHom(P_work, P, indets(P));
    PPMonoidElem ElimIndsProd = LPP(OldToNew(monomial(P,inds)));
    if (IsFractionFieldOfGCDDomain(CoeffRing(P)))
      CoCoA_THROW_ERROR2(ERR::NYI, "Only for FFp");
    else
    {
      GRingInfo GRI(P_work,IsHomogGrD0(G),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
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


  // ------- ComputeSyz --------------------------------

  VectorList ComputeSyz(const FreeModule& SyzFM, const VectorList& G,
                        const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsHomogGrD0(G))
      CoCoA_THROW_ERROR2(ERR::NYI, "Not Yet Tested for non-homog input");
    if (G.empty())  return G;
    bool IsSatAlg=false;
    ModOrdTypeForcing MOType = (IsHomogGrD0(G) ? WDegPosTO : PosWDegTO);
    const FreeModule FM = owner(G[0]);
    const SparsePolyRing P(RingOf(FM));
    const SparsePolyRing P_work(MakeNewPRingFromModule(FM,MOType));
    // Note: the GRI should build itself SyzFM and P_work from the data and deduce FM and P.
    //       All the embedding/deembedding functions should be members of GRI.
    GRingInfo GRI(P_work,P,FM,SyzFM,IsHomogGrD0(G),IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, SyzEmbedVectorList(G,GRI));
    GBR.myDoGBasis();
    return DeEmbedPolyList(GBR.myExportGBasis(), GRI, NumCompts(FM));
  }


  VectorList ComputeSyz(const FreeModule& SyzFM, const PolyList& G,
                        const CpuTimeLimit& CheckForTimeout)
  {
    if (G.empty()) return VectorList{zero(SyzFM)};
    bool IsSatAlg=false;
    ModOrdTypeForcing MOType = (IsHomogGrD0(G) ? WDegPosTO : PosWDegTO);
    const SparsePolyRing P(owner(G[0]));
    const SparsePolyRing P_work(MakeNewPRingForP2(P,MOType));
    GRingInfo GRI(P_work, P, SyzFM, IsHomogGrD0(G),
                  IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, SyzEmbedPolyList(G,GRI));
    GBR.myDoGBasis();
    return DeEmbedPolyList(GBR.myExportGBasis(), GRI, 1);
  }


  // ------- ComputeIntersection --------------------------------

  VectorList ComputeIntersection(const VectorList& G1, const VectorList& G2,
                                 const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    const FreeModule FM = owner(G1[0]);
    bool HomogInput = IsHomogGrD0(G1) && IsHomogGrD0(G2);
    const SparsePolyRing P_work(MakeNewPRingFromModule_PosTO(FM,HomogInput));
    GRingInfo GRI(P_work, RingOf(FM), FM, FM, HomogInput,
                  IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, IntEmbedVectorLists(G1, G2, GRI));
    GBR.myDoGBasis();// homog input standard alg interred
    return DeEmbedPolyList(GBR.myExportGBasis(), GRI, NumCompts(FM));
  }


  PolyList ComputeIntersection(const PolyList& G1, const PolyList& G2,
                               const CpuTimeLimit& CheckForTimeout)
  {
    bool IsSatAlg=false;
    if (G1.empty())  return G1;
    const SparsePolyRing P(owner(G1[0]));
    const bool HomogInput = IsHomogGrD0(G1) && IsHomogGrD0(G2);
    const SparsePolyRing P_work(MakeNewPRingForP2_PosTO(P, HomogInput));
    GRingInfo GRI(P_work, P, HomogInput,
                  IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
    GReductor GBR(GRI, IntEmbedPolyLists(G1, G2, GRI));
    GBR.myDoGBasis();// homog input standard alg interred
    return DeEmbedPolyListToPL(GBR.myExportGBasis(), GRI, 1);
  }


  // ------- ComputeColon --------------------------------

  namespace
  { // namespace // anonymous ----------------------------

    PolyList ComputeColonByPrincipal(const VectorList& G1, const ModuleElem& v,
                                     const CpuTimeLimit& CheckForTimeout)
    {
      CoCoA_ASSERT(!G1.empty()); // guaranteed by ComputeColon
      bool IsSatAlg=false;
      if (!IsHomogGrD0(G1))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
      if (!IsHomogGrD0(v))  CoCoA_THROW_ERROR2(ERR::NYI, "non-homog");
      const FreeModule FM = owner(G1[0]);
      const SparsePolyRing P_work(MakeNewPRingFromModule(FM, WDegPosTO));
      GRingInfo GRI(P_work, RingOf(FM), FM, FM, IsHomogGrD0(G1)&&IsHomogGrD0(v),
                    IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, ColonEmbedVectorLists(G1, v, GRI));
      GBR.myDoGBasis();// homog input standard alg interred
      return DeEmbedPolyListToPL(GBR.myExportGBasis(), GRI, NumCompts(FM));
    }
    

    PolyList ComputeColonByPrincipal(const PolyList& G, ConstRefRingElem f,
                                     const CpuTimeLimit& CheckForTimeout)
    {
      CoCoA_ASSERT(!G.empty()); // guaranteed by ComputeColon
      PolyList tmpColonResult;
      if (IsZero(f))  return std::vector<RingElem>{one(owner(f))};
      bool IsSatAlg=false;
      const SparsePolyRing P(owner(G[0]));
      const SparsePolyRing P_work(MakeNewPRingForP2_PosTO(P,IsHomogGrD0(G)&&IsHomogGrD0(f)));
      GRingInfo GRI(P_work, P, IsHomogGrD0(G) && IsHomogGrD0(f),
                    IsSatAlg,NewDivMaskEvenPowers(), CheckForTimeout);
      GReductor GBR(GRI, ColonEmbedPolyLists(G, f, GRI));
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
    if (G2.empty())  return std::vector<RingElem>{one(RingOf(owner(G1[0])))};
    if (G1.empty())  return std::vector<RingElem>{};
    PolyList aux;
    VectorList::const_iterator it=G2.begin();
    PolyList tmpColon = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
    for (++it; it!=G2.end(); ++it)
    {
      aux = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
      tmpColon = ComputeIntersection(tmpColon, aux, CheckForTimeout);
      CoCoA_ASSERT_ALWAYS(!tmpColon.empty()); // just paranoia -- remove ALWAYS
      if (tmpColon.empty()) break;
    }
    return tmpColon;
  }


  PolyList ComputeColon(const PolyList& G1, const PolyList& G2,
                        const CpuTimeLimit& CheckForTimeout)
  {
    if (G1.empty() && G2.empty())
      CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "both lists are empty");
    if (G2.empty()) return std::vector<RingElem>{one(owner(G1[0]))};
    if (G1.empty()) return G1;
    PolyList aux;
    PolyList::const_iterator it=G2.begin();
    PolyList tmpColon = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
    for (++it; it!=G2.end(); ++it)
      if (!IsZero(*it))
      {
        aux = ComputeColonByPrincipal(G1, *it, CheckForTimeout);
        tmpColon = ComputeIntersection(tmpColon, aux, CheckForTimeout);
        CoCoA_ASSERT_ALWAYS(!tmpColon.empty()); // just paranoia -- remove ALWAYS
      }
    return tmpColon;
  }


  // ------- ComputeSaturation --------------------------------

  namespace
  { // namespace // anonymous ---------------------------------
    // These procedures are for EqualLPPs used by Saturation
  
    // PolyList::const_iterators are ordered according to LPP of their polys
    bool ByLPP(const PolyList::const_iterator& it,
               const PolyList::const_iterator& it2)
    {
      return (SparsePolyRingPtr(owner(*it))->myCmpLPP(raw(*it),raw(*it2)) == -1);
    }
    
    
    // Useful when you have I\subset J and you want to check I==J
    bool AreEqualLPPs(const PolyList& G1, const PolyList& G2)
    {
      if (len(G1) != len(G2))  return false;
      std::vector<PolyList::const_iterator> v1,v2;
      v1.reserve(len(G1));
      v2.reserve(len(G2));
      for (auto it=G1.begin(); it!=G1.end(); ++it)
        v1.push_back(it);
      for (auto it=G2.begin(); it!=G2.end(); ++it)
        v2.push_back(it);
      stable_sort(v1.begin(), v1.end(), ByLPP);
      stable_sort(v2.begin(), v2.end(), ByLPP);
      const long len1 = len(v1);
      for (int i=0; i!=len1; ++i)
        if (LPP(*(v1[i])) != LPP(*(v2[i])))  return false;
      return true;
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
        PushBack(ans, coeff(it), PP(it)/t);
      return ans;
    }


    PolyList ComputeSaturationHomogByPP(const PolyList& G,
                                        ConstRefPPMonoidElem t,
                                        const CpuTimeLimit& CheckForTimeout)
    {
      VerboseLog VERBOSE("ComputeSaturationHomogByPP");
      VERBOSE(99) << "-- called --" << std::endl;
      //    if (IsZero(I))  return ideal(zero(Ph));
      const ring& P = owner(G[0]);
      const ring& K = CoeffRing(P); // CoCoA_ASSERT(IsField(k));
      const long GrDim = GradingDim(P);
      const long NumX = NumIndets(P);
      std::vector<symbol> names;  names.reserve(NumX);
      for (long i=0; i<NumX; ++i)  names.push_back(IndetSymbol(P,i));
      VERBOSE(90) << "GrDim = " << GrDim << "   W = " << GradingMat(P) << std::endl;
      matrix W0 = NewDenseMat(ConcatVer(GradingMat(P), ZeroMat(RingZZ(), 1, NumIndets(P))));
      PolyList tmpGens = G;
      std::vector<long> expv = exponents(t);
      for (long i=0; i<NumIndets(P); ++i)
        if (expv[i] != 0)
        {
          VERBOSE(95) << "doing indet " << i << " = " << indet(P,i) << std::endl;
          SetEntry(W0, GrDim,i, -1); // W0[GrDim,i] = -1
          ring PzDegRevLex = NewPolyRing(K, names, NewMatrixOrdering(MakeTermOrdMat(W0), GrDim));
          SetEntry(W0, GrDim,i, 0); // back to W0[GrDim,i] = 0
          RingHom phi = PolyAlgebraHom(owner(tmpGens[0]), PzDegRevLex, indets(PzDegRevLex));
          tmpGens = ComputeGBasis(phi(tmpGens), CheckForTimeout);
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
      const SparsePolyRing P = owner(G[0]);
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
    if (G2.empty())  return std::vector<RingElem>{one(owner(G1[0]))};
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


  // ------- ComputeHomogenization --------------------------------

  VectorList ComputeHomogenization(const VectorList& G, const PolyList& /*theIndets*/,
                             const CpuTimeLimit& /*CheckForTimeout*/)
  {
    CoCoA_THROW_ERROR2(ERR::NYI, "for modules");
    return G;  // just to keep the compiler quiet
  }


  std::vector<RingElem> ComputeHomogenization(const std::vector<RingElem>& G,
                                              const std::vector<RingElem>& HomogIndets,
                                              const CpuTimeLimit& CheckForTimeout)
  {
    if (G.empty()) return G;
    RingElem IndetProd(product(HomogIndets));
    std::vector<RingElem> tmpHGens;  tmpHGens.reserve(len(G));
    for (const RingElem& g: G)
      tmpHGens.push_back(homog(g, HomogIndets));
    return ComputeSaturationByPrincipal(tmpHGens,IndetProd,CheckForTimeout);
  }


  // ------- RadicalMembership --------------------------------

  // WARN: it supposes ComputeSaturationByPrincipal returns a GB
  bool RadicalMembership(const PolyList& G, ConstRefRingElem f,
                         const CpuTimeLimit& CheckForTimeout)
  {
    PolyList tmpGens;
    tmpGens = ComputeSaturationByPrincipal(G, f, CheckForTimeout);
    if (len(tmpGens) != 1) return false;
    return IsInvertible(tmpGens.front());
  }

  
} // end of namespace cocoa
