//   Copyright (c)  2007-2018  John Abbott and Anna M Bigatti

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
#include "TmpPoincareCPP.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/RingDenseUPolyClean.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/VectorOps.H" // for debugging only
#include "CoCoA/degree.H" // for debugging only
#include "CoCoA/ring.H" // for debugging only
#include "IVectors.h"
#include "TermList.h"
#include "eterms.h"
//#include "unipoly.h"
#include "CoCoA/verbose.H"  // for VerboseLog

#include <vector>
using std::vector;

namespace CoCoA
{

  /**********************************************************/
  void poincare_init(const DenseUPolyRing& HSRing, int IndetsNo);
  /**********************************************************/


  void MakeHPPowerList(std::vector<RingElem>& PL, const DenseUPolyRing& HSRing, size_t MaxDeg)
  {
    vector<RingElem> ResPowerList;
    RingElem x = indet(HSRing, 0);
    RingElem f = one(HSRing);
    
    ResPowerList.push_back(f);
    for (size_t d=1; d<=MaxDeg; ++d)
    {
      HSRing->myMulBy1MinusXExp(raw(f), 1);
      ResPowerList.push_back(f); 
    }
    swap(PL, ResPowerList);
  }


  DenseUPolyRing poincare_init2()
  {
    static int initialized = 0;
    if ( !initialized )
    {
      initialized = 1;      
      //    unipoly_rum_init(PoincareMaxPower);
      MakeGlobalHPPowerList(NewPolyRing_DUP(RingZZ())); //AMB 2010-10-25
    }
    return owner(HPPowerList(1));
  }


//   void poincare_init_without_unipoly(int IndetsNo)
//   {
//     static int initialized = 0;
//     if ( !initialized )
//     {
//       initialized = 1;
//       rum_init_all();
//       //------- disable RUM ----------------------------------
//       rum_init(12,2000); /* pointers */
//       rum_init(200*sizeof(eterm), 120);
//       rum_init(sizeof(TermListStruct), RUM_STD_SIZE);
//       //------- disable RUM ----------------------------------
//       GlobalTermList_init(IndetsNo);
//     }
//   }


  DenseUPolyRing StartPoincareQQt(int IndetsNo)
  {
    ::StartPoincare(IndetsNo);
    return poincare_init2();
  }


//   void StartPoincare(int IndetsNo)
//   {
//     poincare_init(IndetsNo);
//     if ( GlobalTermList_size()!=IndetsNo )
//     {
//       GlobalTermList_free();
//       GlobalTermList_init(IndetsNo);
//     }
//   }

#define malloc_TList() (TermList)mymalloc(sizeof(TermListStruct),"TList")
#define free_TList(p) {myfree(sizeof(TermListStruct),p,"TList");}

  namespace
  { // temporary: for multigraded Hilbert
    
    degree wdeg(const PPMonoid& PPM, const eterm t)
    {
      PPMonoidElem pp(PPM);
      for ( long i=1 ; i<=eterm_get_indetsNo(t) ; ++i )
        if ( t[i] != 0 )
          PPM->myMulIndetPower(raw(pp), i-1, t[i]);
      return wdeg(pp);
    }

    RingElem XExp(const SparsePolyRing& P, degree d)
    {
      RingElem f(one(P));
      for ( long i=0 ; i<NumIndets(P) ; ++i )
        if ( d[i] != 0 )
        {
          f *= IndetPower(P, i, d[i]);
        }
      return f;
    }


    RingElem OneMinusXExp(const SparsePolyRing& P, degree d)
    { return 1 - XExp(P, d); }
    
  } // temporary: for multigraded Hilbert


  /****************************************************************
   *************************  BASE CASES  *************************
   ****************************************************************/

  RingElem SPPoincare(const DenseUPolyRing& P, TermList theTList)
  {
    eterm   SPL = SPList(theTList);
    ints OccInd = Indets(SPL);
    int i = IntsGetLen(OccInd);
    
    if ( i==eterm_degree(SPL) )
    {
      EraseAndFreeTList(theTList);
      RingElem res(P);
      CopyHPPower(res, i);
      return res;
    }
    //    else
    RingElem res(one(P));  //UnipolyOne(eterm_degree(SPL));
    for ( ; i>0; --i) P->myMulBy1MinusXExp(raw(res), (unsigned long)SPL[OccInd[i]]);
    EraseAndFreeTList(theTList);
    return res;
  }


  RingElem SPPoincare(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList)
  {
    eterm   SPL = SPList(theTList);
    ints OccInd = Indets(SPL);
    
    RingElem res(one(P));
    for (int i=IntsGetLen(OccInd) ; i>0; --i)
      res *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, SPL[OccInd[i]])));
    EraseAndFreeTList(theTList);
    return res;
  }


  RingElem LenOnePoincare(const DenseUPolyRing& P, TermList theTList)
  {
    eterm    SPL = SPList(theTList), myTerm = (MTList(theTList))[1];
    ints     OccInd = Indets(SPL);
    int   TDeg =eterm_degree(myTerm), i =IntsGetLen(OccInd), exp;
    RingElem p1(P), p2(P);

    if (eterm_coprime(myTerm, SPL))
    {
      if ( i==eterm_degree(SPL) )
        //        p1 = HPPowerList(i);
        CopyHPPower(p1, i);
      else
      {
        p1 = one(P); //UnipolyOne(eterm_degree(SPL)+TDeg);
        for ( ; i>0 ; --i )  P->myMulBy1MinusXExp(raw(p1), SPL[OccInd[i]]);
      }
      P->myMulBy1MinusXExp(raw(p1), TDeg);
      EraseAndFreeTList(theTList);
      return p1;
    }
    p1 = one(P); //UnipolyOne(eterm_degree(SPL)+TDeg);
    p2 = one(P); //UnipolyOne(eterm_degree(SPL));
    for ( ; i>0 ; --i )
    {
      exp = SPL[OccInd[i]];
      P->myMulBy1MinusXExp(raw(p1), exp);
      P->myMulBy1MinusXExp(raw(p2), exp-myTerm[OccInd[i]]);
    }
    EraseAndFreeTList(theTList);
    //P1MinusXExpP2(p1, TDeg, p2);
    P->myAddMulLM(raw(p1), raw(-one(CoeffRing(P))), TDeg, raw(p2));
    return p1;
  }
  
  RingElem OneTermAndSPPoincare(const DenseUPolyRing& P, eterm theTerm, eterm SPL)
  {
    ints  OccInd = Indets(SPL);
    int   TDeg = eterm_degree(theTerm), i = IntsGetLen(OccInd), exp;
    RingElem p1(P), p2(P);

    if (eterm_coprime(theTerm, SPL))
    {
      if ( i==eterm_degree(SPL) )
        //        p1 = HPPowerList(i);
        CopyHPPower(p1, i);
      else
      {
        p1 = one(P); //UnipolyOne(eterm_degree(SPL)+TDeg);
        for ( ; i>0 ; --i )  P->myMulBy1MinusXExp(raw(p1), SPL[OccInd[i]]);
      }
      P->myMulBy1MinusXExp(raw(p1), eterm_degree(theTerm));
      eterm_free(theTerm);
      eterm_free(SPL);
      return p1;
    }
    p1 = one(P); //UnipolyOne(eterm_degree(SPL)+TDeg);
    p2 = one(P); //UnipolyOne(eterm_degree(SPL));
    for ( ; i>0 ; --i )
    {
      exp = SPL[OccInd[i]];
      P->myMulBy1MinusXExp(raw(p1), exp);
      P->myMulBy1MinusXExp(raw(p2), exp-theTerm[OccInd[i]]);
    }
    // P1MinusXExpP2(p1, TDeg, p2);
    P->myAddMulLM(raw(p1), raw(-one(CoeffRing(P))), TDeg, raw(p2));
    eterm_free(theTerm);
    eterm_free(SPL);
    return p1;
  }


  RingElem LenOnePoincare(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList)
  {
    eterm    SPL = SPList(theTList), myTerm = (MTList(theTList))[1];
    ints     OccInd = Indets(SPL);
    int   exp;
    RingElem p1(one(P)), p2(one(P));

    if (eterm_coprime(myTerm, SPL))
    {
      for (int i=IntsGetLen(OccInd) ; i>0; --i)
        p1 *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, SPL[OccInd[i]])));
      p1 *= OneMinusXExp(P, wdeg(PPM, myTerm));
      EraseAndFreeTList(theTList);
      return p1;
    }
    for (int i=IntsGetLen(OccInd) ; i>0; --i)
    {
      exp = SPL[OccInd[i]];
      p1 *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, exp)));
      p2 *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, exp-myTerm[OccInd[i]])));
    }
    p1 += -XExp(P,wdeg(PPM,myTerm)) * p2;
    EraseAndFreeTList(theTList);
    return p1;
  }
  
  RingElem OneTermAndSPPoincare(const SparsePolyRing& P, const PPMonoid& PPM, eterm theTerm, eterm SPL)
  {
    ints  OccInd = Indets(SPL);
    int   exp;
    RingElem p1(one(P)), p2(one(P));

    if (eterm_coprime(theTerm, SPL))
    {
      for (int i=IntsGetLen(OccInd) ; i>0; --i)
        p1 *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, SPL[OccInd[i]])));
      p1 *= OneMinusXExp(P, wdeg(PPM, theTerm));
      eterm_free(theTerm);
      eterm_free(SPL);
      return p1;
    }
    for (int i=IntsGetLen(OccInd) ; i>0; --i)
    {
      exp = SPL[OccInd[i]];
      p1 *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, exp)));
      p2 *= OneMinusXExp(P, wdeg(IndetPower(PPM, OccInd[i]-1, exp-theTerm[OccInd[i]])));
    }
    p1 += -XExp(P, wdeg(PPM, theTerm)) * p2;
    eterm_free(theTerm);
    eterm_free(SPL);
    return p1;
  }

  /*************************************************************
   ************************    SPLIT    ************************
   *************************************************************/

  RingElem TotalSplitPoincare(const DenseUPolyRing& P, TermList theTList)
  {
    RingElem res(P);
    eterm  myTerm, SPL =SPList(theTList), mySPL;
    MixedTermList  MTL =MTList(theTList);
    int MTLLen =MTListLen(theTList), n =MTLLen, DegSum =0;
    VerboseLog VERBOSE("TotalSplitPoincare");

    for ( ; n>0 ; --n)
    {
      MoveNotCoprimeSP(SPL,
                       ( mySPL  =eterm_init(TListIndetsNo(theTList)) ),
                       ( myTerm =MTL[n] ) );
      if ( eterm_get_OccIndNo(mySPL)!=0 )
      {
        MTLMoveLastToNth(MTL, MTLLen, n);
        if ( IsZero(res) )   res = OneTermAndSPPoincare(P, myTerm,mySPL);
        else  res *= OneTermAndSPPoincare(P, myTerm,mySPL);
      }
      else
      {
        DegSum += eterm_degree(myTerm);
        eterm_free(mySPL);
      }
    }
    n = IntsGetLen(Indets(SPL));
    VERBOSE(90) << "SPL len " << n << std::endl;
    if ((n==eterm_degree(SPL)) && IsZero(res) )
      //      res = HPPowerList(n);
      CopyHPPower(res, n);
    else
    {
      if ( !IsZero(res) )
      {
        //if ( UPSize(res) < UPDeg(res) + eterm_degree(SPL) + DegSum )
        //  res = UnipolyChangeSize(res, UPDeg(res) + eterm_degree(SPL) + DegSum);
      }
      else
        res = one(P); //UnipolyOne(eterm_degree(SPL) + DegSum);
      for ( ; n>0 ; --n )
        P->myMulBy1MinusXExp(raw(res), SPL[(Indets(SPL))[n]]);
    }
    SetMTListLen(theTList,MTLLen);
    for ( n=MTLLen ; n>0 ; --n)   P->myMulBy1MinusXExp(raw(res), eterm_degree( MTL[n]));
    EraseAndFreeTList(theTList);
 
    return res;
  }


  RingElem TotalSplitPoincare(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList)
  {
    RingElem res(P);
    eterm  myTerm, SPL =SPList(theTList), mySPL;
    MixedTermList  MTL =MTList(theTList);
    int MTLLen =MTListLen(theTList), n =MTLLen, DegSum =0;
 
    for ( ; n>0 ; --n)
    {
      MoveNotCoprimeSP(SPL,
                       ( mySPL  =eterm_init(TListIndetsNo(theTList)) ),
                       ( myTerm =MTL[n] ) );
      if ( eterm_get_OccIndNo(mySPL)!=0 )
      {
        MTLMoveLastToNth(MTL, MTLLen, n);
        if ( IsZero(res) )   res = OneTermAndSPPoincare(P, PPM, myTerm,mySPL);
        else  res *= OneTermAndSPPoincare(P, PPM, myTerm,mySPL);
      }
      else
      {
        DegSum += eterm_degree(myTerm);
        eterm_free(mySPL);
      }
    }
    {
      if ( IsZero(res) )
        res = one(P); //UnipolyOne(eterm_degree(SPL) + DegSum);
      for (int n = IntsGetLen(Indets(SPL)) ; n>0 ; --n )
        res *= OneMinusXExp(P, wdeg(IndetPower(PPM, (Indets(SPL))[n]-1, SPL[(Indets(SPL))[n]])));
    }
    SetMTListLen(theTList,MTLLen);
    for ( n=MTLLen ; n>0 ; --n)
      res *= OneMinusXExp(P, wdeg(PPM, MTL[n]));
    EraseAndFreeTList(theTList);
 
    return res;
  }


  RingElem BigRecPoincare(const DenseUPolyRing& P, TermList theTList);
  RingElem RecPoincare(const DenseUPolyRing& P, TermList theTList);
  RingElem RecPoincareCoeffNth(const DenseUPolyRing& P, TermList theTList, int n);

  RingElem RecPoincare(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList);
  RingElem RecPoincareCoeffNth(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList, int n);

  RingElem SplitPoincare(const DenseUPolyRing& P, TermList theTList, TermList SplitterTList)
  {
    RingElem res(P);
    int SplTLLen =MTListLen(SplitterTList), n, NewMTLLen =0, DegSum =0;
    eterm   SPL =SPList(theTList);
    TermList NewTL =NewTList(MTListLen(theTList),TListIndetsNo(theTList)),
      *auxTLs =(TermList*)mymalloc((SplTLLen+1)*sizeof(TermList), "*TermList");
    MixedTermList   SMTL =MTList(SplitterTList), NewMTL =MTList(NewTL);
    VerboseLog VERBOSE("SplitPoincare1");

    VERBOSE(90) << " SplitterTList len = " << MTListLen(SplitterTList) << std::endl;
    for ( n =MTListLen(SplitterTList) ; n>0 ; --n )
    {
      auxTLs[n] = NewTList(MTListLen(theTList),TListIndetsNo(theTList));
      MoveNotCoprime(theTList, auxTLs[n], SMTL[n]);
      eterm_free(SMTL[n]);
    }
    SetMTListLen(SplitterTList, 0);
    for ( n =SplTLLen ; n>0 ; --n )
      if ((MTListLen(auxTLs[n])!=1) ||(IntsGetLen(Indets(SPList(auxTLs[n])))!=0) )
        if ( IsZero(res) )
          res = BigRecPoincare(P, auxTLs[n]);
        else
          res *= BigRecPoincare(P, auxTLs[n]);
      else
      {
        DegSum += eterm_degree((MTList(auxTLs[n]))[1]);
        MTLPutLast(NewMTL, NewMTLLen,(MTList(auxTLs[n]))[1]);
        SetMTListLen(auxTLs[n],0);
        EraseAndFreeTList(auxTLs[n]);
      }
    myfree((SplTLLen+1)*sizeof(TermList), auxTLs, "*TermList");
    //if ( UPSize(res) < UPDeg(res) + eterm_degree(SPL) + DegSum )
    //  res = UnipolyChangeSize(res, UPDeg(res) + eterm_degree(SPL) + DegSum);
    for ( n = IntsGetLen(Indets(SPL)) ; n>0 ; --n )
      P->myMulBy1MinusXExp(raw(res), SPL[(Indets(SPL))[n]]);
    EraseAndFreeTList(theTList);
    for ( n=NewMTLLen ; n>0 ; --n)
      P->myMulBy1MinusXExp(raw(res), eterm_degree((NewMTL)[n]));
    SetMTListLen(NewTL,NewMTLLen);
    EraseAndFreeTList(NewTL);
 
    return res;
  }


  RingElem SplitPoincare(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList, TermList SplitterTList)
  {
    RingElem res(P);
    int SplTLLen =MTListLen(SplitterTList), n, NewMTLLen =0, DegSum =0;
    eterm   SPL =SPList(theTList);
    TermList NewTL =NewTList(MTListLen(theTList),TListIndetsNo(theTList)),
      *auxTLs =(TermList*)mymalloc((SplTLLen+1)*sizeof(TermList), "*TermList");
    MixedTermList   SMTL =MTList(SplitterTList), NewMTL =MTList(NewTL);
    VerboseLog VERBOSE("SplitPoincare2");

    VERBOSE(90) << " SplitterTList len = " << MTListLen(SplitterTList) << std::endl;
    for ( n =MTListLen(SplitterTList) ; n>0 ; --n )
    {
      auxTLs[n] = NewTList(MTListLen(theTList),TListIndetsNo(theTList));
      MoveNotCoprime(theTList, auxTLs[n], SMTL[n]);
      eterm_free(SMTL[n]);
    }
    SetMTListLen(SplitterTList, 0);
    for ( n =SplTLLen ; n>0 ; --n )
      if ((MTListLen(auxTLs[n])!=1) ||(IntsGetLen(Indets(SPList(auxTLs[n])))!=0) )
        if ( IsZero(res) )
          res = RecPoincare(P, PPM, auxTLs[n]);
        else
          res *= RecPoincare(P, PPM, auxTLs[n]);
      else
      {
        DegSum += eterm_degree((MTList(auxTLs[n]))[1]);
        MTLPutLast(NewMTL, NewMTLLen,(MTList(auxTLs[n]))[1]);
        SetMTListLen(auxTLs[n],0);
        EraseAndFreeTList(auxTLs[n]);
      }
    myfree((SplTLLen+1)*sizeof(TermList), auxTLs, "*TermList");
    for ( n = IntsGetLen(Indets(SPL)) ; n>0 ; --n )
      res *= OneMinusXExp(P, wdeg(IndetPower(PPM, (Indets(SPL))[n]-1, SPL[(Indets(SPL))[n]])));
    EraseAndFreeTList(theTList);
    for ( n=NewMTLLen ; n>0 ; --n)
      res *= OneMinusXExp(P, wdeg(PPM, NewMTL[n]));
    SetMTListLen(NewTL,NewMTLLen);
    EraseAndFreeTList(NewTL);
 
    return res;
  }


  RingElem BigRecPoincare(const DenseUPolyRing& P, TermList theTList)
  {
    int MTLLen = MTListLen(theTList), PivotDeg, PivotIndex;
    TermList DivTList, SplitterTL;
    VerboseLog VERBOSE("BigRecPoincare_DUP");

    if (MTLLen == 0) return    SPPoincare(P, theTList);
    if (MTLLen == 1) return    LenOnePoincare(P, theTList);
    BigPivotOf(theTList, &PivotIndex, &PivotDeg);
    VERBOSE(80) << "  PivotIndex " << PivotIndex
                << "  len = " << MTLLen
                << "  SPL len = " << IntsGetLen(Indets(SPList(theTList)))
                << std::endl;
    if ( PivotIndex == 0 )  return    TotalSplitPoincare(P, theTList);
    if ((MTLLen < TListIndetsNo(theTList)) && (MTLLen > 4) )
      if ((SplitterTL=SplitIndets(theTList))!= NULL )
        return SplitPoincare(P, theTList, SplitterTL);
    ReduceAndDivideBySimplePower(theTList, &DivTList, PivotIndex, PivotDeg);
    RingElem f(BigRecPoincare(P, theTList));
    P->myAddMulLM(raw(f), raw(one(CoeffRing(P))), PivotDeg, raw(BigRecPoincare(P,DivTList)));
    return f;
  }

  RingElem RecPoincare(const DenseUPolyRing& P, TermList theTList)
  {
    int MTLLen = MTListLen(theTList), PivotDeg;
    eterm Pivot;
    TermList DivTList, SplitterTL;
    VerboseLog VERBOSE("RecPoincare_DUP");
    
    if (MTLLen == 0) return    SPPoincare(P, theTList);
    if (MTLLen == 1) return    LenOnePoincare(P, theTList);
    VERBOSE(80) << " GCD3PivotOf len =" << MTLLen << std::endl;
    if ((Pivot =GCD3PivotOf(theTList)) == NULL )
      return    TotalSplitPoincare(P, theTList);
    if ((MTLLen < TListIndetsNo(theTList)) && (MTLLen > 4) )
      if ((SplitterTL=SplitIndets(theTList))!= NULL )
      {
        eterm_free(Pivot);
        return     SplitPoincare(P, theTList, SplitterTL);
      }
    PivotDeg = eterm_degree(Pivot);
    ReduceAndDivideByPivot(theTList, &DivTList, Pivot);
  
    //   return P1PlusXExpP2(RecPoincare(P, theTList),  PivotDeg, 
    // 		      RecPoincare(P, DivTList));
    return RecPoincare(P, theTList) +
      IndetPower(P,0,PivotDeg)*RecPoincare(P, DivTList);
    /*
      return P1TimesP2(unipoly_dup(PowerList[LPSNo]),
      P1PlusXExpP2(RecPoincare(theTList), PivotDeg, 
      RecPoincare(DivTList)),
      P_R);
    */
  }

  RingElem RecPoincare(const SparsePolyRing& P, const PPMonoid& PPM, TermList theTList)
  {
    int MTLLen = MTListLen(theTList);
    eterm Pivot;
    TermList DivTList, SplitterTL;
    VerboseLog VERBOSE("RecPoincare_SP");
    
    if (MTLLen == 0) return    SPPoincare(P, PPM, theTList);
    if (MTLLen == 1) return    LenOnePoincare(P, PPM, theTList);
    VERBOSE(80) << " GCD3PivotOf len =" << MTLLen << std::endl;
    if ((Pivot =GCD3PivotOf(theTList)) == NULL )
      return    TotalSplitPoincare(P, PPM, theTList);
    if ((MTLLen < TListIndetsNo(theTList)) && (MTLLen > 4) )
      if ((SplitterTL=SplitIndets(theTList))!= NULL )
      {
        eterm_free(Pivot);
        return     SplitPoincare(P, PPM, theTList, SplitterTL);
      }
    degree PivotDeg = wdeg(PPM, Pivot);
    ReduceAndDivideByPivot(theTList, &DivTList, Pivot);
  
    //   return P1PlusXExpP2(RecPoincare(P, theTList),  PivotDeg, 
    // 		      RecPoincare(P, DivTList));
    return RecPoincare(P, PPM, theTList) +
      XExp(P,PivotDeg)*RecPoincare(P, PPM, DivTList);
    /*
      return P1TimesP2(unipoly_dup(PowerList[LPSNo]),
      P1PlusXExpP2(RecPoincare(theTList), PivotDeg, 
      RecPoincare(DivTList)),
      P_R);
    */
  }

  RingElem TLPoincareNumeratorCPP(const DenseUPolyRing& HSRing, TermList theTList)
  {
    return BigRecPoincare(HSRing, theTList);
  }

  RingElem TLPoincareNumeratorCPP(const SparsePolyRing& HSRing, const PPMonoid& PPM, TermList theTList)
  {
    return RecPoincare(HSRing, PPM, theTList);
  }

}
