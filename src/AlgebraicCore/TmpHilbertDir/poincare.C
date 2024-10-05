//   Copyright (c)  2006  Anna Bigatti,  John Abbott

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

#include "IVectors.h"
#include "eterms.h"
#include "poincare.h"
#include "unipoly.h"

/**********************************************************/

void poincare_init(int IndetsNo)
{
  static int initialized = 0;
  if ( !initialized )
  {
    initialized = 1;

    rum_init_all();
    //------- disable RUM ----------------------------------
    rum_init(12,2000); /* pointers */
    rum_init(200*sizeof(eterm), 120);
    rum_init(sizeof(TermListStruct), RUM_STD_SIZE);

    unipoly_rum_init(PoincareMaxPower);

    rum_init(eterm_size(IndetsNo), RUM_STD_SIZE);
    //------- disable RUM ----------------------------------

    PowerList = MakePowerList(PoincareMaxPower);
    GlobalTermList_init(IndetsNo);
  }
}


void StartPoincare(int IndetsNo)
{
  poincare_init(IndetsNo);
  if ( GlobalTermList_size()!=IndetsNo )
  {
    GlobalTermList_free();
    GlobalTermList_init(IndetsNo);
  }
}


void EndPoincare(int /*IndetsNo*/)
{
  FreeUnipolyList(PowerList, PoincareMaxPower);
  GlobalTermList_free();
}

/****************************************************************
 *************************  BASE CASES  *************************
 ****************************************************************/

unipoly UnipolySPPoincare(TermList theTList)
{
  eterm   SPL = SPList(theTList);
  ints OccInd = Indets(SPL);
  int i = IntsGetLen(OccInd);
  unipoly res;

  if ( i==eterm_degree(SPL) )    res = PowerExtDup(i, i);
  else
  {
    res = UnipolyOne(eterm_degree(SPL));
    for ( ; i>0 ; --i )   MultByOneMinusXExp(res, SPL[OccInd[i]]);
  }
  EraseAndFreeTList(theTList);
  return res;
}

unipoly UnipolyLenOnePoincare(TermList theTList)
{
  eterm    SPL = SPList(theTList), myTerm = (MTList(theTList))[1];
  ints     OccInd = Indets(SPL);
  int   TDeg =eterm_degree(myTerm), i =IntsGetLen(OccInd), exp;
  unipoly p1, p2;

  if ( eterm_coprime(myTerm, SPL) )
  {
    if ( i==eterm_degree(SPL) )   p1 = PowerExtDup(i, i+TDeg);
    else
    {
      p1 = UnipolyOne(eterm_degree(SPL)+TDeg);
      for ( ; i>0 ; --i )  MultByOneMinusXExp(p1, SPL[OccInd[i]]);
    }
    MultByOneMinusXExp(p1, TDeg);
    EraseAndFreeTList(theTList);
    return p1;
  }
  else 
  {
    p1 = UnipolyOne(eterm_degree(SPL)+TDeg);
    p2 = UnipolyOne(eterm_degree(SPL));
    for ( ; i>0 ; --i )
    {
      exp = SPL[OccInd[i]];
      MultByOneMinusXExp(p1, exp);
      MultByOneMinusXExp(p2, exp-myTerm[OccInd[i]]);
    }   
    EraseAndFreeTList(theTList);
    return P1MinusXExpP2(p1, TDeg, p2);
  }
}

unipoly UnipolyOneTermAndSPPoincare(eterm theTerm, eterm SPL)
{
  ints  OccInd = Indets(SPL);
  int   TDeg = eterm_degree(theTerm), i = IntsGetLen(OccInd), exp;
  unipoly p1, p2 ;

  if (eterm_coprime(theTerm, SPL))
  {
    if ( i==eterm_degree(SPL) )   p1 = PowerExtDup(i, i+TDeg);
    else
    {
      p1 = UnipolyOne(eterm_degree(SPL)+TDeg);
      for ( ; i>0 ; --i )  MultByOneMinusXExp(p1, SPL[OccInd[i]]);
    }
    MultByOneMinusXExp(p1, eterm_degree(theTerm));
    eterm_free(theTerm);
    eterm_free(SPL);
    return p1;
  }

  p1 = UnipolyOne(eterm_degree(SPL)+TDeg);
  p2 = UnipolyOne(eterm_degree(SPL));
  for ( ; i>0 ; --i )
  {
    exp = SPL[OccInd[i]];
    MultByOneMinusXExp(p1, exp);
    MultByOneMinusXExp(p2, exp-theTerm[OccInd[i]]);
  }
  eterm_free(theTerm);
  eterm_free(SPL);
  return P1MinusXExpP2(p1, TDeg, p2);
}

/* TermLists.c */
/* TermLists.c: END */

/*************************************************************
 ************************    SPLIT    ************************
 *************************************************************/
 
unipoly TotalSplitPoincare(TermList theTList)
{
  unipoly   res = NULL;
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
      if ( res == NULL )   res = UnipolyOneTermAndSPPoincare(myTerm,mySPL);
      else  res = P1TimesP2(res, UnipolyOneTermAndSPPoincare(myTerm,mySPL));
    }
    else
    {
      DegSum += eterm_degree(myTerm);
      eterm_free(mySPL);
    }
  }
  n = IntsGetLen(Indets(SPL));
  if ((n==eterm_degree(SPL)) && res == NULL )
    res = PowerExtDup(n, n + DegSum);
  else
  {
    if ( res != NULL )
    {
      if ( UPSize(res) < UPDeg(res) + eterm_degree(SPL) + DegSum )
        res = UnipolyChangeSize(res, UPDeg(res) + eterm_degree(SPL) + DegSum);
    }
    else
      res = UnipolyOne(eterm_degree(SPL) + DegSum);
    for ( ; n>0 ; --n )  MultByOneMinusXExp(res, SPL[(Indets(SPL))[n]]);
  }
  SetMTListLen(theTList,MTLLen);
  for ( n=MTLLen ; n>0 ; --n)   MultByOneMinusXExp(res, eterm_degree( MTL[n]));
  EraseAndFreeTList(theTList);
 
  return res;
}


unipoly BigRecPoincare(TermList theTList);
unipoly RecPoincare(TermList theTList);
unipoly RecPoincareCoeffNth(TermList theTList, int n);

unipoly SplitPoincare(TermList theTList, TermList SplitterTList)
{
  unipoly   res;
  int SplTLLen =MTListLen(SplitterTList), n, NewMTLLen =0, DegSum =0;
  eterm   SPL =SPList(theTList);
  TermList NewTL =NewTList(MTListLen(theTList),TListIndetsNo(theTList)),
           *auxTLs =(TermList*)mymalloc((SplTLLen+1)*sizeof(TermList), "*TermList");
  MixedTermList   SMTL =MTList(SplitterTList), NewMTL =MTList(NewTL);
  
  res = NULL;
  for ( n =MTListLen(SplitterTList) ; n>0 ; --n )
  {
    auxTLs[n] = NewTList(MTListLen(theTList),TListIndetsNo(theTList));
    MoveNotCoprime(theTList, auxTLs[n], SMTL[n]);
    eterm_free(SMTL[n]);
  }
  SetMTListLen(SplitterTList, 0);
  for ( n =SplTLLen ; n>0 ; --n )
    if ((MTListLen(auxTLs[n])!=1) ||(IntsGetLen(Indets(SPList(auxTLs[n])))!=0) )
      if ( res == NULL )
        res = BigRecPoincare(auxTLs[n]);
      else
        res = P1TimesP2(res, BigRecPoincare(auxTLs[n]));
    else
    {
      DegSum += eterm_degree((MTList(auxTLs[n]))[1]);
      MTLPutLast(NewMTL, NewMTLLen,(MTList(auxTLs[n]))[1]);
      SetMTListLen(auxTLs[n],0);
      EraseAndFreeTList(auxTLs[n]);
    }
  myfree((SplTLLen+1)*sizeof(TermList), auxTLs, "*TermList");
  if ( UPSize(res) < UPDeg(res) + eterm_degree(SPL) + DegSum )
    res = UnipolyChangeSize(res, UPDeg(res) + eterm_degree(SPL) + DegSum);
  for ( n = IntsGetLen(Indets(SPL)) ; n>0 ; --n )
    MultByOneMinusXExp(res, SPL[(Indets(SPL))[n]]);
  EraseAndFreeTList(theTList);
  for ( n=NewMTLLen ; n>0 ; --n)
    MultByOneMinusXExp(res, eterm_degree((NewMTL)[n]));
  SetMTListLen(NewTL,NewMTLLen);
  EraseAndFreeTList(NewTL);
 
  return res;
}


unipoly BigRecPoincare(TermList theTList)
{
  int MTLLen = MTListLen(theTList), PivotDeg, PivotIndex;
  TermList DivTList, SplitterTL;

  if (MTLLen == 0) return    UnipolySPPoincare(theTList);
  if (MTLLen == 1) return    UnipolyLenOnePoincare(theTList);
  BigPivotOf(theTList, &PivotIndex, &PivotDeg);
  if ( PivotIndex == 0 )  return    TotalSplitPoincare(theTList);
  if ((MTLLen < TListIndetsNo(theTList)) && (MTLLen > 4) )
    if ((SplitterTL=SplitIndets(theTList))!= NULL )
      return SplitPoincare(theTList, SplitterTL);
  ReduceAndDivideBySimplePower(theTList, &DivTList, PivotIndex, PivotDeg);
  
  return P1PlusXExpP2(BigRecPoincare(theTList),  PivotDeg, 
		      BigRecPoincare(DivTList));
}

unipoly RecPoincare(TermList theTList)
{
  int MTLLen = MTListLen(theTList), PivotDeg;
  eterm Pivot;
  TermList DivTList, SplitterTL;
    
  if (MTLLen == 0) return    UnipolySPPoincare(theTList);
  if (MTLLen == 1) return    UnipolyLenOnePoincare(theTList);
  if ((Pivot =GCD3PivotOf(theTList)) == NULL )  
    return    TotalSplitPoincare(theTList);
  if ((MTLLen < TListIndetsNo(theTList)) && (MTLLen > 4) )
    if ((SplitterTL=SplitIndets(theTList))!= NULL )
    {
      eterm_free(Pivot);
      return     SplitPoincare(theTList, SplitterTL);
    }
  PivotDeg = eterm_degree(Pivot);
  ReduceAndDivideByPivot(theTList, &DivTList, Pivot);
  
  return P1PlusXExpP2(RecPoincare(theTList),  PivotDeg, 
		      RecPoincare(DivTList));
/*
  return P1TimesP2(unipoly_dup(PowerList[LPSNo]),
		   P1PlusXExpP2(RecPoincare(theTList), PivotDeg, 
				RecPoincare(DivTList)),
		   P_R);
*/
}

unipoly TLPoincareNumerator(TermList theTList)
{
  //  StartPoincare(TListIndetsNo(theTList));
  
  return BigRecPoincare(theTList);
}
