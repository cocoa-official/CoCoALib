//   Copyright (c)  2006-2007  Anna Bigatti,  John Abbott

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
#include "TermList.h"


ints MFI_Occurrences, MFI_Indets;
int Ints120[121];
MixedTermList  MTL120[121];
TermList GlobalSplitterTList;

/**********************************************************/

#define malloc_MTList(len) \
  (MixedTermList)mymalloc( (len+1)*sizeof(eterm), "MTList" )
#define free_MTList(len,p)     myfree((len+1)*sizeof(eterm),p,"MTList")
#define realloc_MTList(p,oldsz,newsz) \
  (MixedTermList)myrealloc(p,(oldsz+1)*sizeof(eterm),\
  (newsz+1)*sizeof(eterm),"MTList" )

/********************    list of terms    ********************/
void GlobalTermList_init(int IndetsNo)
{
  //------- disable RUM ----------------------------------
  rum_init(eterm_size(IndetsNo), RUM_STD_SIZE);
  //------- disable RUM ----------------------------------
  MFI_Occurrences = ints_init(IndetsNo);
  MFI_Indets = ints_init(IndetsNo);
  GlobalSplitterTList = NewTList(IndetsNo, IndetsNo);
}

void GlobalTermList_free()
{
  ints_free(MFI_Occurrences);
  ints_free(MFI_Indets);
  EraseAndFreeTList(GlobalSplitterTList);
}

int GlobalTermList_size()
{
  return IntsGetLen(MFI_Occurrences);
}

/*************************************************************/


#define malloc_TList() (TermList)mymalloc(sizeof(TermListStruct),"TList")
#define free_TList(p) {myfree(sizeof(TermListStruct),p,"TList");}


/* TermLists.c */

#define Random(len)( 1 )

TermList NewTList(int len, int indetsNo)
{
  TermList  aux = malloc_TList();

  if ( len <= 199 ) len = 199;
  SetMTListSize (aux, len);
  SetMTListLen  (aux, 0);
  MTList(aux) = malloc_MTList(len);
  SPList(aux) = eterm_init(indetsNo);
  
  return aux;
}


void InsInTList(TermList TL, eterm t, int *NewLen)
{
  ints OccInd = Indets(t);
  if ( IntsGetLen(OccInd) == 1 )
  {
    InsInSPList (OccInd[1], eterm_get_nth(t,OccInd[1]), SPList(TL) );
    eterm_free(t);   
  }
  else
    MTLPutLast (MTList(TL), *NewLen, t);
}


void EraseAndFreeTList(TermList theTList)
/* all the terms must be still allocated!! */
{
  int i =MTListLen(theTList);
  MixedTermList MTL =MTList(theTList);
  
  for ( ; i>=0 ; --i )  eterm_free(MTL[i]);
  free_MTList(MTListSize(theTList), MTL);
  free_TList(theTList);
}

TermList TListMake(eterm theSPL, MixedTermList theMTL, int Size, int Len)
{
  TermList aux;

  aux = malloc_TList();
  MTList(aux) = theMTL;
  SPList(aux) = theSPL;
  SetMTListSize(aux, Size);
  TListReduceSize(aux, Len);

  return aux;
}

void TListChangeSize (TermList theTList, int NewSize)
{
  MTList(theTList) = realloc_MTList(MTList(theTList),
				    MTListLen(theTList), NewSize);
  SetMTListSize(theTList, NewSize);
}

void MTLAppend(MixedTermList MTL1, int *MTL1Len, 
	       MixedTermList MTL2, int MTL2Len)
/* allocated memory for MTL1 enough to contain MTL1+MTL2 */
{
  int i;

  for ( i=MTL2Len ; i>0 ; --i )    MTLPutLast(MTL1, *MTL1Len, MTL2[i]);
}

void TListAppendMTL(TermList TL, MixedTermList MTL2, int MTL2Len)
/* allocated memory for MTL1 enough to contain MTL1+MTL2 */
{
  int i, MTL1Len =MTListLen(TL);
  MixedTermList MTL1=MTList(TL);
  
  for ( i=MTL2Len ; i>0 ; --i )  MTLPutLast(MTL1, MTL1Len, MTL2[i]);
  SetMTListLen(TL, MTL1Len);
}

/********    Simple Power List    ********/

void InsInSPList(size_t Index, unsigned Exp, eterm SPL)
{
  size_t CurrentExp =SPL[Index];
  
  if (CurrentExp == 0)
  {  
    eterm_put_non0_nth(SPL, Index, Exp);
    IntsPutLast(Indets(SPL), Index);
  }
  else if (CurrentExp > Exp)
    eterm_put_non0_nth(SPL, Index, Exp);
}

#define SPLDividesTerm(SPL,t)  sp_SmallMult(t,SPL) 

coc_bool TListSimplePowersDivide(TermList theTList, eterm t)
{
  eterm SPL =SPList(theTList);
  int   i;
  
  for ( i=eterm_get_indetsNo(SPL) ; i>0 ; --i )
    if ( SPL[i]!=0 && SPL[i]<=t[i] )   return TRUE;
  return FALSE;
}

/************************    SPLIT    ************************/
 
void DelLinearSimplePowers(TermList theTList, int *LPSNo)
{
  eterm SPL = SPList(theTList);
  ints  OccInd = Indets(SPL);
  size_t   n = IntsGetLen(OccInd), sum = 0, index;
  
  for ( ; n>0 ; --n )
    if ( SPL[(index=OccInd[n])] == 1 )
    {
      ++sum;
      eterm_put0_nth(SPL, index);
      IntsMoveLastToNth(OccInd, n); 
    }
  *LPSNo = sum;  
}

eterm GetCoprimeSPList(TermList theTList)
{
  int    MTLLen =MTListLen(theTList), *OccInd;
  MixedTermList   MTL =MTList(theTList);
  size_t    i, n =MTLLen, index;
  eterm res, theSPL; 
  shortbitset CoprimeSupp;

  if ( IntsGetLen(Indets(SPList(theTList))) == 0 )   return NULL;
  CoprimeSupp = SqFr(SPList(theTList));  
  for ( ; n > 0 ; --n)
  {
    CoprimeSupp &=(~(SqFr(MTL[n]) ) );
    if ( CoprimeSupp == 0 )   return NULL;
  }
  theSPL = SPList(theTList);
  OccInd = Indets(theSPL);  
  res = eterm_init(eterm_get_indetsNo(theSPL));  
  for ( i=IntsGetLen(OccInd) ; i>0 ; --i )
  {
    if ( bs_get_n(CoprimeSupp, (size_t)OccInd[i]) )
    {
      index = OccInd[i];
      eterm_put_non0_nth(res, index, theSPL[index]);
      IntsPutLast(Indets(res), index);
      eterm_put0_nth(theSPL, index);      
      IntsMoveLastToNth(OccInd, i);
    }
  }  
  return res;
}


void MoveNotCoprimeSP(eterm FromSPList, eterm ToSPList, eterm theTerm)
{
  ints OccInd = Indets(FromSPList);
  size_t   i = IntsGetLen(OccInd), index;
  
  for ( ; i>0 ; --i )
    if ( theTerm[(index = OccInd[i])] != 0 )
    {
      eterm_put_non0_nth(ToSPList, index, FromSPList[index]);
      IntsPutLast(Indets(ToSPList), index);
      eterm_put0_nth(FromSPList, index);
      IntsMoveLastToNth(OccInd, i);
    }
}

void MoveNotCoprime(TermList FromTList, TermList ToTList, eterm theTerm)
{
  int       FromMTLLen =MTListLen(FromTList),  ToMTLLen=0;
  MixedTermList   FromMTL =MTList(FromTList),  ToMTL =MTList(ToTList);
  int     n =FromMTLLen;
  
  for ( ; n > 0 ; --n)
    if ( !eterm_coprime(FromMTL[n], theTerm) ) 
    { 
      MTLPutLast(ToMTL, ToMTLLen, FromMTL[n]);
      MTLMoveLastToNth(FromMTL, FromMTLLen, n);
    }
  TListReduceSize(FromTList, FromMTLLen);
  TListReduceSize(ToTList, ToMTLLen);
  
  MoveNotCoprimeSP(SPList(FromTList), SPList(ToTList), theTerm);
}

TermList SplitIndets(TermList theTList)
{
  TermList      res =GlobalSplitterTList;
  int           resLen = 0;
  MixedTermList MTL = MTList(theTList),  resMTL =MTList(res);
  int   n = MTListLen(theTList),  m;
  eterm tn, tm;

  for ( ; n > 0 ; --n)
  {
    tn = MTL[n];
    for ( m=resLen ; m>0 ; --m)
    {
      tm = resMTL[m];
      if ( !eterm_coprime(tn, tm)  )
      {
        eterm_union_and_assign(tm, tn);
        if ( tn != MTL[n] ) eterm_free(tn);
        tn = tm;
        MTLMoveLastToNth(resMTL, resLen, m);
      }  
    }
    if ( tn != MTL[n] )    MTLPutLast(resMTL, resLen, tn);
    else MTLPutLast(resMTL, resLen, eterm_dup(tn));
  }
  if ( resLen == 1 )
  {
    eterm_free(resMTL[1]);
    SetMTListLen(res, 0);
    return NULL;
  }
  SetMTListLen(res, resLen);
  
  return res;
}

/************************    INTERREDUCE    ************************/

void MTLOrderByDecrDegree(MixedTermList theMTList, int Len, int MaxDeg)
{
  eterm t;
  MixedTermList *MTLDeg=(MixedTermList*)mycalloc
   ((MaxDeg+1),sizeof(MixedTermList),"*MixedTermList");
  int *MTLLenDeg=(int*)mycalloc((MaxDeg+1), sizeof(int),"*int"), auxLen =0;
  int i, d;

  for ( d=MaxDeg ; d>=2 ; --d)    MTLDeg[d] = malloc_MTList(Len);

  for ( i=Len ; i>0 ; --i)
  {
    d = eterm_degree((t = theMTList[i]) ); 
    MTLPutLast( MTLDeg[d], MTLLenDeg[d], t);
  }
  for ( d=MaxDeg ; d>1 ; --d)
  {
    if (MTLLenDeg[d] != 0)
      for ( i=MTLLenDeg[d]; i>0; --i)
        MTLPutLast(theMTList, auxLen, (MTLDeg[d])[i]);
    free_MTList(Len, MTLDeg[d]);
  }
  myfree((MaxDeg+1)*sizeof(MixedTermList), MTLDeg, "*MixedTermList");  
  myfree((MaxDeg+1)*sizeof(int), MTLLenDeg, "*int");

}

coc_bool MTLDividesTerm(MixedTermList MTL, int MTLLen, eterm T)
{
  int i = MTLLen;

  for ( ; i>0 ; --i ) if ( eterm_divides(MTL[i], T) )  return TRUE;

  return FALSE;
}

void InterreduceTList(TermList theTList)
{
  MixedTermList   MTL =MTList(theTList);
  int      MTLLen =MTListLen(theTList);
  int      n =MTLLen, j, MaxDeg =0;
  eterm    auxT;

  for ( ; n> 0 ; --n)
    if ( TListSimplePowersDivide(theTList,(auxT=MTL[n]) )  )
    {
      eterm_free(auxT);
      MTLMoveLastToNth(MTL, MTLLen, n);
    }
    else
      MaxDeg = MAX(MaxDeg, eterm_degree(auxT));
  if (MTLLen > 1)
  {
    MTLOrderByDecrDegree(MTL, MTLLen, MaxDeg);
    for ( n=MTLLen-1 ; n>0; --n)
    {
      auxT =  MTL[n];
      for ( j = MTLLen ; j>n ; --j)
        if ( eterm_divides(MTL[j], auxT)  )
        {
          eterm_free(auxT);
          MTLMoveLastToNth(MTL, MTLLen, n);
          break;
        }
    }
  }
  SetMTListLen(theTList, MTLLen);
}

/************************  PIVOT  ************************/

int MostFrequentIndet(TermList theTList)
{
  MixedTermList MTL =MTList(theTList);
  int IndNo = TListIndetsNo(theTList), exp, MFIndNo = 1, i, j;
  ints OccInd;
  
  for ( j=IndNo ; j>0 ; --j )    MFI_Occurrences[j] = 0;
  for ( i=MTListLen(theTList) ; i>0 ; --i)
  {
    OccInd = Indets(MTL[i]);
    for ( j=IntsGetLen(OccInd) ; j>0 ; --j )   MFI_Occurrences[OccInd[j]]++;
  }
  MFI_Indets[1] = (j=IndNo);
  i = MFI_Occurrences[j];
  for ( --j ; j>0 ; --j )
    if ((exp=MFI_Occurrences[j]) >= i)
    {
      if (exp > i)  { MFIndNo = 0;  i = exp; }
      MFI_Indets[++MFIndNo] = j;
    }  
 
  if ( i == 1 )       return 0;
  if ( MFIndNo == 1 ) return MFI_Indets[1];
  else                return MFI_Indets[MFIndNo/2];
}

MixedTermList GetTermsWithIndet(TermList theTList, int Indet, int *resLen)
{
  int i =MTListLen(theTList);
  eterm t;
  MixedTermList res =malloc_MTList(i), MTL =MTList(theTList);
  
  *resLen = 0;
  for ( ; i>0 ; --i)
    if ( (t =MTL[i])[Indet] )   MTLPutLast(res, *resLen, t);
  return res;
}

eterm VarPivotOf(TermList theTList)
{
  eterm res;
  size_t MFIndet;

  if ((MFIndet=MostFrequentIndet(theTList)) == 0)   return NULL;
 
  res = eterm_init(TListIndetsNo(theTList));
  eterm_put_non0_nth(res, MFIndet, 1);
  IntsPutLast(Indets(res), MFIndet);

  return res;
}

eterm GCD3PivotOf(TermList theTList)
{
  eterm res, t;
  int MTLLen, MFIndet, r;
  MixedTermList MTL;
  eterm t1, t2;

  if ((MFIndet=MostFrequentIndet(theTList)) == 0)   return NULL;
 
  MTL = GetTermsWithIndet(theTList, MFIndet, &MTLLen);

  t1 =  MTL[(r =Random(MTLLen)) ];
  MTLMoveLastToNth(MTL, MTLLen, r);
  t2 =  MTL[(r =Random(MTLLen)) ];
  MTLMoveLastToNth(MTL, MTLLen, r);
  res = eterm_gcd(t1, t2);
  if (MTLLen != 0)
  {
    res = eterm_gcd((t=res),  MTL[Random(MTLLen)]);
    eterm_free(t);
  }
  free_MTList(MTListLen(theTList),MTL);

  return res;
}

void BigPivotOf(TermList theTList, int *PivotIndex, int *PivotDeg)
{  
  int MTLLen, MFIndet, r, e1, e2;  
  MixedTermList MTL;

  MFIndet = MostFrequentIndet(theTList);
  *PivotIndex = MFIndet;
  if ( MFIndet == 0 ) return;
  
  MTL = MTList(theTList);
  MTLLen = MTListLen(theTList);
  r = 1; /* fix macro */
  while((e1 =(MTL[++r])[MFIndet]) == 0 ) {}
  r = MTLLen;
  while((e2 =(MTL[r--])[MFIndet]) == 0 ) {}
 
  *PivotDeg   = MIN(e1,e2);
}

/*
#define PivotOf(theTList,MTLLen)\
    ((MTLLen<3)?GCD3PivotOf(theTList):BigPivotOf(theTList))
*/

#define PivotOf(theTList,MTLLen) BigPivotOf(theTList)

eterm eterm_colon_SP(eterm theTerm, size_t index, unsigned exp)
{
  unsigned OldExp;
 
  if ((OldExp =theTerm[index]) != 0 )
  {
    if ( OldExp <= exp ) printf("OldExp <= exp");
    eterm_put_non0_nth(theTerm, index, OldExp-exp );
  }
  return theTerm;
}

void OccIndDel(eterm T, int index)
{
  ints    OccInd = Indets(T);
  int     n = IntsGetLen(OccInd);

  while ( OccInd[n] != index )  --n;
  IntsMoveLastToNth(OccInd, n);
}

void ReduceAndDivideBySimplePower( TermList theTList, TermList *DivTList,
                                   size_t PIndex, unsigned PExp)
{
  int MTLLen = MTListLen(theTList), auxLen = MAX(MTListLen(theTList),199);
  int DivMTLLen, MTLLenEe, i;
  unsigned TExp, e, index;
  ints DivMTLLenExp;
  eterm  T;
  MixedTermList   *DivMTLExp, MTL =MTList(theTList), DivMTL, MTLEe;
  eterm   DivSPL;
  
  if ( PExp>119 )
  {    
    DivMTLLenExp = (int*)mymalloc((PExp+2)*sizeof(int), "*int");
    DivMTLExp = (MixedTermList*)mymalloc((PExp+2)*sizeof(MixedTermList),
                                         "*MixedTermList");
  }
  else
  {  DivMTLLenExp = Ints120;     DivMTLExp = MTL120;   }
  for ( i=PExp+1 ; i>=0 ; --i )
  {  DivMTLExp[i] = malloc_MTList(auxLen);    DivMTLLenExp[i] = 0;  }
  
  /* DivTList */
  DivSPL = eterm_colon_SP(eterm_dup(SPList(theTList)), PIndex, PExp);
  /* theTList */
  InsInSPList(PIndex, PExp, SPList(theTList));
 
  for ( i=MTLLen ; i > 0 ; --i)
  {
    TExp = (T=MTL[i])[PIndex];
    if ( TExp > PExp )
    {
      MTLMoveLastToNth(MTL, MTLLen, i);
      eterm_put_non0_nth(T, PIndex, TExp-PExp);
      MTLPutLast(DivMTLExp[PExp+1], DivMTLLenExp[PExp+1], T);
    }
    else if ( TExp == PExp )
    {
      MTLMoveLastToNth(MTL, MTLLen, i);
      if ( eterm_get_OccIndNo(T) == 2 )
      {
	if ((index = (Indets(T))[1]) == PIndex )  index = (Indets(T))[2];
        InsInSPList(index, T[index], DivSPL);
	eterm_free(T);
      }
      else
      {
	eterm_put0_nth(T, PIndex);
	OccIndDel(T, PIndex);
        MTLPutLast(DivMTLExp[TExp], DivMTLLenExp[TExp], T);
      }
    }
    else if ( TExp != 0 )
    {
      if ( eterm_get_OccIndNo(T) == 2 )
      {
	if ((index = (Indets(T))[1]) == PIndex )  index = (Indets(T))[2];
        InsInSPList(index, T[index], DivSPL);
      }
      else
      {
	eterm_put0_nth(T, PIndex);
	OccIndDel(T, PIndex);
        MTLPutLast(DivMTLExp[TExp], DivMTLLenExp[TExp], T);
      }
    }
    else /* TExp == 0 */ MTLPutLast(DivMTLExp[TExp], DivMTLLenExp[TExp], T);
  }
 
  /* DivTList(Interreduction) */
  DivMTL    = DivMTLExp[PExp];
  DivMTLLen = DivMTLLenExp[PExp];
  for ( e=PExp-1 ; e>0 ; --e )
  {
    MTLEe    = DivMTLExp[e];
    MTLLenEe = DivMTLLenExp[e];
    for ( i=MTLLenEe ; i>0 ; --i )
    {
      T = MTLEe[i];
      if ( SPLDividesTerm(DivSPL, T) )  MTLMoveLastToNth(MTLEe, MTLLenEe, i);
      else if ( MTLDividesTerm(DivMTL, DivMTLLen, T) )
        MTLMoveLastToNth(MTLEe, MTLLenEe, i);
      else
        MTLPutNth(MTLEe, i, eterm_dup(T));
      eterm_put_non0_nth(T, PIndex, e);
      IntsPutLast(Indets(T), PIndex);
    }
    if ( MTLLenEe!=0 )    MTLAppend(DivMTL, &DivMTLLen, MTLEe, MTLLenEe);
    free_MTList(auxLen, MTLEe);
  }
  MTLEe    = DivMTLExp[0];
  MTLLenEe = DivMTLLenExp[0];
  for ( i=MTLLenEe ; i>0 ; --i )
  {
    T = MTLEe[i];
    if ( SPLDividesTerm(DivSPL, T) )   MTLMoveLastToNth(MTLEe, MTLLenEe, i);
    else if ( MTLDividesTerm(DivMTL, DivMTLLen, T) )
      MTLMoveLastToNth(MTLEe, MTLLenEe, i);
    else
      MTLPutNth(MTLEe, i, eterm_dup(T));
  }
  if ( MTLLenEe!=0 )    MTLAppend(DivMTL, &DivMTLLen, MTLEe, MTLLenEe);
  free_MTList(auxLen, MTLEe);
  
  MTLAppend(DivMTL, &DivMTLLen, DivMTLExp[PExp+1], DivMTLLenExp[PExp+1]);
  free_MTList(auxLen, DivMTLExp[PExp+1]);
  *DivTList = TListMake(DivSPL, DivMTL, auxLen, DivMTLLen);
 
  /* theTList */
  TListReduceSize(theTList, MTLLen);
  if ( PExp>119 )
  {
    myfree((PExp+2)*sizeof(int), DivMTLLenExp, "*int");
    myfree((PExp+2)*sizeof(MixedTermList), DivMTLExp, "*MixedTermList");
  }
}

void ReduceAndDivideByMixedTerm( TermList theTList, TermList *DivTList,
                                 eterm Pivot)
{
  int  MTLLen =MTListLen(theTList), OldMTLLen =MTLLen,
       BMLen =0, CMTLLen =0, DivMTLLen =0;
  int index, TDeg;
  eterm     DivT,  t, auxSPL;
  int     i = MTLLen;
  MixedTermList   MTL =MTList(theTList),  BigMultMTL, CoprimeMTL, DivMTL;
  eterm   DivSPL;
 
  *DivTList = NewTList(MTLLen, TListIndetsNo(theTList));
  DivMTL = MTList(*DivTList);
  DivSPL = eterm_colon(eterm_dup(SPList(theTList)), Pivot);
  auxSPL =(SPList(*DivTList));  /* eterm_init() */

  BigMultMTL = malloc_MTList(MTLLen);
  CoprimeMTL = malloc_MTList(MTLLen);
 
  for ( ; i > 0 ; --i )
  {
    t=MTL[i];
    TDeg = eterm_degree(t);

    if ( eterm_divides(Pivot,t) )
    {
      /* theTList */
      MTLMoveLastToNth(MTL, MTLLen, i);
      /* DivTList */
      if (sp_BigMult(t,Pivot))
	MTLPutLast(BigMultMTL, BMLen, eterm_colon(t,Pivot));
      else
      {
	DivT = eterm_colon(t, Pivot);
	if ( IntsGetLen(Indets(DivT)) == 1)
	{
	  index =  (Indets(DivT))[ 1];
	  InsInSPList(index, DivT[index], DivSPL);
	  InsInSPList(index, DivT[index] + Pivot[index], auxSPL);
	  eterm_free(DivT);
	}
	else
	  MTLPutLast(DivMTL, DivMTLLen, DivT);
      }
    }
    else if ( SPLDividesTerm(auxSPL, t) )
      ;
    else if (TDeg==eterm_degree((DivT=eterm_colon(eterm_dup(t),Pivot) )) )
      MTLPutLast(CoprimeMTL, CMTLLen, DivT);
    else if ( IntsGetLen(Indets(DivT)) == 1 )
    {
      index =  (Indets(DivT))[ 1];
      InsInSPList(index, DivT[index], DivSPL);
      InsInSPList(index, DivT[index] + Pivot[index], auxSPL);
      eterm_free(DivT);
    }
    else
      MTLPutLast(DivMTL, DivMTLLen, DivT);
  }
  /* theTList */
  MTLPutLast(MTL, MTLLen, Pivot);
  TListReduceSize(theTList, MTLLen);
  /* DivTList */
  eterm_colon(auxSPL, Pivot);
  SetMTListLen(*DivTList, DivMTLLen);
  InterreduceTList(*DivTList);
  DivMTLLen = MTListLen(*DivTList);

  if ( DivMTLLen != 0 )

    for ( i=CMTLLen ; i>0 ; --i )
    {
      if ( SPLDividesTerm(DivSPL,(t=(CoprimeMTL)[i]) ) )
      {
	eterm_free(t);
	MTLMoveLastToNth(CoprimeMTL, CMTLLen, i);
      }
      else if ( MTLDividesTerm(DivMTL, DivMTLLen, t) )
      {
	eterm_free(t);
	MTLMoveLastToNth(CoprimeMTL, CMTLLen, i);
      }
    }
  else
    for ( i=CMTLLen ; i>0 ; --i )
      if ( SPLDividesTerm(DivSPL,(t=(CoprimeMTL)[i]) ) )
      {
	eterm_free(t);
	MTLMoveLastToNth(CoprimeMTL, CMTLLen, i);
      }
  SPList(*DivTList) = DivSPL;
  if (CMTLLen!=0)  MTLAppend(DivMTL, &DivMTLLen, CoprimeMTL, CMTLLen);
  if (BMLen!=0)    MTLAppend(DivMTL, &DivMTLLen, BigMultMTL, BMLen);
  free_MTList(OldMTLLen, CoprimeMTL);
  free_MTList(OldMTLLen, BigMultMTL);
  eterm_free(auxSPL);
  TListReduceSize(*DivTList, DivMTLLen);
}

void ReduceAndDivideByPivot(TermList theTList, TermList *DivTList, eterm Pivot)
{
  if ( eterm_get_OccIndNo(Pivot) == 1)
  {
    int index = (Indets(Pivot))[1], exp = eterm_degree(Pivot);
    eterm_free(Pivot);
    ReduceAndDivideBySimplePower(theTList, DivTList, index, exp);
  }
  else
    ReduceAndDivideByMixedTerm(theTList, DivTList, Pivot);
}

/* TermLists.c: END */

