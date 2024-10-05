//   Copyright (c)  2011  Anna Bigatti,  John Abbott

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
/* $Id: toric.c, 2000/05/04 */

#include <time.h>
#include "toric.h"
#include "IVectors.h"
#include "eterms.h"
#include "poincare.h"
#include "unipoly.h"
#include "../TmpFactorDir/linalg/Zkernel.h"

#include "CoCoA/verbose.H"  // for VerboseLog
#include "CoCoA/error.H"  // for CoCoA_THROW_ERROR

#include <iostream>

// #ifndef ANNA

// extern object math_error;

// #include "coctypes.h"
// #include "cocdefs.h"

// #include "arpoly.h"
// #include "coeffs.h"
// #include "engine.h"
// #include "indet.h"
// #include "monomial.h"
// #include "monopoli.h"
// #include "olist.h"
// #include "omat.h"
// #include "object.h"
// #include "ordering.h"
// #include "polys.h"
// #include "ring.h"
// #include "terms.h"
// #include "weights.h"

#define Print(X) { std::cout << X << std::endl; }


// /*  ELIM  SEQUENTIAL  HDRIVEN  */

// #else  /* ANNA */

// #include "AnnaIO.h"

// #define SEQUENTIAL
// #endif /* ANNA */

/*---------------------------------------------------------------------*/

#define _ANNA_TIME

#define HILB 3000000
#define B_MONITOR 0
#define _PPI
#define _STAT

#define GM
#define SatGM

/*---------------------------------------------------------------------*/

/**************************  global  **************************/

ints INDICES;
ivec ElimVect;

/* order */
static int ( * VECT_WHICH_GT ) (ivec Vect, ivec_elem Index);
static int * COST;

/* to avoid leaks */
static BitermList ReducedGBElems;
/* Statistics */
static int MultiplePairs;
static int InsertedPairs;
static int CoprimePairs;
static int GMPairs;
static int HDPairs;
static int NonSatPairs;

ToricAlgType ToricAlg;

/******************  IVEC  ******************/

static ivec ivec_this_inv(ivec v)
{
  int i;
  for (i=ivec_len(v);i>0;i--) ivec_set_nth(v,i,-ivec_nth(v,i));
  return v;
}

static coc_bool ivec_neg_coprime(ivec v1, ivec v2)
{
  int i;
  for (i=ivec_len(v1) ; i>0 ; i-- )
    if ( ivec_nth(v1,i)<0 && ivec_nth(v2,i)<0 ) return FALSE;
  return TRUE;
}

/******************  BPAIR  ******************/

bpair BPairMake(biterm i, biterm j, int deg)
{
  bpair res =bpair_malloc();

  BPairSet1 (res, i);
  BPairSet2 (res, j);
  BPairSetDeg(res, deg);
  BPairSetLowerDeg(res, !ivec_neg_coprime(i->Vect,j->Vect));

  return res;
}

/******************  BPAIRS  ******************/

bpairs BPsNew(long size)
{
  bpairs res =bpairs_malloc(size);

  BPsSetSize(res, size);
  BPsSetLen(res, 0);

  return res;
}

/******************  BPAIRLIST  ******************/

BPairList BPLNew(int Size)
{
  BPairList res;

  res = BPList_malloc();
  BPLSetLowerMinDegPairs (res, bpairs_malloc(Size));
  BPLSetLowerMinDegSize  (res, Size);
  BPLSetLowerMinDegLen   (res, 0);
  BPLSetMinDegPairs (res, bpairs_malloc(Size));
  BPLSetMinDegSize  (res, Size);
  BPLSetMinDegLen   (res, 0);
  BPLSetPairs (res, bpairs_malloc(Size));
  BPLSetSize  (res, Size);
  BPLSetLen   (res, 0);

  return res;
}

void FreeBPList (BPairList theBPList)
{
  bpairs_free (BPLPairs(theBPList));
  bpairs_free (BPLMinDegPairs(theBPList));
  bpairs_free (BPLLowerMinDegPairs(theBPList));
  BPList_free (theBPList);
}

void EraseMinDegBPList(BPairList theBPList)
{
  int n;
  bpairs BPs1 =BPLMinDegPairs(theBPList),
         BPs2 =BPLLowerMinDegPairs(theBPList);

#ifdef STAT
  for ( n=GetLen(BPs1) ; n>0 ; n--)  { HDPairs++ ; bpair_free(BPs1[n]);}
  for ( n=GetLen(BPs2) ; n>0 ; n--)  { HDPairs++ ; bpair_free(BPs2[n]);}
#else
  for ( n=GetLen(BPs1) ; n>0 ; n--)  bpair_free(BPs1[n]);
  for ( n=GetLen(BPs2) ; n>0 ; n--)  bpair_free(BPs2[n]);
#endif
  BPsSetLen (BPs1, 0);
  BPsSetLen (BPs2, 0);
}

void EraseOtherBPList(BPairList theBPList)
{
  int n;
  bpairs BPs1 =BPLPairs(theBPList);

#ifdef STAT
  for ( n=GetLen(BPs1) ; n>0 ; n--)  { HDPairs++ ; bpair_free(BPs1[n]);}
#else
  for ( n=GetLen(BPs1) ; n>0 ; n--)  bpair_free(BPs1[n]);
#endif
  BPsSetLen (BPs1, 0);
}

void EraseAndFreeBPList(BPairList theBPList)
{
  EraseMinDegBPList(theBPList);
  EraseOtherBPList(theBPList);
  FreeBPList(theBPList);
}

void BPLChangeSize (BPairList theBPList, int NewSize)
{
  BPLSetPairs(theBPList, bpairs_realloc(BPLPairs(theBPList),NewSize) );
  if ( BPLPairs(theBPList)==0 )
  {
    printf("Err: BPLChangeSize, out of memory\n");
    abort();
  }
  BPLSetSize(theBPList, NewSize);
}

void BPLAddBPair (BPairList theBPList, bpair BP)
/* non-minimal degree */
{
  if ( BPLSize(theBPList)==BPLLen(theBPList) )
    BPLChangeSize(theBPList, BPLSize(theBPList)+16);
  /*
  PutLast(BPLPairs(theBPList), BP);
  */
  (BPLPairs(theBPList))[BPLLen(theBPList)+1] = BP;
  BPLSetLen(theBPList, BPLLen(theBPList)+1);
}

void BPLChangeMinDegSize (BPairList theBPList, int NewSize)
{
  BPLSetMinDegPairs(theBPList,
                    bpairs_realloc(BPLMinDegPairs(theBPList),NewSize) );
  if ( BPLMinDegPairs(theBPList)==0 )
  {
    printf("Err: BPLChangeMinDegSize, out of memory\n");
    abort();
  }
  BPLSetMinDegSize(theBPList, NewSize);
}

void BPLChangeLowerMinDegSize (BPairList theBPList, int NewSize)
{
  BPLSetLowerMinDegPairs(theBPList,
                    bpairs_realloc(BPLLowerMinDegPairs(theBPList),NewSize) );
  if ( BPLLowerMinDegPairs(theBPList)==0 )
  {
    //    Print("Err BPLChangeLowerMinDegSize, out of memory\n");
    abort();
  }
  BPLSetLowerMinDegSize(theBPList, NewSize);
}

void BPLAddMinDegBPair(BPairList theBPList, bpair BP)
{
  if ( BPairLowerDeg(BP) )
  {
    if ( BPLLowerMinDegSize(theBPList)==BPLLowerMinDegLen(theBPList) )
      BPLChangeLowerMinDegSize(theBPList, BPLLowerMinDegSize(theBPList)+16);
    BPLSetLowerMinDegLen(theBPList, BPLLowerMinDegLen(theBPList)+1);
    (BPLLowerMinDegPairs(theBPList))[BPLLowerMinDegLen(theBPList)] = BP;
  }
  else
  {
    if ( BPLMinDegSize(theBPList)==BPLMinDegLen(theBPList) )
      BPLChangeMinDegSize(theBPList, BPLMinDegSize(theBPList)+16);
    BPLSetMinDegLen(theBPList, BPLMinDegLen(theBPList)+1);
    (BPLMinDegPairs(theBPList))[BPLMinDegLen(theBPList)] = BP;
  }
}

void BPListResetMinDeg(BPairList theBPList, int OldMinDeg)
{
  bpairs BPs =BPLPairs(theBPList);
  int n, MinDeg, len =BPLLen(theBPList);

  if ( GetLen(BPs)==0 ) return;

  MinDeg = BPairDeg(BPs[len]);
  n = len;
  OldMinDeg++;
  while ( --n>0 && MinDeg!=OldMinDeg)
    if ( BPairDeg(BPs[n]) < MinDeg ) MinDeg = BPairDeg(BPs[n]);
  for ( n=BPLLen(theBPList) ; n>0 ; n-- )
    if ( BPairDeg(BPs[n])==MinDeg )
    {
      BPLAddMinDegBPair (theBPList, BPs[n]);
      BPs[n] = BPs[len--];
    }
  BPsSetLen(BPs,len);
}

void BPLInsertBPair (BPairList theBPList, bpair BP)
{
  int n;
  bpairs MDBPs, BPs;

  if ( BPLMinDegLen(theBPList)==0 && BPLLowerMinDegLen(theBPList)==0 )
  {
    if ( BPLLen(theBPList)!=0 ) Print("Err BPLInsertBPair: BPLMinDegLen==0");
    /* First Insertion */
    if ( BPairLowerDeg(BP) )
    {
      BPLSetLowerMinDegLen(theBPList, BPLLowerMinDegLen(theBPList)+1);
      (BPLLowerMinDegPairs(theBPList))[BPLLowerMinDegLen(theBPList)] = BP;
    }
    else
    {
      BPLSetMinDegLen(theBPList, BPLMinDegLen(theBPList)+1);
      (BPLMinDegPairs(theBPList))[BPLMinDegLen(theBPList)] = BP;
    }
    return;
  }
  if ( BPairDeg(BP)==BPLMinDeg(theBPList) )
  {
    BPLAddMinDegBPair(theBPList, BP);
    return;
  }
  if ( BPairDeg(BP) > BPLMinDeg(theBPList) )
  {
    if ( BPLSize(theBPList)==BPLLen(theBPList) )
      BPLChangeSize(theBPList, BPLSize(theBPList)+16);
    BPLSetLen(theBPList, BPLLen(theBPList)+1);
    (BPLPairs(theBPList))[BPLLen(theBPList)] = BP;
  }
  else  /*  BPairDeg(BP) < MinDeg  */
  {
    if ( BPLSize(theBPList) <
         BPLLen(theBPList)+BPLMinDegLen(theBPList)+
         BPLLowerMinDegLen(theBPList) )
      BPLChangeSize(theBPList, BPLLen(theBPList)+BPLMinDegLen(theBPList)+
                    BPLLowerMinDegLen(theBPList));
    BPs = BPLPairs(theBPList);
    MDBPs = BPLMinDegPairs(theBPList);
    for ( n=GetLen(MDBPs) ; n>0 ; n-- )
    {
      BPsSetLen(BPs, GetLen(BPs)+1);
      BPs[GetLen(BPs)] = MDBPs[n];
    }
    BPsSetLen(MDBPs, 0);
    MDBPs = BPLLowerMinDegPairs(theBPList);
    for ( n=GetLen(MDBPs) ; n>0 ; n-- )
    {
      BPsSetLen(BPs, GetLen(BPs)+1);
      BPs[GetLen(BPs)] = MDBPs[n];
    }
    BPsSetLen(MDBPs, 0);
    if ( BPairLowerDeg(BP) )
    {
      BPLSetLowerMinDegLen(theBPList, BPLLowerMinDegLen(theBPList)+1);
      (BPLLowerMinDegPairs(theBPList))[BPLLowerMinDegLen(theBPList)] = BP;
    }
    else
    {
      BPLSetMinDegLen(theBPList, BPLMinDegLen(theBPList)+1);
      (BPLMinDegPairs(theBPList))[BPLMinDegLen(theBPList)] = BP;
    }
    return;
  }
}

/******************  ANNA_IVEC   ******************/

#if 0
/* ORDERING FUNCTIONS */

/* CYCLIC REVLEX *///---------- inside #if 0 ------------------------
int VectWhichGT(ivec Vect, ivec_elem Index)
{
  int i=Index;

  while ( ivec_nth(Vect,i)==0 )
    if ( i--==1 )
    {
      i = ivec_len(Vect);
      while ( ivec_nth(Vect,i)==0 )
        if ( i--==Index ) return ANSWER_BOTH;
    }
  if ( ivec_nth(Vect,i) > 0) return ANSWER_SECOND;
  /* if ( ivec_nth(Vect,i) < 0) */ return ANSWER_FIRST;
}

/* NON-CYCLIC REVLEX *///---------- inside #if 0 ------------------------
int VectWhichGT(ivec Vect, ivec_elem Index)
{
  int i;

  if ( ivec_nth(Vect,Index) < 0)  return ANSWER_FIRST;
  if ( ivec_nth(Vect,Index) > 0)  return ANSWER_SECOND;

  for ( i=ivec_len(Vect) ; i>0 ; i-- )
  {
    if ( ivec_nth(Vect,i) < 0)  return ANSWER_FIRST;
    if ( ivec_nth(Vect,i) > 0)  return ANSWER_SECOND;
  }
  return ANSWER_BOTH;
}

/* BLOCK-REVLEX *///---------- inside #if 0 ------------------------
int VectWhichGT(ivec Vect, ivec_elem Index)
{
  int i;

  if ( ivec_nth(Vect,Index) < 0)  return ANSWER_FIRST;
  if ( ivec_nth(Vect,Index) > 0)  return ANSWER_SECOND;

  if ( INDICES!=NULL )
    for ( i=IntsGetLen(INDICES) ; i>0 ; i-- )
    {
      if ( ivec_nth(Vect,INDICES[i]) < 0)  return ANSWER_FIRST;
      if ( ivec_nth(Vect,INDICES[i]) > 0)  return ANSWER_SECOND;
    }
  for ( i=ivec_len(Vect) ; i>0 ; i-- )
  {
    if ( ivec_nth(Vect,i) < 0)  return ANSWER_FIRST;
    if ( ivec_nth(Vect,i) > 0)  return ANSWER_SECOND;
  }
  return ANSWER_BOTH;
}

#endif  //---------- end  #if 0 ------------------------

/* BLOCK-REVLEX-LEX */
int VectWhichGT(ivec Vect, ivec_elem Index)
{
  int i;

  //  fprintf(stderr,"\nIndex = %d ",Index);
  if ( ivec_nth(Vect,Index) < 0)  return ANSWER_FIRST;
  if ( ivec_nth(Vect,Index) > 0)  return ANSWER_SECOND;
  //  fprintf(stderr,"\nINDICES = %x ",INDICES);
  if ( INDICES!=NULL )
    for ( i=IntsGetLen(INDICES) ; i>0 ; i-- )
    {
      //      fprintf(stderr,"INDICES[i] = %d ",INDICES[i]);
      if ( ivec_nth(Vect,INDICES[i]) < 0)  return ANSWER_FIRST;
      if ( ivec_nth(Vect,INDICES[i]) > 0)  return ANSWER_SECOND;
    }
  for ( i=ivec_len(Vect) ; i>0 ; i-- )
  {
    //    fprintf(stderr," i = %d ",i);
    if ( ivec_nth(Vect,i) > 0)  return ANSWER_FIRST;
    if ( ivec_nth(Vect,i) < 0)  return ANSWER_SECOND;
  }
  return ANSWER_BOTH;
}

/* COST + REVLEX */
int VectTestSetWhichGT(ivec Vect, ivec_elem /*Index*/)
{
  int i;
  /*  float sp=0, sn=0; */
  int sp=0, sn=0;

  for ( i=ivec_len(Vect) ; i>0 ; i-- )
    if ( ivec_nth(Vect,i) > 0)  sp += ivec_nth(Vect,i)*COST[i];
    else    if ( ivec_nth(Vect,i) < 0)  sn -= ivec_nth(Vect,i)*COST[i];
  if ( sp > sn )  return ANSWER_FIRST;
  if ( sp < sn )  return ANSWER_SECOND;

  for ( i=ivec_len(Vect) ; i>0 ; i-- )
  {
    if ( ivec_nth(Vect,i) < 0)  return ANSWER_FIRST;
    if ( ivec_nth(Vect,i) > 0)  return ANSWER_SECOND;
  }
  return ANSWER_BOTH;
}

/* TestSet XEL */
int VectTestSetXelWhichGT(ivec Vect, ivec_elem /*Index*/)
{
  int i;

  for ( i=ivec_len(Vect) ; i>0 ; i-- )
  {
    if ( ivec_nth(Vect,i) > 0)  return ANSWER_FIRST;
    if ( ivec_nth(Vect,i) < 0)  return ANSWER_SECOND;
  }
  return ANSWER_BOTH;
}

/* TestSet LEX */
int VectTestSetLexWhichGT(ivec Vect, ivec_elem /*Index*/)
{
  int i;

  for ( i=1 ; i<=ivec_len(Vect) ; i-- )
  {
    if ( ivec_nth(Vect,i) > 0)  return ANSWER_FIRST;
    if ( ivec_nth(Vect,i) < 0)  return ANSWER_SECOND;
  }
  return ANSWER_BOTH;
}

int VectDegree(ivec Vect)
{
  int i=ivec_len(Vect), res=0;

  for ( ; i>0 ; i-- )
    if ( ivec_nth(Vect,i)>0 ) res+=ivec_nth(Vect,i);
  return res;
}

/******************  BITERM  ******************/

void BitermSet (biterm B, eterm Lt, ivec Vect)
{
  B->Lt = Lt;
  B->Vect = Vect;
}

biterm biterm_dup(biterm b)
{
  biterm res =biterm_malloc();
  BitermSet (res, eterm_dup(b->Lt), ivec_dup(b->Vect));
  return res;
}

biterm BitermMake(ivec Vect, ivec_elem Index)
{
  biterm res_b;  // 2017-02-02 JAA??? set to NULL to keep compiler quiet???

  switch  ( VECT_WHICH_GT(Vect, Index) )
  {
  case ANSWER_FIRST:
    res_b = biterm_malloc();
    BitermSet (res_b, ivec_pos2eterm(Vect), Vect);
    break;
  case ANSWER_SECOND:
    res_b = biterm_malloc();
    ivec_this_inv(Vect);
    BitermSet (res_b, ivec_pos2eterm(Vect), Vect);
    break;
  case ANSWER_BOTH:
    res_b = NULL;
    ivec_free(Vect);
    break;
  default:
    res_b = NULL;
    CoCoA_THROW_ERROR("Should never get here", "BitermSet");
  }
  return res_b;
}

biterm BitermReduce (biterm B, biterm BReducer, ivec_elem Index)
{
  ivec V;

  V = ivec_sub(B->Vect, BReducer->Vect);
  biterm_free(B);

  return BitermMake(V, Index);
}

biterm BitermTailReduce (biterm B, biterm BReducer, ivec_elem Index)
{
  ivec V;

  V = ivec_sum(B->Vect, BReducer->Vect);
  biterm_free(B);

  return BitermMake(V, Index);
}

biterm BPLRemoveMin(BPairList BPL, ivec_elem Index, coc_bool LowerDegFirst)
{
  biterm res_b;
  bpair BP;

  if ( BPLLowerMinDegLen(BPL)==0 && BPLMinDegLen(BPL)==0 )
    Print("Err ChooseBPair: BPLMinDegLen==0\n");

  if ( LowerDegFirst )
  {
    if ( BPLLowerMinDegLen(BPL)!=0 )
    {
      /*
      BP = GetLast (BPLLowerMinDegPairs(BPL));
      */
      BP = (BPLLowerMinDegPairs(BPL))[BPLLowerMinDegLen(BPL)];
      BPLSetLowerMinDegLen(BPL, BPLLowerMinDegLen(BPL)-1);
    }
    else
    {
      /*
      BP = GetLast (BPLMinDegPairs(BPL));
      */
      BP = (BPLMinDegPairs(BPL))[BPLMinDegLen(BPL)];
      BPLSetMinDegLen(BPL, BPLMinDegLen(BPL)-1);
    }
  }
  else
    if ( BPLMinDegLen(BPL)!=0 )
    {
      /*
      BP = GetLast (BPLMinDegPairs(BPL));
      */
      BP = (BPLMinDegPairs(BPL))[BPLMinDegLen(BPL)];
      BPLSetMinDegLen(BPL, BPLMinDegLen(BPL)-1);
    }
    else
    {
      /*
      BP = GetLast (BPLLowerMinDegPairs(BPL));
      */
      BP = (BPLLowerMinDegPairs(BPL))[BPLLowerMinDegLen(BPL)];
      BPLSetLowerMinDegLen(BPL, BPLLowerMinDegLen(BPL)-1);
    }
  res_b = BitermMake(SBitermVect(BP), Index);
  bpair_free(BP);

  return res_b;
}

biterm BsRemoveMin(biterms Bs)
/* TODO : improve */
{
  int n =1;
  biterm res_b =Bs[1];

  for ( ; n<GetLen(Bs) ; n++ )  Bs[n] = Bs[n+1];
  BsSetLen(Bs, GetLen(Bs)-1);

  return res_b;
}

/******************  BITERMLIST  ******************/

BitermList BLNew(int len, int indetsNo)
{
  BitermList res;
  int d;

  res = BList_malloc();
  SetBiterms     (res, biterms_malloc(len));
  SetReducers    (res, biterms_malloc(len));
  BsSetSize (Reducers(res),len);
  BsSetLen  (Reducers(res),0);
  BListSetSize   (res, len);
  BListSetLen    (res, 0);
  BListSetMaxDeg (res, 0);
  SetReducersLabels (res, ints_init(indetsNo));
  for ( d=indetsNo ; d>=0 ; d--)  SetReducersLabelNth(res, d, 0);

  return res;
}

void EraseAndFreeBList (BitermList theBList)
{
  int i =BListLen(theBList);
  biterm B;
  eterm t;

  for ( ; i>0 ; i-- )
    /*
    biterm_free((Biterms(theBList))[i]);
    */
  {
    B = (Biterms(theBList))[i];
    {
      t = B->Lt;
      if ( t!=NULL)  eterm_free(t);
      if (B->Vect!=NULL)
        ivec_free(B->Vect);
      biterm_struct_free(B);
      B=NULL;
    }
  }

  biterms_free(Biterms(theBList));
  biterms_free(Reducers(theBList));
  ints_free (ReducersLabels(theBList));
  BList_free(theBList);
}

void SetOrderedBiterms (BitermList theBList)
{
  int i, d, Len =BListLen(theBList), MaxIndex, OBsLen =0;
  biterms *BsBucket, OBs =Reducers(theBList), Bs =Biterms(theBList);

  if ( Len==0 )
  {
    fprintf(stderr,"Err: SetOrderedBiterms, Len==0");
    BsSetLen (OBs,0);
    for ( d=biterm_get_indetsNo((Biterms(theBList))[1]) ; d>=0 ; d--)
      SetReducersLabelNth (theBList, d, 0);
    return;
  }
  MaxIndex = biterm_get_indetsNo((Biterms(theBList))[1]);
  BsBucket = (biterms*)flexarray_malloc(MaxIndex, "*biterms");
  SetSize(BsBucket, reinterpret_cast<bitermStruct**>(MaxIndex));

  for ( d=MaxIndex ; d>=0 ; d-- )
  {
    BsBucket[d] = biterms_malloc(Len);
    BsSetSize (BsBucket[d], Len);
    BsSetLen (BsBucket[d], 0);
  }

  for ( i=Len ; i>0 ; i--)
  {
    d = IntsGetLen(Indets(Bs[i]->Lt));
    /*
    PutLast(BsBucket[d], b);
    */
    BsSetLen(BsBucket[d], GetLen(BsBucket[d])+1);
    (BsBucket[d])[GetLen(BsBucket[d])] = Bs[i];
  }

  for ( d=0 ; d<=MaxIndex ; d++)
  {
    for ( i=GetLen(BsBucket[d]); i>0; i-- )  OBs[++OBsLen] = BsBucket[d][i];
    biterms_free (BsBucket[d]);
    SetReducersLabelNth (theBList, d, OBsLen);
  }
  flexarray_free(BsBucket,"*biterms");
  BsSetLen(OBs, OBsLen);
}

void InsOrderedBiterms(BitermList theBList, biterm b)
{
  int i =biterm_get_indetsNo(b), IndNo;
  biterms OBs =Reducers(theBList);
  ints RsLabels = ReducersLabels(theBList);

  BsSetLen(OBs, ++(RsLabels[i]));
  if ( IntsGetLen(Indets(b->Lt)) == i )
    OBs[RsLabels[i]] = b;
  else
  {
    IndNo = IntsGetLen(Indets(b->Lt));
    OBs[RsLabels[i]] = OBs[++(RsLabels[i-1])];
    for ( i-- ; i>IndNo ; i-- )
      OBs[RsLabels[i]] = OBs[++(RsLabels[i-1])];
    OBs[RsLabels[IndNo]] = b;
  }
}

void BListChangeSize  (BitermList theBList, int NewSize)
{
  SetBiterms(theBList, biterms_realloc(Biterms(theBList), NewSize));
  SetReducers(theBList, biterms_realloc(Reducers(theBList), NewSize));

  if ( Biterms(theBList)==0 )
  {
    Print("Err BListChangeSize: out of memory\n");
    abort();
  }
  BListSetSize(theBList, NewSize);
}

void BLOrder(BitermList theBList)
{
  int i, n =BListLen(theBList);
  biterm B;
  biterms Bs =Biterms(theBList);

  for ( ; n>0 ; n-- )
    for ( i=n ; i>0 ; i-- )
      if ( BitermLtDeg(Bs[n]) < BitermLtDeg(Bs[i]) )
      {
        B=Bs[n];  Bs[n]=Bs[i];  Bs[i]=B;
      }
}

void SBitermVectCheck (ivec V)
{
  int n;

  for ( n=ivec_len(V) ; n>0 ; n-- )
    if ( ivec_nth(V,n)!=0 ) return;
  Print("Err SBitermVectCheck: SBiterm is 0");
}

ivec SBitermVect (bpair BP)
{
  ivec aux_v;

  if ( eterm_coprime(BPair1(BP)->Lt, BPair2(BP)->Lt) )
  {
    Print("Err SBitermVect: coprime?!?! here?!?!");
    return NULL;
  }
  aux_v = ivec_sub(BPair1(BP)->Vect, BPair2(BP)->Vect);

  SBitermVectCheck (aux_v);

  return aux_v;
}

/******************  ETERMS  ******************/

ints OrderedEtermsIndicesDegNoBuckets (ints theEtermsDeg,
                                       int BLMaxDeg, ints BucketLen)
/* it does not consider eterms of degree 0 */
{
  int i, d, Len =GetLen(theEtermsDeg);
  ints res;

  for ( d=BLMaxDeg ; d>1 ; d--) BucketLen[d-1] += BucketLen[d];
  res = ints_init(BucketLen[1]);
  IntsSetLen(res, BucketLen[1]);
  for ( i=Len ; i>0 ; i--)
    if ( theEtermsDeg[i]!=0 ) res[(BucketLen[theEtermsDeg[i]])--] = i;

  return res;
}

void BLAddBiterm (BitermList theBList, biterm B)
{
  int BLen =BListLen(theBList);

  BListResize (theBList, BLen+1);
  Biterms(theBList)[++BLen] = B;
  BListSetLen (theBList, BLen);

  if ( BitermLtDeg(B) > BListMaxDeg(theBList) )
    BListSetMaxDeg (theBList, eterm_degree(B->Lt));
  InsOrderedBiterms(theBList, B);
}

/******************         ******************/
/******************  TORIC  ******************/
/******************         ******************/

/*---------------------------------------------------------------------*\

	Part of project CoCoA

\*---------------------------------------------------------------------*/

ints HostenShapiro(BitermList BL)
{
  //  CoCoA::VerboseLog VERBOSE("HostenShapiro");
  int i, n, LenBL=BListLen(BL), BR=0, PosNegMin, Other=0, Pos, Neg,
       ncols =BListIndetsNo(BL);
  biterms Bs = Biterms(BL);
  ints res, BestRow;
  ivec *M;
  coc_bool IsNeg=FALSE;

  M = (ivec*)flexarray_malloc(LenBL, "ivec*");
  SetLen(M, reinterpret_cast<ivec>(LenBL));
  SetSize(M, reinterpret_cast<ivec>(LenBL));
  for ( n=LenBL ; n>0 ; n-- )   M[n] = ivec_dup(Bs[n]->Vect);
  res = ints_init(ncols);

  while ( LenBL!=0 )
  {
    PosNegMin = ncols;
    for ( n=LenBL ; n>0 ; n-- )
    {
      Pos = 0; Neg = 0;
      for (i=ncols;i>0;i--) if (M[n][i]>0) Pos++; else if (M[n][i]<0) Neg++;
      if ( Pos < PosNegMin || (Pos==PosNegMin && Neg>Other) )
      {  BR=n ;  PosNegMin=Pos ;  Other=Neg ;  IsNeg=FALSE; }
      if ( Neg < PosNegMin || (Neg==PosNegMin && Pos>Other) )
      {  BR=n ;  PosNegMin=Neg ;  Other=Pos ;  IsNeg=TRUE; }
    }
    BestRow = M[BR];  M[BR] = M[LenBL];  M[LenBL--] = BestRow;
    for ( i=ncols ; i>0 ; i-- )
      if ( BestRow[i]!=0 )
      {
        for ( n=LenBL ; n>0 ; n-- )  M[n][i]=0;
        if ( (IsNeg && BestRow[i]<0) || ((!IsNeg) && BestRow[i]>0) )
          IntsPutLast(res, i);
      }
  }
  for ( n=GetSize(M) ; n>0 ; n-- ) ivec_free(M[n]);
  flexarray_free(M, "ivec*");

  return res;
}

/******************   REDUCTION   ******************/

biterm LabelEtermSearchReducer(eterm t, BitermList theBList)
{
  int n;
  biterms Bs;

  if ( t == NULL ) return NULL;
  Bs = Reducers(theBList);
  for ( n=ReducersLabelNth(theBList,IntsGetLen(Indets(t))) ; n>0 ; n-- )
    if ( eterm_divides(Bs[n]->Lt, t) )   return Bs[n];
  return NULL;
}

biterm EtermSearchReducer(eterm t, BitermList theBList)
{
  int n;
  biterms Bs =Reducers(theBList);

  for ( n=BListLen(theBList) ; n>0 ; n-- )
    if ( eterm_divides( Bs[n]->Lt , t) )  return Bs[n];
  return NULL;
}

biterm BitermNormalForm(biterm B, BitermList theBList, ivec_elem Index)
{
  biterm Reducer;
  eterm t;
  /* #ifdef ELIM */
  ivec aux_v1, aux_v2;
  biterm aux_b;
  /* #endif   */

  if ( B==NULL ) return NULL;
  while ( (Reducer =LabelEtermSearchReducer(B->Lt,theBList)) != NULL )
  {
    B = BitermReduce (B, Reducer, Index);
    if ( B==NULL ) return NULL;
  }
  t = BitermTl(B);
  while ( (Reducer =LabelEtermSearchReducer(t,theBList)) != NULL )
  {
    B = BitermTailReduce (B, Reducer, Index);
    eterm_free(t);
    if ( B==NULL ) return NULL;
    t = BitermTl(B);
  }
  eterm_free(t);

  /* #ifdef ELIM */
  if ( ToricAlg == EATI )
  {

  if ( B != NULL )
    if ( ivec_nth(B->Vect, Index)!=0 )
    {
      aux_v1 = ivec_sub(B->Vect, ElimVect);
      while (ivec_nth(aux_v1, Index)!=0)
      {
        aux_v2 = ivec_sub(aux_v1, ElimVect);
        ivec_free(aux_v1);
        aux_v1 = aux_v2;
      }
      if ( VECT_WHICH_GT(aux_v1, Index)==ANSWER_FIRST )
      {
        aux_b = BitermMake(aux_v1, Index);
        if ( aux_b!=NULL )
        {
          biterm_free(B);
          B = aux_b;
        }
        else
          biterm_free(aux_b);
      }
      else
        ivec_free(aux_v1);
    }
  }
  /* #endif */

  return B;
}

void CheckMinDegReducedNewGen(biterms RedBs, int MinDeg, ints Indices)
{
  int n =IntsGetLen(Indices), Idx;

  for ( ; n>0 ; n--)
  {
    Idx = Indices[n];
    if ( BitermLtDeg(RedBs[Idx])!=MinDeg )
      fprintf(stderr, "MinDegReducedNewGen WRONG! %d != %d\n",
              BitermLtDeg(RedBs[Idx]), MinDeg);
  }
}

biterms MinDegReducedNewGen(BitermList* GBases, biterm B,
                            ints Indices, int B_Index)
/*  B is freed  */
{
  CoCoA::VerboseLog VERBOSE("MinDegReducedNewGen");
  coc_bool Checking;
  int n, n1, Idx, B_Deg = BitermLtDeg(B);
  biterm MinB =B, aux_B;
  biterms res_Bs;

  res_Bs = biterms_malloc(GetLen(GBases));
  BsSetLen (res_Bs, GetLen(GBases));
  BsSetSize (res_Bs, GetLen(GBases));

  for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
  {
    Idx = Indices[n];
    VERBOSE(120) << "B_Index = " << B_Index << " Idx = " << Idx << std::endl;
    aux_B = NULL;
    if ( Idx!=B_Index )
      aux_B = BitermNormalForm (BitermMake(ivec_dup(MinB->Vect), Idx),
                                GBases[Idx], Idx);
    else
      if ( B_Deg==BitermLtDeg(MinB) )
      {
        aux_B = B;
        B = NULL;
      }
      else
      {
        biterm_free(B);
        B = NULL;
        aux_B = BitermNormalForm (BitermMake(ivec_dup(MinB->Vect), Idx),
                                  GBases[Idx], Idx);
      }
    if ( aux_B==NULL )
    {
      for ( n1=IntsGetLen(Indices) ; n1>n ; n1-- )
        biterm_free(res_Bs[Indices[n1]]);
      if ( B!=NULL ) biterm_free(B);
      biterms_free(res_Bs);
      VERBOSE(120) << "--> 0 " << std::endl;
      return NULL;
    }
    res_Bs[Idx] = aux_B;
    if ( BitermLtDeg(aux_B)<BitermLtDeg(MinB) )
      MinB = aux_B;
  }
  if ( BitermLtDeg(MinB)<B_Deg )
  {
    VERBOSE(120) << "LOWER DEGREE -->  " << BitermLtDeg(MinB) << std::endl;
    Checking = TRUE;
    while ( Checking )
    {
      for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
      {
        Idx = Indices[n];
        VERBOSE(120) << " Idx = " << Idx << std::endl;
        aux_B = NULL;
        if ( BitermLtDeg(res_Bs[Idx])==BitermLtDeg(MinB) )
        {
          Checking = FALSE;
          break;
        }
        aux_B = BitermNormalForm (BitermMake(ivec_dup(MinB->Vect), Idx),
                                  GBases[Idx], Idx);
        biterm_free (res_Bs[Idx]);
        res_Bs[Idx] = aux_B;
        if ( aux_B==NULL )
        {
          for ( n1=IntsGetLen(Indices) ; n1>0 ; n1-- )
            if ( n1 != n ) biterm_free(res_Bs[Indices[n1]]);
          if ( B!=NULL ) biterm_free(B);
          biterms_free(res_Bs);
          VERBOSE(120) << "--> 0 " << std::endl;
          return NULL;
        }
        if ( BitermLtDeg(aux_B)<BitermLtDeg(MinB) )
          MinB = aux_B;
      }
    }
  }
  /*
  CheckMinDegReducedNewGen(res_Bs, BitermLtDeg(MinB), Indices);
  */

  if ( B!=NULL ) biterm_free(B);
  return res_Bs;
}

BitermList PreProcessing(BitermList theBList, ints Indices, int *RedNo)
{
  CoCoA::VerboseLog VERBOSE("PreProcessing(toric)");
  BitermList *GBases, GB;
  long n, IndNo =BListIndetsNo(theBList), Deg_B;
  biterm B, *BsReducedB;

  GBases = (BitermList*)flexarray_malloc(IndNo, "BitermList*");
  SetSize(GBases, reinterpret_cast<BitermList>(IndNo));
  SetLen(GBases, reinterpret_cast<BitermList>(IndNo));
  *RedNo = 0;
  for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
    GBases[Indices[n]] = BLNew(BListLen(theBList), IndNo);
  while ( BListLen(theBList)!=0 )
  {
    B = BsRemoveMin(Biterms(theBList));
    Deg_B = BitermLtDeg(B);
    VERBOSE(95) << "INPUT GENERATOR: deg " << Deg_B << std::endl;
    BsReducedB = MinDegReducedNewGen(GBases, B, Indices, 0);
    if ( BsReducedB != NULL )
    {
      if ( Deg_B > BitermLtDeg(BsReducedB[Indices[1]]) ) (*RedNo)++;
      for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
        BLAddBiterm(GBases[Indices[n]], BsReducedB[Indices[n]]);
      biterms_free(BsReducedB);
    }
  }
  EraseAndFreeBList (theBList);
  for ( n=IntsGetLen(Indices) ; n>1 ; n-- )
    EraseAndFreeBList(GBases[Indices[n]]);
  BLOrder(GBases[Indices[1]]);
  GB = GBases[Indices[1]];
  flexarray_free(GBases, "BitermList*");

  return GB;
}

/******************   NEW GEN   ******************/

bpairs BPsReducedBList (BPairList /*theBPList*/, BitermList theBList,
                        biterm B, ivec_elem Index)
{
  int n, BLen =BListLen(theBList), BPDeg;
  eterm LN =B->Lt;
  biterms Bs =Biterms(theBList);
  biterm aux_b;
  bpair BP;
  bpairs resBPs;

  resBPs = bpairs_malloc(BLen);
  BPsSetSize (resBPs, BLen);
  BPsSetLen (resBPs, 0);
  for ( n=BLen ; n>0 ; n-- )
    if ( eterm_divides(LN,Bs[n]->Lt) )
    {
      BPDeg = BitermLtDeg(Bs[n]);
      aux_b = BitermReduce(biterm_dup(Bs[n]), B, Index);
      if ( aux_b!=NULL )
      {
        BP = BPairMake (B, Bs[n], BPDeg);
        biterm_free (aux_b);
        /*
        PutLast (resBPs, BP);
        */
        resBPs[GetLen(resBPs)+1] = BP;
        BPsSetLen(resBPs, GetLen(resBPs)+1);
      }
      BLAddBiterm (ReducedGBElems, Bs[n]);
      Bs[n] = Bs[BLen--];
      BListSetLen (theBList, BLen);
    }
  if ( GetLen(resBPs)!=0 )  SetOrderedBiterms(theBList);

  return resBPs;
}

coc_bool TailCrit(ivec V1, ivec V2, ints Indices, ivec_elem Index)
{
  int i;

  for ( i=IntsGetLen(Indices) ; i>0 ; i-- )
  /*
  for ( i=1 ; i<=IntsGetLen(Indices) ; i++ )
  */
  {
    if (Indices[i] == Index) return FALSE;
    if ( ivec_nth(V1,Indices[i])<0 && ivec_nth(V2,Indices[i])<0 ) return TRUE;
  }
  return FALSE;
}

void MakeCofactors(BitermList theBList, biterm NewGen,
                   eterms *Cofactor, ints *CofactorDeg, ints *OrdPos)
{
  int n, BLen;
  biterms Bs = Biterms(theBList);
  eterm LN =NewGen->Lt;
  eterms Cofactor2;
  ints Ints, BucketLen, Cofactor2Deg;

  BLen = BListLen(theBList);
  Cofactor2 = eterms_malloc(BLen);
  EsSetLen(Cofactor2, BLen);
  EsSetSize(Cofactor2, BLen);

  BucketLen = ints_init(BListMaxDeg(theBList));
  for ( n=BListMaxDeg(theBList) ; n>0 ; n-- ) BucketLen[n] = 0;
  Cofactor2Deg = ints_init(BLen);
  IntsSetLen(Cofactor2Deg, BLen);
  for ( n=BLen ; n>0 ; n-- )
    if ( eterm_coprime(LN, Bs[n]->Lt) )
    {
#ifdef STAT
      CoprimePairs++;
#endif /* STAT */
      Cofactor2[n] = NULL;
      Cofactor2Deg[n] = 0;
    }
    else
    {
      Cofactor2[n] = special_eterm_colon_dup(Bs[n]->Lt, LN);
      (BucketLen[Cofactor2Deg[n] = eterm_degree(Cofactor2[n])])++;
    }
  Ints = OrderedEtermsIndicesDegNoBuckets(Cofactor2Deg,
                                        BListMaxDeg(theBList), BucketLen);
  ints_free(BucketLen);
  *Cofactor = Cofactor2;
  *CofactorDeg = Cofactor2Deg;
  *OrdPos = Ints;
}

void UpdateBListAndBPList(BitermList theBList, BPairList theBPList,
                          biterm NewGen,
                          ivec_elem Index,  coc_bool IsNewGen)
{
  int n, i, BPLen, BLen;
  bpairs BPs, RedBPs;
  biterms Bs =Biterms(theBList);
  ivec V;
  eterm LN, C2Nn, *Cofactor2;
  ints CofPos, Cofactor2Deg, GMPos;

  BLen = BListLen(theBList);
  if ( BLen!= 0 )
  {
    RedBPs = BPsReducedBList(theBPList, theBList, NewGen, Index);
    MakeCofactors(theBList, NewGen, &Cofactor2, &Cofactor2Deg, &CofPos);
    GMPos = ints_init(BLen);
    for ( n=IntsGetLen(CofPos) ; n>0 ; n-- )
      if ( (C2Nn=Cofactor2[CofPos[n]])!=NULL )
      {
#ifdef GM
        for ( i=IntsGetLen(GMPos) ; i>0 ; i-- )
          if ( eterm_divides( Cofactor2[GMPos[i]],C2Nn ) )
          {
#ifdef STAT
            GMPairs++;
#endif /* STAT */
            break;
          }
        if ( i==0 )  /* ToBeInserted */
#endif /* GM */
        {
          eterm_complete(C2Nn);
          IntsPutLast(GMPos, CofPos[n]);
        }
      }
    for ( n=GetLen(Cofactor2) ; n>0 ; n-- )
      if (Cofactor2[n]!=NULL) eterm_free3(Cofactor2[n]);
    eterms_free(Cofactor2);
    ints_free(CofPos);

    V = NewGen->Vect;
    LN = NewGen->Lt;
    if ( IsNewGen )
    {
      if ( BPLMinDegLen(theBPList)==0 && BPLLowerMinDegLen(theBPList)==0 )
        BPListResetMinDeg(theBPList, 0);
      for ( n=IntsGetLen(GMPos) ; n>0 ; n-- )
        /* #ifdef HDRIVEN */
        if ( ToricAlg==PATI &&
             TailCrit((Bs[GMPos[n]])->Vect, V, INDICES, Index) )
#ifdef STAT
          MultiplePairs++
#endif /* STAT */
          ;
        else
          /* #endif  HDRIVEN */
        {
#ifdef STAT
          InsertedPairs++;
#endif /* STAT */
          BPLInsertBPair(theBPList,
                         BPairMake(NewGen, Bs[GMPos[n]],
                                   Cofactor2Deg[GMPos[n]]+eterm_degree(LN)));
        }
      for ( n=GetLen(RedBPs) ; n>0 ; n-- )
        BPLInsertBPair(theBPList, RedBPs[n]);
    }
    else
    {
      BPLen = BPLLen(theBPList);
      BPListResize(theBPList, BPLen+IntsGetLen(GMPos)+GetLen(RedBPs));
      BPs = BPLPairs(theBPList);
      for ( n=IntsGetLen(GMPos) ; n>0 ; n-- )
        /* #ifdef HDRIVEN */
        if ( ToricAlg==PATI &&
             TailCrit((Bs[GMPos[n]])->Vect, V, INDICES, Index) )
#ifdef STAT
          MultiplePairs++
#endif /* STAT */
          ;
        else
          /* #endif  HDRIVEN */
        {
#ifdef STAT
          InsertedPairs++;
#endif /* STAT */
          BPs[++BPLen] = BPairMake(NewGen, Bs[GMPos[n]],
                                   Cofactor2Deg[GMPos[n]]+eterm_degree(LN));
        }
      for ( n=GetLen(RedBPs) ; n>0 ; n-- ) BPs[++BPLen] = RedBPs[n];
      BPLSetLen(theBPList, BPLen);
    }
    ints_free(GMPos);
    ints_free(Cofactor2Deg);
    bpairs_free(RedBPs);
  }

  BitermSetIsMinimal(NewGen, IsNewGen);
  BLAddBiterm (theBList, NewGen);
}

/* #ifdef HDRIVEN */
TermList BListToTList (BitermList BL)
{
  eterm t;
  TermList auxTL;
  int NewLen =0, n =BListLen(BL);
  biterms Bs =Biterms(BL);
  ints OccInd;

  auxTL = NewTList (n, BListIndetsNo(BL));

  for ( ; n>0 ; n-- )
  {
    t = eterm_dup(Bs[n]->Lt);
    OccInd = Indets(t);
    if ( IntsGetLen(OccInd) == 1 )
    {
      InsInSPList (OccInd[1], eterm_get_nth(t, OccInd[1]), SPList(auxTL) );
      eterm_free (t);
    }
    else
    {
      MTLPutLast(MTList(auxTL), NewLen, t);
      t = NULL;
    }
  }
  TListReduceSize (auxTL, NewLen);
  return auxTL;
}

void NewGenToric(BitermList* GBases, BPairList* theBPLists,
                 biterms BsReducedB, ints Indices)
{
  int n, Idx;

  for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
  {
    Idx = Indices[n];
    UpdateBListAndBPList(GBases[Idx], theBPLists[Idx],
                         BsReducedB[Idx], Idx, TRUE);
  }
  biterms_free(BsReducedB);
}
/* #else */
void NewGenBGrobner(BitermList GBasis, BPairList theBPList, biterm B, int Idx)
{
  UpdateBListAndBPList(GBasis, theBPList, B, Idx, TRUE);
}
/* #endif   HDRIVEN */


/**************************  BGROBNER  **************************/

/***********           TODO: OTTIMIZZAZIONI

memoria:
tempo  :
- implementare heap ? (ordinamento piu' efficiente)

************/

biterm DegGBasis(BPairList theBPList, BitermList GB,
                 int MyDeg,
                 ivec_elem Idx, coc_bool LowerDegFirst)
/****
  Computes the GBasis in degree MyDeg. (pairs. no generators)
  If there is an element of degree < MyDeg, it stops and returns it
****/
{
  biterm B;

#if B_MONITOR
  WriteMonitorDegGBasis(MyDeg,Idx,BPLLowerMinDegLen(theBPList),
                        BPLMinDegLen(theBPList));
#endif
  if ( ( BPLMinDegLen(theBPList)==0 && BPLLowerMinDegLen(theBPList)==0 ))
    Print("Err DegGBasis: MinDegLen(theBPList)==0\n");
  if ( BPLMinDeg(theBPList)!=MyDeg )
    Print("Err DegGBasis: MinDeg(theBPList)!=MyDeg\n");

  while ( BPLMinDegLen(theBPList)!=0 || BPLLowerMinDegLen(theBPList)!=0 )
  {
    B = BPLRemoveMin(theBPList, Idx, LowerDegFirst);
    B = BitermNormalForm (B, GB, Idx);
    if ( B!=NULL )
    {
      if ( BitermLtDeg(B)<MyDeg )  return B;
#if B_MONITOR
      printf("// GENERATOR:  Degree %d  Index %d\n", BitermLtDeg(B), Idx);
#endif
      UpdateBListAndBPList(GB, theBPList, B, Idx, FALSE);
    }
  }
  BPListResetMinDeg(theBPList, MyDeg);

  return NULL;
} /* DegGBasis */

int NextMinDeg(BPairList theBPList, BitermList theBList)
/* TODO: macro */
{
  if ( BPLMinDegLen(theBPList)==0 && BPLLowerMinDegLen(theBPList)==0 )
  {
    if ( BListLen(theBList)==0 )
      return 0;
    else
      return  BitermLtDeg((Biterms(theBList))[1]);
  }
  if ( BListLen(theBList)==0 )
    return  BPLMinDeg(theBPList);
  return MIN(BPLMinDeg(theBPList), BitermLtDeg((Biterms(theBList))[1]));
}

/**************************  ifdef HDRIVEN  **************************/
/* #ifdef HDRIVEN */

biterm DegHDGBasis(BPairList theBPList, BitermList GB,
                   int MyDeg, ints /*Indices*/, ivec_elem Idx,
                   Unipoly PSNum, int /*PSMaxDeg*/)
/****
  Computes the GBasis in degree MyDeg. (pairs. no generators)
  If there is an element of degree < MyDeg, it stops and returns it
****/
{
  int NewGensNo;
  biterm B;
  Unipoly MyPSNum;

  if ( ( BPLMinDegLen(theBPList)==0 && BPLLowerMinDegLen(theBPList)==0 ))
    Print("Err HDDegGBasis: MinDegLen(theBPList)==0\n");
  if ( BPLMinDeg(theBPList)!=MyDeg )
    Print("Err HDDegGBasis: MinDeg(theBPList)!=MyDeg\n");

  /*  MyPSNum = TLTruncPoincareNumerator(BListToTList(GB), MyDeg);*/
  MyPSNum = TLPoincareNumerator(BListToTList(GB));
  NewGensNo = UnipolySubCoeffsOfXd(MyPSNum, PSNum, MyDeg);

  fprintf(stderr, "HD Ord: %d NewGensNo: %d \n",Idx,NewGensNo);

  /* FreeUnipoly(MyPSNum); corretto? */
#if B_MONITOR
  WriteMonitorDegHDGBasis(MyDeg,Idx,NewGensNo,BPLLowerMinDegLen(theBPList),
                          BPLMinDegLen(theBPList));
#endif
  while ( NewGensNo!=0 )
  {
    B = BPLRemoveMin(theBPList, Idx, FALSE);
    B = BitermNormalForm (B, GB, Idx);
    if ( B!=NULL )
    {
      if ( BitermLtDeg(B)<MyDeg )
      {
        if ( NewGensNo==1 ) EraseMinDegBPList (theBPList);
        return B;
      }

#if B_MONITOR
      printf("// GENERATOR:  Degree %d  Index %d\n", BitermLtDeg(B), Idx);
#endif
      UpdateBListAndBPList(GB, theBPList, B, Idx, FALSE);
      NewGensNo--;
    }
  }
  EraseMinDegBPList (theBPList);
  BPListResetMinDeg (theBPList, MyDeg);

  return NULL;
} /* DegHDGBasis */

int NextMinDegAndToDoIndices(BPairList* theBPLists, BitermList theBList,
                             ints Indices, ints ToDoIdxs)
{
  int n =IntsGetLen(Indices), res_mindeg, Idx;

  if ( BListLen(theBList)!=0 )
    res_mindeg = BitermLtDeg((Biterms(theBList))[1]);
  else
    res_mindeg = 0;

  IntsSetLen(ToDoIdxs, 0);
  while ( n>0 &&
          BPLMinDegLen(theBPLists[Indices[n]])==0 &&
          BPLLowerMinDegLen(theBPLists[Indices[n]])==0 )
    n--;
  if ( n==0 ) return res_mindeg;

  Idx = Indices[n];
  if ( res_mindeg==0 )
  {
    IntsPutLast(ToDoIdxs, Idx);
    res_mindeg = BPLMinDeg(theBPLists[Idx]);
  }
  else
    if ( BPLMinDeg(theBPLists[Idx]) <= res_mindeg )
    {
      res_mindeg = BPLMinDeg(theBPLists[Idx]);
      IntsPutLast(ToDoIdxs, Idx);
    }

  n--;
  for ( ; n>0 ; n-- )
  {
    Idx = Indices[n];
    if ( BPLMinDegLen(theBPLists[Idx])!=0 ||
         BPLLowerMinDegLen(theBPLists[Idx])!=0)
    {
      if ( BPLMinDeg(theBPLists[Idx]) == res_mindeg )
        IntsPutLast(ToDoIdxs, Idx);
      else
        if ( BPLMinDeg(theBPLists[Idx]) < res_mindeg )
        {
          res_mindeg = BPLMinDeg(theBPLists[Idx]);
          IntsSetLen(ToDoIdxs, 0);
          IntsPutLast(ToDoIdxs, Idx);
        }
    }
  }

  return res_mindeg;
}

/* from y_n to y_1 */
int GetIdx(ints ToDoIdxs)
{
  return IntsGetLast(ToDoIdxs);
}

/* from y_1 to y_n */
/*
int GetIdx(ints ToDoIdxs)
{
  int i,n;

  i = ToDoIdxs[1];
  for ( n=2 ; n<=IntsGetLen(ToDoIdxs) ; n++ )
    ToDoIdxs[n-1]=ToDoIdxs[n];
  IntsSetLen(ToDoIdxs, IntsGetLen(ToDoIdxs)-1);

  return i;
}
*/

/* #ifdef HDRIVEN */

BitermList HDToric(BitermList theBList, ints Indices)
{
  int MinDeg, n, Idx =0, IndNo, PSMaxDeg =0, PSIdx =0, MinDegLen;
  ints ToDoIdxs;
  biterm B;
  biterms BsReducedB;
  BitermList *GBases, GBasis;
  BPairList *theBPLists;
  Unipoly PSNum =NULL;

/*  fprintf(stderr,"--T--  HDToric %d\n", HILB); */
  MultiplePairs = 0;
  InsertedPairs = 0;
  CoprimePairs = 0;
  GMPairs = 0;
  HDPairs = 0;

  if ( IntsGetLen(Indices)==0 )
  {  ints_free(Indices);  Indices = HostenShapiro(theBList);  }
  INDICES = Indices;

  BLOrder(theBList);
  do  theBList = PreProcessing(theBList, Indices, &n); while ( n!=0 );
  IndNo = BListIndetsNo(theBList);
  StartPoincare(IndNo);

  /*  INIT  */
  GBases    = (BitermList*)flexarray_malloc(IndNo, "BitermList*");
  theBPLists = (BPairList*)flexarray_malloc(IndNo, "BPairList*");
  SetSize(GBases, reinterpret_cast<BitermList>(IndNo));
  SetLen(GBases, reinterpret_cast<BitermList>(IndNo));
  SetSize(theBPLists, reinterpret_cast<BPairList>(IndNo));
  SetLen(theBPLists, reinterpret_cast<BPairList>(IndNo));
  for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
  {
    GBases[Indices[n]]     = BLNew(BListLen(theBList), IndNo);
    theBPLists[Indices[n]] = BPLNew(16*BListLen(theBList));
  }
  ReducedGBElems = BLNew(16, IndNo);
  ToDoIdxs = ints_init(IndNo);

  /*  START  */
  B = NULL;
  MinDeg = BitermLtDeg((Biterms(theBList))[1]);

  while ( MinDeg != 0 )
  {
#if B_MONITOR
    printf("------------------------------------------");
#endif
    if ( IntsGetLen(ToDoIdxs)!=0 && PSMaxDeg < MinDeg )
    {
      Idx = GetIdx(ToDoIdxs);
      B =DegGBasis(theBPLists[Idx], GBases[Idx], MinDeg, Idx, TRUE);
      if ( B==NULL && IntsGetLen(ToDoIdxs)>0 )
      {
        n = IntsGetLen(ToDoIdxs);
        MinDegLen = BPLMinDegLen(theBPLists[ToDoIdxs[n]]);
        for ( ; n>0 ; n-- )
        {
          if ( BPLMinDegLen(theBPLists[ToDoIdxs[n]]) > MinDegLen )
            MinDegLen = BPLMinDegLen(theBPLists[ToDoIdxs[n]]);
          if ( BPLLowerMinDegLen(theBPLists[ToDoIdxs[n]])>MinDegLen )
            MinDegLen = BPLLowerMinDegLen(theBPLists[ToDoIdxs[n]]);
        }
        if ( MinDegLen > HILB )
        {
          if ( PSNum!=NULL ) {FreeUnipoly(PSNum); PSNum=NULL; }
          if (BPLMinDegLen(theBPLists[Idx])==0 &&
              BPLLowerMinDegLen(theBPLists[Idx])==0)
            PSMaxDeg = INT_MAX;
          else
            PSMaxDeg = BPLMinDeg(theBPLists[Idx])-1;
          PSIdx = Idx;
        }
      }
    }
    if ( B==NULL )
      while ( IntsGetLen(ToDoIdxs)!=0 )
      {
        Idx = IntsGetLast(ToDoIdxs);
        if ( PSMaxDeg!=0 &&
             ( BPLLowerMinDegLen(theBPLists[Idx])>HILB ||
               BPLMinDegLen(theBPLists[Idx])>HILB ) )
        {
          if (PSNum==NULL)
          {
            /*
            PSNum=TLTruncPoincareNumerator(BListToTList(GBases[PSIdx]),MinDeg);
            */
            PSNum=TLPoincareNumerator(BListToTList(GBases[PSIdx]));
            fprintf(stderr, "HD Deg: %d\n", MinDeg);
          }
          B = DegHDGBasis(theBPLists[Idx], GBases[Idx], MinDeg,
                          Indices, Idx,
                          PSNum, PSMaxDeg);
        }
        else
          B = DegGBasis(theBPLists[Idx], GBases[Idx], MinDeg,
                        Idx, TRUE);
        if ( B!=NULL )  break;
      }
    if ( B!=NULL )
    {
      if ( PSNum!=NULL ) {FreeUnipoly(PSNum); PSNum=NULL; }
      PSMaxDeg = 0;
#if B_MONITOR
      printf("\n// NEW GENERATOR: Degree %d\n", BitermLtDeg(B));
#endif
      BsReducedB = MinDegReducedNewGen(GBases, B, Indices, Idx);
      if ( BsReducedB==NULL )
        Print("BsRedB==NULL: NEW GENERATOR is in the ideal\n");
      NewGenToric(GBases, theBPLists, BsReducedB, Indices);
    }
    else
    {
      while ( BListLen(theBList)!=0 &&
              BitermLtDeg((Biterms(theBList))[1])==MinDeg )
      {
        B = BsRemoveMin(Biterms(theBList));
#if B_MONITOR
        printf("\n// INPUT GENERATOR: Degree %d\n", BitermLtDeg(B));
#endif
        BsReducedB = MinDegReducedNewGen(GBases, B, Indices, 0);
        if ( BsReducedB != NULL )
        {
          if ( BitermLtDeg(BsReducedB[Indices[1]]) < MinDeg )
          {
            NewGenToric(GBases, theBPLists, BsReducedB, Indices);
            break;
          }
          else
            NewGenToric(GBases, theBPLists, BsReducedB, Indices);
        }
      }
    }
    B = NULL;
#ifdef BB_MONITOR
    for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
    {
      Idx = Indices[n];
      printf("\n/**************/\nGBASES[%d] := [\n", Idx);
      WriteBList(GBases[Idx]);
      printf("0];");
    }
#endif
    MinDeg = NextMinDegAndToDoIndices(theBPLists, theBList, Indices, ToDoIdxs);
  }

#ifdef STAT
  fprintf(stderr, "--T-- CoprimePairs:  %d\n", CoprimePairs);
  fprintf(stderr, "--T-- GMPairs:       %d\n", GMPairs);
  fprintf(stderr, "--T-- TailPairs:     %d\n", MultiplePairs);
  fprintf(stderr, "--T-- InsertedPairs: %d\n", InsertedPairs);
  fprintf(stderr, "--T-- HDPairs:       %d\n", HDPairs);
#endif

#ifdef PPI
  BLOrder(GBases[Indices[1]]);
#endif

  GBasis = GBases[Indices[1]];
  /*  FREE  */
  for ( n=IntsGetLen(Indices) ; n>0 ; n-- ) FreeBPList(theBPLists[Indices[n]]);
  flexarray_free(theBPLists,"BPairList*");
  EraseAndFreeBList (theBList);
  ints_free(ToDoIdxs);
  for ( n=IntsGetLen(Indices) ; n>1 ; n-- )
    EraseAndFreeBList(GBases[Indices[n]]);
  flexarray_free(GBases,"BitermList*");
#if B_MONITOR
  printf("\n// ReducedGBElems: %d\n", BListLen(ReducedGBElems));
#endif
  EraseAndFreeBList (ReducedGBElems);
  ints_free(Indices);

  return GBasis;
}

/**************************  ifndef HDRIVEN  **************************/
/* #else */

BitermList BGrobner(BitermList theBList, ivec_elem Index)
{
  int MinDeg, IndNo;
  biterm B, aux_B;
  BitermList GBasis;
  BPairList theBPList;

  BLOrder(theBList);
  IndNo = BListIndetsNo(theBList);

  GBasis = BLNew (BListLen(theBList), IndNo);
  theBPList = BPLNew (16*BListLen(theBList));

  MinDeg = BitermLtDeg((Biterms(theBList))[1]);

  while ( MinDeg != 0 )
  {
    B = NULL;
    if ( ( BPLMinDegLen(theBPList) || BPLLowerMinDegLen(theBPList) ) &&
         BPLMinDeg(theBPList)==MinDeg )
      B = DegGBasis(theBPList, GBasis, MinDeg, Index, TRUE);
    if ( B!=NULL )
    {
#if B_MONITOR
      printf("\n// NEW GENERATOR: Degree %d\n", BitermLtDeg(B));
#endif
      NewGenBGrobner(GBasis, theBPList, B, Index);
    }
    else
    {
      while ( BListLen(theBList)!=0 &&
              BitermLtDeg((Biterms(theBList))[1])==MinDeg )
      {
        B = BsRemoveMin(Biterms(theBList));
#if B_MONITOR
        printf("\n// INPUT GENERATOR: Degree %d\n", BitermLtDeg(B));
#endif
        aux_B = BitermNormalForm (BitermMake(ivec_dup(B->Vect), Index),
                                  GBasis, Index);
        biterm_free(B);
        if ( aux_B!=NULL )
        {
          NewGenBGrobner(GBasis, theBPList, aux_B, Index);
          if ( BitermLtDeg(aux_B)<MinDeg )  break;
        }
      }
    }

    MinDeg = NextMinDeg(theBPList, theBList);
  }
  FreeBPList (theBPList);

  return GBasis;
}

/* #ifdef ELIM */

BitermList ElimToric(BitermList theBList, ints Indices, ints Weights)
{
  int n, s = 0, IndNo;
  ivec aux_v;
  BitermList res_GBasis;
#ifdef ANNA_TIME
  clock_t t1, t2;

  get_cputime(&t1);
  printf("-->> Elim indets %d\n", IntsGetLen(Indices));
#endif /* ANNA_TIME */

  MultiplePairs = 0;
  InsertedPairs = 0;
  CoprimePairs = 0;
  GMPairs = 0;
  HDPairs = 0;

  if ( BListLen(theBList) == 0 )  { ints_free(Indices); return theBList; }
  if ( IntsGetLen(Indices)==0 )
  {  ints_free(Indices);  Indices = HostenShapiro(theBList); }
  IndNo = BListIndetsNo(theBList);
  aux_v = ivec_init(IndNo);
  for ( n=IndNo ; n>0 ; n-- ) ivec_set_nth(aux_v, n, 0);
  for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
  {
    ivec_set_nth(aux_v, Indices[n], Weights[Indices[n]]);
    s += Weights[Indices[n]];
  }
  ivec_set_nth (aux_v, IndNo, -s);
  ElimVect = ivec_dup(aux_v);
  BLAddBiterm (theBList, BitermMake(aux_v, IndNo));
#ifdef ANNA_TIME
  printf("-->> Elim indets %d\n", IntsGetLen(Indices));
#endif /* ANNA_TIME */
  INDICES = ints_dup(Indices);
  IntsSetLen(Indices,1);
  Indices[1] = IndNo;
  Print("--T-- Elim\n");
  BLOrder(theBList);
  do  theBList = PreProcessing(theBList, Indices, &n);  while ( n!=0 );
  ReducedGBElems = BLNew(16, BListIndetsNo(theBList));

  res_GBasis = BGrobner(theBList, IndNo);

#if B_MONITOR
  printf("-->> Len(res_GBasis) = %d", BListLen(res_GBasis));
#endif

#ifdef STAT
  fprintf(stderr, "--T-- CoprimePairs:  %d\n", CoprimePairs);
  fprintf(stderr, "--T-- GMPairs:       %d\n", GMPairs);
  fprintf(stderr, "--T-- InsertedPairs: %d\n", InsertedPairs);
#endif

  /*  FREE  */
  EraseAndFreeBList(theBList);
  EraseAndFreeBList(ReducedGBElems);
  ivec_free(ElimVect);
  ints_free(INDICES);
  ints_free(Indices);
#ifdef ANNA_TIME
  get_cputime(&t2);
  printf("-->> Elim time = %7.2f\n", (double)(cputime_diff(t2,t1))/1000);
#endif /* ANNA_TIME */
  return res_GBasis;
}

/* #endif   ELIM */
/* #ifdef SEQUENTIAL */

BitermList SequentialToric(BitermList theBList, ints Indices)
{
  int n;
  BitermList res_GBasis =0;
#ifdef ANNA_TIME
#if !defined(__MINGW32__) && !defined(_WINDOWS)
  double d1;
#ifdef TIME_SECONDS
  time_t t1, t2;
#else
  clock_t t1, t2;
#endif
#endif
#endif

  CoCoA::VerboseLog VERBOSE("SequentialToric(BL)");
  MultiplePairs = 0;
  InsertedPairs = 0;
  CoprimePairs = 0;
  GMPairs = 0;

  if ( BListLen(theBList) == 0 )
  {
    ints_free(Indices);
    return theBList;
  }
  if ( IntsGetLen(Indices)==0 )
  { 
    ints_free(Indices);
    Indices = HostenShapiro(theBList);  
    VERBOSE(90) << IntsGetLen(Indices) << " indices:" << std::endl;
    for (long i=1; i<=IntsGetLen(Indices); ++i)
      VERBOSE(90) << Indices[i] << ", " << std::endl;
    VERBOSE(90) << "------------" << std::endl;
  }
  INDICES = Indices;  /*  INDICES=NULL prima di terminare */
  BLOrder(theBList);
  do
  {
    theBList = PreProcessing(theBList, Indices, &n);
    VERBOSE(95) << "PreProcessing: n = " << n << std::endl;
  }
  while ( n!=0 );
  for ( n=IntsGetLen(Indices) ; n>0 ; n-- )
  {
    ReducedGBElems = BLNew(16, BListIndetsNo(theBList));
#ifdef ANNA_TIME
#if !defined(__MINGW32__) && !defined(_WINDOWS)
#ifdef TIME_SECONDS
  t1 = time(&t1);
#else
  t1 = clock();
#endif
#endif
#endif
    res_GBasis = BGrobner(theBList, Indices[n]);
#ifdef ANNA_TIME
#if !defined(__MINGW32__) && !defined(_WINDOWS)
#ifdef TIME_SECONDS
    t2 = time(&t2);  d1 = difftime(t2,t1);
#else
    t2 = clock();  t1 = t2-t1; d1=(double)t1/1000000;// DENOM SHOULD BE CLOCKS_PER_SEC
#endif
#if B_MONITOR
    fprintf(stderr,"// time: %4.2f  sec  \t", d1);
    fprintf(stderr,"|    %d :  len = %d (%d)\n",
            Indices[n], BListLen(res_GBasis), MinimalLen(res_GBasis));
    printf("// ReducedGBElems: %d\n", BListLen(ReducedGBElems));
#endif
#endif
#endif
    /*  FREE  */
    EraseAndFreeBList (ReducedGBElems);
    EraseAndFreeBList (theBList);

    theBList = res_GBasis;
  }

#ifdef STAT
  fprintf(stderr, "--T-- CoprimePairs:  %d\n", CoprimePairs);
  fprintf(stderr, "--T-- GMPairs:       %d\n", GMPairs);
  fprintf(stderr, "--T-- InsertedPairs: %d\n", InsertedPairs);
#endif

#ifdef PPI
  BLOrder(res_GBasis);
#endif

  /*  FREE  */
  ints_free(Indices);
  INDICES = NULL;

  return res_GBasis;
}

/* #endif   SEQUENTIAL */
/* #endif   HDRIVEN */


int null_space(int ***basis, int **M, int nrows, int ncols)
{
  int trouble;
  int dim, **ans;
  Zmat Mbig, Zbasis;
  int i, j;

  Mbig = Zmat_ctor(nrows, ncols);
  for (i=0; i < nrows; i++)
    for (j=0; j< ncols; j++)
      mpz_set_si(Mbig->entry[i][j], M[i][j]);

  dim = Zmat_left_kernel_Zbasis(&Zbasis, Mbig);

  /* Convert answer to small integers, set trouble if conversion "failed" */
  trouble = 0;
  ans = (int**)malloc(dim*sizeof(int*));
  for (i=0; i < dim; i++)
  {
    ans[i] = (int*)malloc(nrows*sizeof(int));
    for (j=0; j < nrows; j++)
    {
      ans[i][j] = mpz_get_si(Zbasis->entry[i][j]);
      if (mpz_sizeinbase(Zbasis->entry[i][j], 2) > 8*sizeof(int)-1) trouble = 1;
    }
  }

  Zmat_dtor(Mbig);
  Zmat_dtor(Zbasis);
  *basis = ans;
  if (trouble) return -dim;
  return dim;
}

// #ifndef ANNA

// biterm PPs2Binom(term t1, term t2)
// {
//   ivec Vect;
//   ivec_elem i=term_get_indetsNo(t1);
//   int ord;

//   if ( ToricAlg == EATI )
//   {
//     Vect = ivec_init(i+1);
//     ivec_set_nth(Vect, i+1, 0);
//     ord = i+1;
//   }
//   else
//   {
//     Vect = ivec_init(i);
//     ord = 1;
//   }
//   for ( ; i>0 ; i-- )
//     ivec_set_nth(Vect, i, term_get_nth(t1,i)-term_get_nth(t2,i));
//   if ( term_degree(t1)!=term_degree(t2) )
//   {
//     ivec_free(Vect);
//     return NULL;
//   }
//   return BitermMake(Vect, ord);
// }

// BitermList olist2BList(olist L, ring R)
// { /* assumed L to be a list of binomials */
//   int n, i, MaxDeg = 0, BLLen = 0;
//   monopoli aux_m;
//   BitermList BL;
//   biterms Bs;
//   arpoly aux_p;
//   object aux_o;
//   coeff c1, c2, c;
// #ifdef ANNA_TIME
//   clock_t t1, t2;

//   get_cputime(&t1);
// #endif /* ANNA_TIME */

//   INDICES = NULL;
//   n = olist_len(L);

//   if ( ToricAlg == EATI )
//     BL = BLNew(n, R->indetsNo+1);
//   else
//     BL = BLNew(n, R->indetsNo);
//   Bs = Biterms(BL);
//   for ( i=1 ; i<=n ; i++ )
//   {
//     aux_o = olist_nth(L,i);
//     aux_p = object_cast_arpoly(object_link(aux_o), R);
//     if ( math_error != object_NULL )
//     {
//       EraseAndFreeBList(BL);
//       return NULL;
//     }
//     aux_m = this_arpoly2monopoli(aux_p,R);
//     if ( aux_m != NULL )
//     {
//       if ( monomial_next(aux_m) == NULL
//            || monomial_next(monomial_next(aux_m)) != NULL )
//       {
//         EraseAndFreeBList(BL);
//         monopoli_free(aux_m,R);
//         math_error = cocoa_error("BAD_PARAMS");
//         return NULL;
//       }
//       c1 = monomial_coeff(aux_m);
//       c2 = monomial_coeff(monomial_next(aux_m));
//       c = coeff_sum(c1, c2, R);
//       if ( !(coeff_is_zero(c, R) && (coeff_is_one(c1,R)||coeff_is_one(c2,R))))
//       {
//         EraseAndFreeBList(BL);
//         monopoli_free(aux_m,R);
//         coeff_free(c, R);
//         math_error = cocoa_error("BAD_PARAMS");
//         return NULL;
//       }
//       coeff_free(c, R);
//       Bs[++BLLen] = PPs2Binom(monomial_term(aux_m),
//                               monomial_term(monomial_next(aux_m)));
//       monopoli_free(aux_m, R);
//       if ( Bs[BLLen]==NULL )
//       {
//         EraseAndFreeBList(BL);
//         math_error = cocoa_error("BAD_PARAMS");
//         return NULL;
//       }
//       BListSetLen(BL, BLLen);
//       if ( BitermLtDeg(Bs[BLLen]) > MaxDeg)  MaxDeg = BitermLtDeg(Bs[BLLen]);
//     }
//   }
//   BListSetMaxDeg (BL, MaxDeg);
// #ifdef ANNA_TIME
//   get_cputime(&t2);
//   printf("-->> olist2BList time = %7.2f\n",(double)(cputime_diff(t2,t1))/1000);
// #endif /* ANNA_TIME */

//   return BL;
// }

// ints ring_weilist2Ints(ring R)
// {
//   int n = R->indetsNo;
//   weight w = weilist_nth(ring_weilist(R),1);
//   ints res_w = ints_init(n);

//   IntsSetLen(res_w,n);
//   for (  ; n>0 ; n-- )  res_w[n] = weight_nth(w,n);

//   return res_w;
// }

// ints ring_2ordrow2Ints(ring R)
// {
//   int n = R->indetsNo;
//   ordering o = ring_ordering(R);
//   ints res = ints_init(n);

//   IntsSetLen(res,n);
//   for (  ; n>0 ; n-- )  res[n] = ordering_elem(o,2,n);

//   return res;
// }

// ints olist2Ints(olist L, ring R)
// {
//   int i, idx;
//   arpoly aux_p;
//   object aux_o;
//   ints Indices = ints_init(R->indetsNo);
//   for ( i=1 ; i<=olist_len(L) ; i++ )
//   {
//     aux_o = olist_nth(L,i);
//     aux_p = object_cast_arpoly(object_link(aux_o), R);
//     if ( math_error != object_NULL )
//     {
//       ints_free(Indices);
//       return NULL;
//     }
//     if ( aux_p != NULL)
//     {
//       idx = arpoly2indet(aux_p,R);
//       arpoly_free(aux_p,R);
//       if ( idx == unknown_NULL )
//       {
//         ints_free(Indices);
//         math_error = cocoa_error("BAD_PARAMS");
//         return NULL;
//       }
//       else
//         IntsPutLast(Indices, idx+1);
//     }
//   }
//   return Indices;
// }

// arpoly WBinom2arpoly(biterm B, ints Weights, ints ring_weights, ring R)
// {
//   ivec Vect = B->Vect;
//   ivec_elem i = ivec_len(Vect);
//   term t1, t2;
//   arpoly aux_p;
//   t1 = term_init(R->indetsNo);
//   t2 = term_init(R->indetsNo);

//   if ( ToricAlg == EATI )
//     i--;
//   for ( ; i>0 ; i-- )
//     if ( ivec_nth(Vect,i)>0 )
//       term_put_nth(t1, i, ring_weights[i]*ivec_nth(Vect,i)/Weights[i]);
//     else
//       term_put_nth(t2, i, -ring_weights[i]*ivec_nth(Vect,i)/Weights[i]);
//   aux_p = arpoly_this_term2arpoly(t1,R);
//   aux_p = arpoly_sub_these(aux_p,arpoly_this_term2arpoly(t2,R),R);

//   return aux_p;
// }

// arpoly Binom2arpoly(biterm B, ring R)
// {
//   ivec Vect = B->Vect;
//   ivec_elem i = ivec_len(Vect);
//   term t1, t2;
//   arpoly aux_p;

//   t1 = term_init(R->indetsNo);
//   t2 = term_init(R->indetsNo);

//   if ( ToricAlg == EATI )
//     i-- ;
//   for ( ; i>0 ; i-- )
//     if ( ivec_nth(Vect,i)>0 )
//       term_put_nth(t1, i, ivec_nth(Vect,i));
//     else
//       term_put_nth(t2, i, -ivec_nth(Vect,i));
//   aux_p = arpoly_this_term2arpoly(t1,R);
//   aux_p = arpoly_sub_these(aux_p,arpoly_this_term2arpoly(t2,R),R);

//   return aux_p;
// }

// olist this_BList2olist(BitermList BL, ints Weights, ring R)
// {
//   biterms Bs = Biterms(BL);
//   int i, j=0, n = BListLen(BL), IndNo;
//   olist res_ol = olist_init(n);
//   arpoly p;

//   if ( n != 0 )
//   {
//     IndNo = BListIndetsNo(BL);

//     for ( i=1 ; i<=n ; i++ )
//       if ( ToricAlg != EATI || ivec_nth(Bs[i]->Vect, IndNo)==0 )
//       {
//         j++;
//         p = Binom2arpoly(Bs[i], R);
//         olist_set_nth(res_ol, j, arpoly2object(p,R));
//       }
//     olist_set_len(res_ol, j);
//   }
//   EraseAndFreeBList(BL);
//   if ( Weights != NULL ) ints_free(Weights);

//   return res_ol;
// }

// olist this_WBList2olist(BitermList BL, ints Weights, ring R)
// {
//   biterms Bs = Biterms(BL);
//   int i, j=0, n = BListLen(BL), IndNo;
//   olist res_ol = olist_init(n);
//   arpoly p;
//   ints ring_weights;

//   if ( n != 0 )
//   {
//     ring_weights = ring_weilist2Ints(R);
//     IndNo = BListIndetsNo(BL);
//     for ( i=1 ; i<=n ; i++ )
//       if ( ToricAlg != EATI || ivec_nth(Bs[i]->Vect, IndNo)==0 )
//       {
//         j++;
//         p = WBinom2arpoly(Bs[i], Weights, ring_weights, R);
//         olist_set_nth(res_ol, j, arpoly2object(p,R));
//       }
//     olist_set_len(res_ol, j);
//     ints_free(ring_weights);
//   }
//   EraseAndFreeBList(BL);
//   ints_free(Weights);

//   return res_ol;
// }


// int** smallint_matrix(int kerdim, int nrows, mpz_t **Mbig, int *trouble)
// {
//   int **Msmall, i, j;

//   /* Convert answer to small integers, set trouble if conversion "failed" */
//   *trouble = 0;
//   Msmall = (int**)malloc(kerdim*sizeof(int*));
//   for (i=0; i < kerdim; i++)
//   {
//     Msmall[i] = (int*)malloc(nrows*sizeof(int));
//     for (j=0; j < nrows; j++)
//     {
//       Msmall[i][j] = mpz_get_si(Mbig[i][j]);
//       if (mpz_sizeinbase(Mbig[i][j], 2) > 8*sizeof(int)-1)  *trouble = 1;
//     }
//   }
//   return Msmall;
// }

// #if 0   /* test speed of null_space before deleting this */

// mpz_t** matrix_copy(mpz_t** matrix, int nrows, int ncols)
// {
//   mpz_t **res;
//   int i,j;

//   res = (mpz_t**)MALLOC(nrows*sizeof(mpz_t*));
//   for (i=0; i < nrows; i++)
//   {
//      res[i] = (mpz_t*)MALLOC(ncols*sizeof(mpz_t));
//      for (j=0 ; j<ncols ; j++)  mpz_init_set(res[i][j], matrix[i][j]);
//   }
//   return res;
// }

//   /*
//   aux_matrix = matrix_copy(matrix, nrows, ncols);
//   dim = Zmat_left_kernel_Zbasis(&basis, aux_matrix, nrows, ncols);
//   small_basis = smallint_matrix(dim, nrows, basis, &trouble);
//   for (i=0 ; i<nrows ; i++)
//   {
//     for (j=0;j<nrows;j++) mpz_clear(basis[i][j]);  free(basis[i]);
//     for (j=0;j<ncols;j++) mpz_clear(aux_matrix[i][j]); free(aux_matrix[i]);
//   }
//   free(aux_matrix);
//   free(basis);
//   */

// #endif

// void omat_ker2BList(omat M, BitermList *BL, ints *Indices, ints *Weights)
// {
//   mpz_t **matrix;
//   mpz_t D, gcd, tmp2;
//   mpq_t tmpq;
//   arbrat tmp;
//   int dim, nrows, ncols, i, j;
//   int trouble;
//   int **small_matrix, **small_basis, matrices =8, CurrLen =0, IndNo;
//   ints Indices1, weights;
//   biterm B, MaxDegB;
//   biterms Bs;
//   ivec aux_V;
//   int Index, MaxDeg=0;

//   /* We shall automatically transpose the input matrix!! */
//   nrows = omat_ncol(M);
//   ncols = omat_nrow(M);

//   if ( ToricAlg == EATI )
//     IndNo = nrows+1;
//   else
//     IndNo = nrows;

//   /*  TODO CONTROLLA!!!!!!!!!!!!!!!!!!!!!!!!!!!
//   rum_init(sizeof(expType)*(IndNo+1), RUM_STD_SIZE);
//   rum_init(sizeof(expType)*(IndNo+2), RUM_STD_SIZE);
//   */

//   /* INDICES */
//   Indices1 = ints_malloc(ncols);
//   IntsSetSize(Indices1, ncols);
//   IntsSetLen(Indices1, 0);
//   *Indices = Indices1;
//   INDICES = Indices1;

//   /* OMAT TO MAT */
//   mpz_init_set_ui(D, 1);
//   mpz_init(tmp2);
//   mpz_init(gcd);
//   mpq_init(tmpq);
//   for (i=0; i < nrows; i++)
//   {
//     for (j=0; j < ncols; j++)
//     {
//       tmp = object_cast_arbrat(object_link(omat_ij(M, j+1, i+1)));
//       if (math_error != object_NULL) return; /* LEAK!!! */
//       mpq_set(tmpq, tmp);
//       mpq_clear(tmp);
//       free_verylongrat(tmp);
//       mpz_gcd(gcd, D, mpq_denref(tmpq));
//       mpz_divexact(tmp2, mpq_denref(tmpq), gcd);
//       mpz_mul(D, D, tmp2);
//     }
//   }
//   matrix = (mpz_t**)MALLOC(nrows*sizeof(mpz_t*));
//   for (i=0; i < nrows; i++)
//   {
//     matrix[i] = (mpz_t*)MALLOC(ncols*sizeof(mpz_t));
//     for (j=0; j < ncols; j++)
//     {
//       tmp = object_cast_arbrat(object_link(omat_ij(M, j+1, i+1)));
//       mpq_set(tmpq, tmp);
//       mpq_clear(tmp);
//       free_verylongrat(tmp);
//       mpz_init(matrix[i][j]);
//       mpz_set(matrix[i][j], mpq_numref(tmpq));
//       mpz_divexact(tmp2, D, mpq_denref(tmpq));
//       mpz_mul(matrix[i][j], matrix[i][j], tmp2);
//     }
//   }
//   mpz_clear(tmp2);
//   mpq_clear(tmpq);
//   mpz_clear(D);
//   mpz_clear(gcd);
//   /* WEIGHTS */
//   weights = ints_init(nrows);
//   IntsSetLen(weights, nrows);

//   small_matrix = smallint_matrix(nrows, ncols, matrix, &trouble);
//   for (i=0 ; i<nrows ; i++)
//   { for (j=0 ; j<ncols ; j++) mpz_clear(matrix[i][j]);  free(matrix[i]); }
//   free(matrix);
//   for (i=0 ; i<nrows ; i++)
//   {  weights[i+1]=0; for (j=0;j<ncols;j++) weights[i+1]+=small_matrix[i][j]; }
//   i = 1;
//   while ( weights[i+1]==weights[i] ) if ((++i)==nrows) break;
//   if ((i)==nrows)
//     for (j=1; j <= nrows; j++) weights[j] =1;
//   else
//   {
//     i = 0;  while ( small_matrix[i][0]!=0 ) if ((++i)==nrows) break;
//     if ((i)==nrows)  for (j=0; j<nrows; j++) weights[j+1] = small_matrix[j][0];
//   }
//   *Weights = weights;

//   /* NULL_SPACE TO BLIST */
//   dim = null_space( &small_basis, small_matrix, nrows, ncols);
//   *BL = BLNew (matrices*dim, IndNo);
//   BListSetLen (*BL, matrices*dim);
//   Bs = Biterms(*BL);
//   for ( i=0 ; i<dim ; i++)
//   {
//     aux_V =ivec_init(IndNo);
//     for ( Index=0 ; Index<IndNo ; Index++ )
//       if ( Index>=nrows )
//         ivec_set_nth(aux_V, Index+1, 0);
//       else
//         ivec_set_nth(aux_V, Index+1, weights[Index+1]*small_basis[i][Index]);
//     B = BitermMake(aux_V,1);
//     Bs[++CurrLen] = B;
//     if ( BitermLtDeg(B) > MaxDeg) { MaxDeg = BitermLtDeg(B); MaxDegB = B; }
//   }
//   for (i=0; i<dim; i++) free(small_basis[i]);
//   free(small_basis);
//   for ( matrices-- ; matrices>0 ; matrices-- )
//   {
//     dim = null_space( &small_basis, small_matrix, nrows, ncols);
//     for ( i=0 ; i<dim ; i++ )
//     {
//       aux_V =ivec_init(IndNo);
//       for ( Index=0 ; Index<IndNo ; Index++ )
//         if ( Index>=nrows )
//           ivec_set_nth(aux_V, Index+1, 0);
//         else
//           ivec_set_nth(aux_V, Index+1, weights[Index+1]*small_basis[i][Index]);
//       B = BitermMake(aux_V,1);
//       Bs[++CurrLen] = B;
//       if ( BitermLtDeg(B) > MaxDeg) { MaxDeg = BitermLtDeg(B); MaxDegB = B; }
//     }
//     for (i=0; i<dim; i++) free(small_basis[i]);
//     free(small_basis);
//   }
//   for (i=0; i<nrows ; i++) free(small_matrix[i]);
//   free(small_matrix);
//   BListSetMaxDeg (*BL, MaxDeg);

//   /*
//   if (trouble) return -dim;
//   return dim;
//   */
// }

// #endif

/******************  PROGRAMMAZIONE INTERA  ******************/

biterm SaturatedBitermNormalForm(biterm B, BitermList theBList)
{
  biterm Reducer;
  eterm t;

  if ( B==NULL ) return NULL;
  while ( (Reducer =LabelEtermSearchReducer(B->Lt,theBList)) != NULL )
  {
    if (!ivec_neg_coprime(Reducer->Vect,B->Vect))
    {
      biterm_free(B);
      return NULL;
    }
    B = BitermReduce (B, Reducer, 1);
    if ( B==NULL ) return NULL;
    /*
    if ( eterm_divides(Reducer->Lt, B->Lt) )
    {
      fprintf(stderr, " %d ", IntsGetLen(Indets(Reducer->Lt)));
      B = BitermReduce (B, Reducer, 1);
      if ( B==NULL ) return NULL;
      while ( eterm_divides(Reducer->Lt, B->Lt) )
      {
        Print(".");
        B = BitermReduce (B, Reducer, 1);
        if ( B==NULL ) return NULL;
      }
    }
    */
  }
  t = BitermTl(B);
  while ( (Reducer =LabelEtermSearchReducer(t,theBList)) != NULL )
  {
    B = BitermTailReduce (B, Reducer, 1);
    eterm_free(t);
    if ( B==NULL ) return NULL;
    t = BitermTl(B);
  }
  eterm_free(t);

  return B;
}

coc_bool TailsNotCoprime(ints NegIdxs, ivec v)
{
  int i;

  for (i=IntsGetLen(NegIdxs) ; i>0 ; i-- )
    if ( ivec_nth(v, NegIdxs[i]) < 0 ) return TRUE;

  return FALSE;
}

void SaturatedMakeCofactors(BitermList theBList, biterm NewGen,
                            eterms *Cofactor, ints *CofactorDeg, ints *OrdPos)
{
  int n, BLen;
  biterms Bs = Biterms(theBList);
  eterm LN =NewGen->Lt;
  ivec VN =NewGen->Vect;
  eterms Cofactor2;
  ints Ints, BucketLen, Cofactor2Deg, NegIdxs;

  NegIdxs = ints_init(ivec_len(VN));
  for ( n=ivec_len(VN) ; n>0 ; n-- )
    if ( ivec_nth(VN, n) < 0 ) IntsPutLast(NegIdxs, n);

  BLen = BListLen(theBList);
  Cofactor2 = eterms_malloc(BLen);
  EsSetLen(Cofactor2, BLen);
  EsSetSize(Cofactor2, BLen);

#ifdef GM
  BucketLen = ints_init(BListMaxDeg(theBList));
  for ( n=BListMaxDeg(theBList) ; n>0 ; n--) BucketLen[n] = 0;
#endif
  Cofactor2Deg = ints_init(BLen);
  IntsSetLen(Cofactor2Deg, BLen);
  for ( n=BLen ; n>0 ; n-- )
    if ( eterm_coprime(LN,Bs[n]->Lt) || TailsNotCoprime(NegIdxs,Bs[n]->Vect) )
    {
#ifdef STAT
      if ( eterm_coprime(LN, Bs[n]->Lt))
        CoprimePairs++;
      else
        NonSatPairs++;
#endif
      Cofactor2[n] = NULL;
      Cofactor2Deg[n] = 0;
    }
    else
    {
      Cofactor2[n] = special_eterm_colon_dup(Bs[n]->Lt, LN);
#ifdef GM
      (BucketLen[Cofactor2Deg[n] = eterm_degree(Cofactor2[n])])++;
#endif
    }
#ifdef GM
  Ints = OrderedEtermsIndicesDegNoBuckets(Cofactor2Deg,
                                          BListMaxDeg(theBList), BucketLen);
  ints_free(BucketLen);
#endif
  ints_free(NegIdxs);
  *Cofactor = Cofactor2;
  *CofactorDeg = Cofactor2Deg;
#ifdef GM
  *OrdPos = Ints;
#endif
}

void SaturatedUpdateBListAndBPList(BitermList theBList,
                                   BPairList theBPList,
                                   biterm NewGen, coc_bool IsNewGen)
{
  int n, i, BPLen, BLen;
  bpairs BPs;
  biterms Bs =Biterms(theBList);
//  ivec V;  // 2022-03-04 assigned but not used?!?!
  eterm LN, C2Nn, *Cofactor2;
  ints CofPos, Cofactor2Deg, GMPos;

  BLen = BListLen(theBList);
  if ( BLen!= 0 )
  {
    SaturatedMakeCofactors(theBList, NewGen,
                           &Cofactor2, &Cofactor2Deg, &CofPos);
    GMPos = ints_init(BLen);
#ifdef GM
    for ( n=IntsGetLen(CofPos) ; n>0 ; n-- )
      if ( (C2Nn=Cofactor2[CofPos[n]])!=NULL )
      {
        for ( i=IntsGetLen(GMPos) ; i>0 ; i-- )
          if ( eterm_divides( Cofactor2[GMPos[i]],C2Nn ) )
          {
#ifdef STAT
            GMPairs++;
#endif /* STAT */
            break;
          }
        if ( i==0 )  /* ToBeInserted */
        {
          eterm_complete(C2Nn);
          IntsPutLast(GMPos, CofPos[n]);
        }
      }
    ints_free(CofPos);
#else  /*  ifndef GM  */
    for ( n=GetLen(Cofactor2) ; n>0 ; n-- )
      if ( (C2Nn=Cofactor2[n])!=NULL )
      {        eterm_complete(C2Nn); IntsPutLast(GMPos, n); }
#endif  /*  GM  */
    for ( n=GetLen(Cofactor2) ; n>0 ; n-- )
      if (Cofactor2[n]!=NULL) eterm_free3(Cofactor2[n]);
    eterms_free(Cofactor2);

//    V  = NewGen->Vect;  // 2022-03-04 assigned but not used?!?!
    LN = NewGen->Lt;
    if ( IsNewGen )
    {
      if ( BPLLowerMinDegLen(theBPList)!=0 ) fprintf(stderr, "ERR1");
      if ( BPLMinDegLen(theBPList)==0 ) BPListResetMinDeg(theBPList, 0);
      for ( n=IntsGetLen(GMPos) ; n>0 ; n-- )
      {
#ifdef STAT
        InsertedPairs++;
#endif /* STAT */
        BPLInsertBPair(theBPList,
                       BPairMake(NewGen, Bs[GMPos[n]],
                                 Cofactor2Deg[GMPos[n]]+eterm_degree(LN)));
      }
    }
    else
    {
      BPLen = BPLLen(theBPList);
      BPListResize(theBPList, BPLen+IntsGetLen(GMPos));
      BPs = BPLPairs(theBPList);
      for ( n=IntsGetLen(GMPos) ; n>0 ; n-- )
      {
#ifdef STAT
        InsertedPairs++;
#endif /* STAT */
        BPs[++BPLen] = BPairMake(NewGen, Bs[GMPos[n]],
                                 Cofactor2Deg[GMPos[n]]+eterm_degree(LN));
      }
      BPLSetLen(theBPList, BPLen);
    }
    ints_free(GMPos);
    ints_free(Cofactor2Deg);
  }
  if ( BPLLowerMinDegLen(theBPList)!=0 ) fprintf(stderr, "ERR2");
  if ( BPLMinDegLen(theBPList)==0 ) BPListResetMinDeg(theBPList, 0);
  BLAddBiterm (theBList, NewGen);
} /* SaturatedUpdateBListAndBPList */

void SaturatedDegGBasis(BPairList theBPList, BitermList GB, int MyDeg)
/****
  Computes the GBasis in degree MyDeg. (pairs. no generators)
****/
{
  biterm B;

#if B_MONITOR
  WriteMonitorDegGBasis(MyDeg,1, 0, BPLMinDegLen(theBPList));
#endif
  if ( BPLLowerMinDegLen(theBPList)!=0 ) fprintf(stderr, "ERR3");
  if ( BPLMinDegLen(theBPList)==0 )
    Print("Err DegGBasis: MinDegLen(theBPList)==0\n");
  if ( BPLMinDeg(theBPList)!=MyDeg )
    Print("Err DegGBasis: MinDeg(theBPList)!=MyDeg\n");

  /* BPLMinDeg might change in this loop */
  while ( BPLMinDegLen(theBPList) != 0 && BPLMinDeg(theBPList)==MyDeg )
  {
    B = BPLRemoveMin(theBPList, 1, FALSE);
    B = SaturatedBitermNormalForm (B, GB);
    if ( B!=NULL )
    {
#if B_MONITOR
      printf("// GENERATOR:  Degree %d \n", BitermLtDeg(B));
#endif
      SaturatedUpdateBListAndBPList(GB, theBPList, B, FALSE);
    }
  }
  if (BPLMinDegLen(theBPList) == 0)  BPListResetMinDeg(theBPList, MyDeg);
} /* SaturatedDegGBasis */

BitermList SaturatedBGrobner(BitermList theBList, int trunc)
{
  int MinDeg, IndNo;
  biterm B, aux_B;
  BitermList GBasis;
  BPairList theBPList;

  IndNo = BListIndetsNo(theBList);

  GBasis = BLNew (BListLen(theBList), IndNo);
  theBPList = BPLNew (16*BListLen(theBList));

  MinDeg = BitermLtDeg((Biterms(theBList))[1]);

  while ( MinDeg != 0 )
  {
    if ( ( MinDeg > trunc ) && ( trunc != 0 ) )
    { EraseAndFreeBPList(theBPList); return GBasis; }
    B = NULL;
    if ( BPLLowerMinDegLen(theBPList)!=0 ) fprintf(stderr, "ERR4");
    if ( BPLMinDegLen(theBPList)!=0 && BPLMinDeg(theBPList)==MinDeg )
      SaturatedDegGBasis(theBPList, GBasis, MinDeg);
    while ( BListLen(theBList)!=0 &&
            BitermLtDeg((Biterms(theBList))[1])==MinDeg )
    {
      B = BsRemoveMin(Biterms(theBList));
#if B_MONITOR
      printf("\n// INPUT GENERATOR: Degree %d\n", BitermLtDeg(B));
#endif
      aux_B = SaturatedBitermNormalForm (BitermMake(ivec_dup(B->Vect), 1),
                                         GBasis);
      biterm_free(B);
      if ( aux_B!=NULL )
        SaturatedUpdateBListAndBPList(GBasis, theBPList, aux_B, TRUE);
    }

    MinDeg = NextMinDeg(theBPList, theBList);
  }
  FreeBPList (theBPList);

  return GBasis;
} /* SaturatedBGrobner */


BitermList TestSet(BitermList theBList, int trunc)
{
  /* assumiamo theBList omogeneo */
  BitermList res_GBasis =0;
#ifdef ANNA_TIME
  clock_t t1, t2;

  get_cputime(&t1);
  printf("-->> TestSet Gens %d\n", BListLen(theBList));
#endif /* ANNA_TIME */

  Print("--T--  TestSet\n");
  InsertedPairs = 0;
  CoprimePairs = 0;
  GMPairs = 0;
  NonSatPairs = 0;

  INDICES = NULL;

  BLOrder(theBList);

  /* non serve */  ReducedGBElems = BLNew(16, BListIndetsNo(theBList));

  res_GBasis = SaturatedBGrobner(theBList, trunc);

  /*  FREE  */
  /* non serve */  EraseAndFreeBList (ReducedGBElems);
  EraseAndFreeBList (theBList);

  theBList = res_GBasis;

#ifdef STAT
  fprintf(stderr, "--T-- NonSatPairs :  %d\n", NonSatPairs);
  fprintf(stderr, "--T-- CoprimePairs:  %d\n", CoprimePairs);
  fprintf(stderr, "--T-- GMPairs:       %d\n", GMPairs);
  fprintf(stderr, "--T-- InsertedPairs: %d\n", InsertedPairs);
#endif

#ifdef ANNA_TIME
  get_cputime(&t2);
  printf("-->> TestSet time = %7.2f\n", (double)(cputime_diff(t2,t1))/1000);
#endif /* ANNA_TIME */

  return res_GBasis;
}

/******************  START  ******************/

void toric_init()
{
  static int initialized = 0;

    ToricAlg = SATI;
  /*
    ToricAlg = EATI;
    #endif
    #ifdef SEQUENTIAL
    ToricAlg = SATI;
    #endif
    #ifdef HDRIVEN
    ToricAlg = PATI;
    #endif
    */
  if ( !initialized )
  {
    initialized = 1;

    //    poincare_init(0);
    StartPoincare(0);

    //------- disable RUM ----------------------------------
    rum_init(sizeof(BPairStruct), 5000);
    rum_init(sizeof(bitermStruct), 2000);
    //------- disable RUM ----------------------------------
  }
}

void StartTestSet(int* C, int IndetsNo)
{
  rum_init(sizeof(iv_elem)*(IndetsNo+1),5000);  /* ivec */
  rum_init(eterm_size(IndetsNo), RUM_STD_SIZE);

  VECT_WHICH_GT = ( & VectTestSetWhichGT );
  COST = C;
}

void StartTestSetXel(int IndetsNo)
{
  rum_init(sizeof(iv_elem)*(IndetsNo+1),5000);  /* ivec */
  rum_init(eterm_size(IndetsNo), RUM_STD_SIZE);

  VECT_WHICH_GT = ( & VectTestSetXelWhichGT );
}

void StartTestSetLex(int IndetsNo)
{
  rum_init(sizeof(iv_elem)*(IndetsNo+1),5000);  /* ivec */
  rum_init(eterm_size(IndetsNo), RUM_STD_SIZE);

  VECT_WHICH_GT = ( & VectTestSetLexWhichGT );
}

void StartToric(int IndetsNo)
{
  int n = IndetsNo;
  toric_init();

  if ( ToricAlg == EATI )
    n++;

  //------- disable RUM ----------------------------------
  //   rum_init(sizeof(iv_elem)*(n+1),2000);  /* ivec */
  //   rum_init(eterm_size(n), RUM_STD_SIZE);
  //------- disable RUM ----------------------------------

  VECT_WHICH_GT = ( & VectWhichGT );
}

#ifndef ANNA
/******************  NONNEG.C  ******************/

void nonneg1(int *used, int nrows, int ncols, int **M, int *rhs)
{
  int i, j, sum, col;

  for (col=0; col < ncols && used[col] >= 0; col++) {};
  if (col == ncols) return;
  for (i = 0;; i++)
  {
    used[col] = i;
    sum = 0;
    for (j=0; j < nrows; j++) sum += rhs[j];
    if (sum == 0)
    {
      printf("soln: ");
      for (j=0; j < ncols; j++) printf("%3d ", used[j]);
      printf("\n");
    }
    nonneg1(used, nrows, ncols, M, rhs);
    for (j=0; j < nrows; j++) if (M[j][col] > rhs[j]) goto unscramble;

    for (j=0; j < nrows; j++) rhs[j] -= M[j][col];
  }
  unscramble:
  for (j=0; j < nrows; j++) rhs[j] += i*M[j][col];
  used[col] = -1;
}

int* nonneg2(int *used, int nrows, int ncols, int **M, int *rhs)
{
  int i, j, sum, col;
  int* sol;

  for ( col=0 ; col < ncols && used[col] >= 0; col++) {};
  if (col == ncols) return 0;
  for (i = 0;; i++)
  {
    used[col] = i;
    sum = 0;
    for (j=0; j < nrows; j++) sum += rhs[j];
    if (sum == 0)
    {
      //      sol = (int*)MALLOC(ncols*sizeof(int));
      sol = (int*)malloc(ncols*sizeof(int));
      for (j=0; j < ncols; j++) sol[j] = MAX(used[j], 0);
      return sol;
    }
    sol = nonneg2(used, nrows, ncols, M, rhs);
    if (sol != 0) return sol;
    for (j=0; j < nrows; j++) if (M[j][col] > rhs[j]) goto unscramble;
    for (j=0; j < nrows; j++) rhs[j] -= M[j][col];
  }
  unscramble:
  for (j=0; j < nrows; j++) rhs[j] += i*M[j][col];
  used[col] = -1;
  return 0;
}

int* nonneg(int nrows, int ncols, int **M, int* rhs)
{
  int i;
  int *used, *ans;

  used = (int*)malloc(ncols*sizeof(int));
  for (i=0; i < ncols; i++) used[i] = -1;
  /* The loop below eliminates zero columns from consideration */
  for (i=0; i < ncols; i++)
  {
    int j;
    for (j=0; j < nrows; ++j) if (M[j][i] != 0) break;
    if (j == nrows) used[i] = 0;
  }
  ans = nonneg2(used, nrows, ncols, M, rhs);
  free(used);
  return ans;
}

// object ip_solve(const omat M, const olist L)
// {
//   olist answer;
//   int i, j, nrows, ncols;
//   int **matrix, *soln, *vector;
//   int tmp;

//   nrows = omat_nrow(M);
//   ncols = omat_ncol(M);
//   if (olist_len(L) != nrows)
//     return error2object(my_strdup("NonNegSol: incompatible matrix/vector dimensions"));

//   /* Do a sanity check on the input matrix */
//   for (i=0; i < nrows; i++)
//     for (j=0; j < ncols; j++)
//     {
//       int junk = object_cast_int(object_link(omat_ij(M, i+1, j+1)));
//       if (math_error != object_NULL) return math_error;
//       if (junk < 0)
//         return error2object(my_strdup("NonNegSol: matrix entries must be non-negative integers"));
//     }

//   /* Do a sanity check on the input vector */
//   for (i=0; i < nrows; i++)
//   {
//     int junk = object_cast_int(object_link(olist_nth(L, i+1)));
//     if (math_error != object_NULL) return math_error;
//   }

//   /* OK, input is sane, so proceed with the computation */
//   matrix = (int**)MALLOC(nrows*sizeof(int*));
//   for (i=0; i < nrows; i++)
//   {
//     matrix[i] = (int*)MALLOC(ncols*sizeof(int));
//     for (j=0; j < ncols; j++)
//       matrix[i][j] = object_cast_int(object_link(omat_ij(M, i+1, j+1)));
//   }
//   vector = (int*)MALLOC(nrows*sizeof(int*));
//   for (i=0; i < nrows; i++)
//   {
//     tmp = object_cast_int(object_link(olist_nth(L, i+1)));
//     vector[i] = tmp;
//   }
//   soln = nonneg(nrows, ncols, matrix, vector);
//   if ( soln == NULL ) answer = olist_init(0);
//   else
//   {
//     answer = olist_init(ncols);
//     for (i=0; i < ncols; i++)
//       olist_set_nth(answer, i+1, int2object(soln[i]));
//   }

//   for (i=0; i < nrows; i++)  free(matrix[i]);
//   free(vector); free(matrix);
//   free(soln);

//   return olist2object(answer);
// }

#if 0
int main()
{
  int i, j, nrows, ncols;
  int **M, *rhs;

  printf("nrows= ");
  scanf("%d", &nrows);
  printf("ncols= ");
  scanf("%d", &ncols);
  printf("Enter matrix, row by row.\n");
  M = (int**)malloc(nrows*sizeof(int*));
  for (i=0; i < nrows; i++)
  {
    M[i] = (int*)malloc(ncols*sizeof(int));
    for (j=0; j < ncols; j++)
      scanf("%d", &M[i][j]);
  }

  printf("Enter RHS.\n");
  rhs = (int*)malloc(nrows*sizeof(int));
  for (i=0; i < nrows; i++)
    scanf("%d", &rhs[i]);
  nonneg(nrows, ncols, M, rhs);

  return 0;
}
#endif

#endif  /* ANNA */
