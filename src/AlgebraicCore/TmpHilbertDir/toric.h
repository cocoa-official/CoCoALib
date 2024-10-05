#ifndef TORIC_H
#define TORIC_H

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
/*-----------------------------------------------------------------------------
	Part of project CoCoA
-----------------------------------------------------------------------------*/

/* $Id: toric.h,v 1.2 2022/02/23 21:44:40 abbott Exp $ */

#ifndef BITERMS_H
#define BITERMS_H

#include "AnnaUtils.h"

/* #ifndef ANNA */

/* #ifndef AUXTYPES_H */
/* # include "auxtypes.h" */
/* #endif */

/* #ifndef COCTYPES_H */
/* # include "coctypes.h" */
/* #endif */

/* #else */

/* #ifndef AnnaUTILS_H */
/* # include "AnnaUtils.h" */
/* #endif */

/* #endif  */

/* -----------[ Algorithms for Toric Ideals ]-------------- */

typedef enum 
{
  EATI,
  PATI,
  SATI
} ToricAlgType;

extern ToricAlgType ToricAlg;

/* -----------[ biterm ]-------------- */

/*
typedef eterm *biterm;
*/

typedef struct bitermStruct  * biterm;

typedef struct bitermStruct
{
  ivec Vect;
  eterm Lt;
  coc_bool bitermStruct_MINIMAL;
}
bitermStruct;

/* -----------[ biterms ]-------------- */

typedef biterm * biterms;

/* -----------[ BList ]-------------- */

typedef struct BitermListStruct * BitermList;
typedef struct BitermListStruct
{
  biterms  BitermListStructBiterms;
  biterms  BitermListStructReducers;
  int      BitermListStructMaxDeg;
  int    * BitermListStructLabels;  
  int    * BitermListStructInfo_Buckets;  
  int    * BitermListStructInfo_BucketLens;  
} BitermListStruct;

/* -----------[ bpair ]-------------- */

typedef struct BPairStruct * bpair;
typedef struct BPairStruct
{
  biterm   BPairStruct1;
  biterm   BPairStruct2;
  int   BPairStructDeg;
  coc_bool  BPairStructLowerDeg;
} BPairStruct;

/* -----------[ bpairs ]-------------- */

typedef   bpair   * bpairs;

/* -----------[ BPList - double array]-------------- */

typedef struct BPairListStruct * BPairList;
typedef struct BPairListStruct
{
  bpairs BPairListStructMinDegBPairs;
  bpairs BPairListStructLowerMinDegBPairs;
  bpairs BPairListStructOtherBPairs;
} BPairListStruct;

/**************************  extern  **************************/

extern ints INDICES;
extern ivec ElimVect;

/*********  BPAIRS  (flexarray)  *********/

#define bpairs_malloc(sz)         (bpairs)flexarray_malloc(sz, "bpairs")
#define bpairs_free(BPs)          flexarray_free(BPs, "bpairs")
#define bpairs_realloc(BPs,newsz) \
                           (bpairs)flexarray_realloc(BPs,newsz,"bpairs")

/*********  BITERMS (flexarray)  *********/

#define biterms_malloc(len)    (biterms)flexarray_malloc(len, "biterms")
#define biterms_free(Bs)       flexarray_free(Bs, "biterms")
#define biterms_realloc(Bs,newsz) \
                          (biterms)flexarray_realloc(Bs,newsz,"biterms")

/*********  ETERMS  (flexarray)  *********/

//----------------------------------------------------------------------
#define EsSetLen(A,L)    SetLen(A,reinterpret_cast<eterm>(L))
#define EsSetSize(A,S)   SetSize(A,reinterpret_cast<eterm>(S))
//----------------------------------------------------------------------

#define eterms_malloc(len)       (eterms)flexarray_malloc(len, "eterms")
#define eterms_free(I)           flexarray_free(I, "eterms")
#define eterms_realloc(I,newsz) \
                             (eterms)flexarray_realloc(I,newsz,"eterms")

/* -----------[ eterms ]-------------- */

typedef eterm * eterms;

/*********  BITERM  *********/

#define biterm_malloc()  (biterm)mymalloc(sizeof(bitermStruct), "bitermSt")
#define biterm_struct_free(t) myfree(sizeof(bitermStruct),t,"bitermSt")
#define biterm_free(t) \
    {if (t->Lt!=NULL) eterm_free(t->Lt);\
     if (t->Vect!=NULL) ivec_free(t->Vect);\
     biterm_struct_free(t); /*debug* t=NULL;*/}

#define BitermIsMinimal(t) t->bitermStruct_MINIMAL
#define BitermSetIsMinimal(t,s) BitermIsMinimal(t)=s

#define biterm_get_indetsNo(b)  eterm_get_indetsNo((b)->Lt)
#define BitermLtDeg(b)          eterm_degree((b)->Lt)

#define BitermTl(b) ivec_neg2eterm(b->Vect)

void   BitermSet(biterm B, eterm Lt, ivec Vect);
biterm BitermMake(ivec Vect, ivec_elem Index);
biterm biterm_dup(biterm t);
biterm BitermReduce(biterm B, biterm BReducer, ivec_elem Index);

/*********  BPAIR  *********/

#define bpair_malloc()  (bpair)mymalloc(sizeof(BPairStruct), "bpairSt")
#define bpair_free(P)    myfree(sizeof(BPairStruct),P,"bpairSt")

#define BPair1(P)    P->BPairStruct1
#define BPair2(P)    P->BPairStruct2
#define BPair2Index(P)    P->BPairStruct2Index
#define BPairLowerDeg(P)  P->BPairStructLowerDeg
#define BPairDeg(P)  P->BPairStructDeg

#define BPairSet1(P,n)    BPair1(P)=n
#define BPairSet2(P,n)    BPair2(P)=n
#define BPairSet2Index(P,n)    BPair2Index(P)=n
#define BPairSetLowerDeg(P,n)  BPairLowerDeg(P)=n
#define BPairSetDeg(P,n)  BPairDeg(P)=n

/*********  BPAIRS LISTS  *********/

//----------------------------------------------------------------------
#define BPsSetLen(A,L)    SetLen(A,reinterpret_cast<BPairStruct*>(L))
#define BPsSetSize(A,S)   SetSize(A,reinterpret_cast<BPairStruct*>(S))
//----------------------------------------------------------------------

#define BPList_malloc() \
             (BPairList)mymalloc(sizeof(BPairListStruct),"BPairList")
#define BPList_free(p)  myfree(sizeof(BPairListStruct),p,"BPairList")

#define BPLLowerMinDegPairs(BPL) ((BPL)->BPairListStructLowerMinDegBPairs)
#define BPLLowerMinDegSize(BPL)   GetSize(BPLLowerMinDegPairs(BPL))
#define BPLLowerMinDegLen(BPL)    GetLen(BPLLowerMinDegPairs(BPL))

#define BPLMinDegPairs(BPL) ((BPL)->BPairListStructMinDegBPairs)
#define BPLMinDegSize(BPL)   GetSize(BPLMinDegPairs(BPL))
#define BPLMinDegLen(BPL)    GetLen(BPLMinDegPairs(BPL))

#define BPLPairs(BPL)  ((BPL)->BPairListStructOtherBPairs)
#define BPLSize(BPL)   GetSize(BPLPairs(BPL))
#define BPLLen(BPL)    GetLen(BPLPairs(BPL))

#define BPLMinDeg(BPL)  ( (BPLLowerMinDegLen(BPL)!=0) ? \
			  BPairDeg((BPLLowerMinDegPairs(BPL))[1]) :\
			  BPairDeg((BPLMinDegPairs(BPL))[1]) )

#define BPLSetLowerMinDegPairs(BPL,BPs) BPLLowerMinDegPairs(BPL)=BPs
#define BPLSetLowerMinDegSize(BPL,s)    BPsSetSize(BPLLowerMinDegPairs(BPL),s)
#define BPLSetLowerMinDegLen(BPL,l)     BPsSetLen(BPLLowerMinDegPairs(BPL),l)
#define BPLSetMinDegPairs(BPL,BPs) BPLMinDegPairs(BPL)=BPs
#define BPLSetMinDegSize(BPL,s)    BPsSetSize(BPLMinDegPairs(BPL),s)
#define BPLSetMinDegLen(BPL,l)     BPsSetLen(BPLMinDegPairs(BPL),l)
#define BPLSetPairs(BPL,BPs) BPLPairs(BPL)=BPs
#define BPLSetSize(BPL,s)    BPsSetSize(BPLPairs(BPL),s)
#define BPLSetLen(BPL,l)     BPsSetLen(BPLPairs(BPL),l)

#define BPListResize(BPL,newsz)  \
        if (BPLSize(BPL)<newsz) BPLChangeSize(BPL,2*(newsz))

BPairList BPLNew(int Size);
void BPLChangeSize (BPairList theBPList, int NewSize);
void FreeBPList (BPairList theBPList);

/*********  BITERMS LISTS  *********/

//----------------------------------------------------------------------
#define BsSetLen(A,L)    SetLen(A,reinterpret_cast<bitermStruct*>(L))
#define BsSetSize(A,S)   SetSize(A,reinterpret_cast<bitermStruct*>(S))
//----------------------------------------------------------------------

#define BList_malloc() \
        (BitermList)mymalloc(sizeof(BitermListStruct),"BitermList")
#define BList_free(p)  {myfree(sizeof(BitermListStruct),p,"BitermList");}

#define Biterms(BL)      (BL)->BitermListStructBiterms
#define BListSize(BL)    GetSize(Biterms(BL))
#define BListLen(BL)     GetLen(Biterms(BL))
#define BListMaxDeg(BL)  (BL)->BitermListStructMaxDeg
#define Reducers(BL)     (BL)->BitermListStructReducers
#define ReducersLabels(BL)     (BL)->BitermListStructLabels
#define ReducersLabelNth(BL,i) ((BL)->BitermListStructLabels)[i]
#define BListIndetsNo(BL)  eterm_get_indetsNo((Biterms(BL))[1]->Lt)

#define BListSetSize(BL,size)   BsSetSize(Biterms(BL),size)
#define BListSetLen(BL,len)     BsSetLen(Biterms(BL),len)
#define BListSetMaxDeg(BL,deg)  BListMaxDeg(BL) = deg

#define SetBiterms(BL,l)        Biterms(BL) = l
#define SetReducers(BL,l)       Reducers(BL) = l
#define SetReducersLabels(BL,l) ReducersLabels(BL) = (l)
#define SetReducersLabelNth(BL,i,l) ReducersLabelNth(BL,i) = l

#define BListResize(BL,newsz)  \
        if (BListSize(BL)<newsz) BListChangeSize(BL,2*(newsz))

BitermList BLNew (int len, int indetsNo);
void EraseAndFreeBList (BitermList theBList);
void BListChangeSize (BitermList theBList, int NewSize);

ivec SBitermVect (bpair BP);
biterm SBiterm(biterms BL, int n1, int n2, ivec_elem Index);
biterm BitermNormalForm(biterm B, BitermList theBList, ivec_elem Index);
void SetOrderedBiterms (BitermList theBList);
void BLAddBiterm (BitermList theBList, biterm B);

/******************  TORIC  ******************/

ints HostenShapiro(BitermList BL);
BitermList  BGrobner(BitermList theBList, ivec_elem Index);
BitermList  ElimToric(BitermList theBList, ints Indices, ints Weights);
BitermList  SequentialToric(BitermList theBList, ints Indices);
BitermList  HDToric(BitermList theBList, ints Indices);
BitermList* Toric(BitermList theBList, ints Indices);
BitermList  TestSet(BitermList theBList, int trunc);

void toric_init();
/* void StartTestSet(float* C, int IndetsNo); */
void StartTestSet(int* Cost, int IndetsNo);
void StartTestSetLex(int IndetsNo);
void StartTestSetXel(int IndetsNo);
void StartToric(int IndetsNo);

int null_space(int ***basis, int **M, int nrows, int ncols);

/* #ifndef ANNA */

/* BitermList olist2BList(olist L, ring R); */
/* ints ring_weilist2Ints(ring R); */
/* ints ring_2ordrow2Ints(ring R); */
/* ints olist2Ints(olist L, ring R); */
/* olist this_BList2olist(BitermList BL, ints Weights, ring R); */
/* olist this_WBList2olist(BitermList BL, ints Weights, ring R); */
/* void omat_ker2BList(omat M, BitermList *BL, ints *Indices, ints *Weights); */

/* #endif */

/******************  NONNEG.H  ******************/

/* #ifndef ANNA */

/* #include <stdio.h> */
/* #include <stdlib.h> */

/* /\******************\/ */
/* /\* CoCoA includes *\/ */
 
/* #include "coctypes.h" */
/* #include "ring.h" */
/* #include "polys.h" */
/* #include "utils.h" */
/* #include "coeffs.h" */
/* #include "terms.h" */
/* #include "ivectors.h" */
/* #include "options.h" */
/* #include "object.h" */
/* #include "olist.h" */
/* #include "omat.h" */
/* #include "monomial.h" */
/* #include "mathcoc.h" */
/* #include "ck.h" */
/* #include "ca.h" */
/* /\******************\/ */

/* object ip_solve(omat M, olist L); */
/* #endif  /\* ANNA *\/ */

#endif  /* BITERMS_H */

#endif /* TORIC_H */
