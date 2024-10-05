#ifndef TERMLIST_H
#define TERMLIST_H

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

#include "AnnaUtils.h"

/**********************************************************/
typedef eterm      * MixedTermList;
typedef struct TermListStruct * TermList;
 
typedef struct TermListStruct
{
  MixedTermList   MixedTerms;
  int MixedTermListSize;
  int MixedTermListLen;
} TermListStruct; 

#define TListIndetsNo(theTList) eterm_get_indetsNo(SPList(theTList))
#define SPList(theTList)      ((theTList)->MixedTerms)[0]

void EraseAndFreeTList(TermList theTList);
void InterreduceTList (TermList theTList);

void MoveNotCoprimeSP(eterm FromSPList, eterm ToSPList, eterm theTerm);
void MoveNotCoprime(TermList FromTList, TermList ToTList, eterm theTerm);
TermList SplitIndets(TermList theTList);
void ReduceAndDivideBySimplePower( TermList theTList, TermList *DivTList,
                                   size_t PIndex, unsigned PExp);
void ReduceAndDivideByPivot(TermList theTList, TermList *DivTList, eterm Pivot);
void BigPivotOf(TermList theTList, int *PivotIndex, int *PivotDeg);
eterm GCD3PivotOf(TermList theTList);
void GlobalTermList_init(int IndetsNo);
void GlobalTermList_free();
int GlobalTermList_size();

/********************    list of terms    ********************/

TermList NewTList(int len, int indetsNo);
void InsInTList(TermList TL, eterm t, int *NewLen);

void TListChangeSize(TermList theTList, int NewLen);
#define TListReduceSize(TL,newsz) \
        {  SetMTListLen(TL,newsz);\
	   if ( MTListSize(TL)!=199 && MTListSize(TL) > newsz*16)\
	     TListChangeSize(TL, newsz);}

/********************    list of mixed terms    ********************/

#define MTList(theTList)              (theTList)->MixedTerms
#define MTListSize(theTList)          (theTList)->MixedTermListSize
#define MTListLen(theTList)           (theTList)->MixedTermListLen
#define SetMTListSize(theTList,size)  (theTList)->MixedTermListSize = size
#define SetMTListLen(theTList,len)    (theTList)->MixedTermListLen = len

#define MTLPutNth(MTList,N,t)           (MTList)[N] = t
#define MTLPutLast(MTList,len,t)        (MTList)[++len] = t
#define MTLMoveLastToNth(MTList,len,N)  (MTList)[N] = (MTList)[len--]

/********************    list of simple powers    ********************/
/*
    SPList(theTList) is represented with an eterm ( MTL[0] )
    which is the product of the simple powers.
*/
void InsInSPList(size_t Index, unsigned Exp, eterm SPL);


/**********************************************************************/

#endif /* TERMLIST_H */
