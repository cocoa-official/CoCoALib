#ifndef UNIPOLY_H
#define UNIPOLY_H

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

#include "AnnaUtils.h"



/**************************  extern  **************************/

#define Unipoly unipoly
#define PoincareMaxPower 100

extern unipoly *PowerList;

/*********  UNIPOLY-RUM  *********/
typedef struct unipoly_rum_aux {
  unipoly * slots;
  int MaxDeg;
  int top;
} unipoly_rum_stack;

extern unipoly_rum_stack UNIPOLY_RUM;

void unipoly_rum_init(int n);
void unipoly_rum_reset(int n);
unipoly unipoly_rum_malloc();
void    unipoly_rum_free(unipoly p);

/*********    ALLOC-FREE   *********/
/* #define unipoly_malloc(sz) \ */
/*    ( ( sz == UNIPOLY_RUM.MaxDeg ) ? \ */
/*      (unipoly_rum_malloc(sz)+3) : (unipoly)mymalloc(UPRealSize(sz),"UP")+3 ) */
#define unipoly_malloc(sz) \
   ( ( sz == UNIPOLY_RUM.MaxDeg ) ? \
     (unipoly_rum_malloc()+3) : (unipoly)mymalloc(UPRealSize(sz),"UP")+3 )
#define unipoly_free(p) \
   { if ( UPSize(p) == UNIPOLY_RUM.MaxDeg ) \
     unipoly_rum_free((p)-3); else myfree(UPRealSize(UPSize(p)),(p)-3,"UP");}
#define unipoly_realloc(p,newsz) \
   (unipoly)myrealloc((p)-3, UPRealSize(UPSize(p)), UPRealSize(newsz), "UP")+3

/**************************    **************************/

#define UPRealSize(sz) ((sz+4)*sizeof(UPCoeff))
#define UPFirst  -3
#define UPLast(n) n

#define UPDeg(P)    (P)[-1].i
#define UPMax(P)    (P)[-2].i
#define UPSize(P)   (P)[-3].i
#define UPSetDeg(P,D)  UPDeg(P)  = (D)
#define UPSetMax(P,M)  UPMax(P)  = (M)
#define UPSetSize(P,S) UPSize(P) = (S)

/* TODO remove */
#define UnipolyDegree(P)      UPDeg(P)
/* TODO remove end */

unipoly unipoly_dup(unipoly p);

unipoly NewUnipoly(int size);
unipoly UnipolyOne(int size);
void FreeUnipoly(unipoly p);

unipoly UnipolyChangeSize(unipoly p, int NewSize);

/*  every function frees all the arguments  */
unipoly P1PlusXExpP2  (unipoly P1, int exp, unipoly P2);
unipoly P1MinusXExpP2 (unipoly P1, int exp, unipoly P2);
unipoly P1TimesP2     (unipoly P1, unipoly P2);
void MultByOneMinusXExp (unipoly p, int exp);

unipoly* MakePowerList (int n);
unipoly PowerExtDup (int n, int size);

int UnipolySubCoeffsOfXd(unipoly p1, unipoly p2, int d);
int UnipolyAddCoeffsOfXd(unipoly p1, unipoly p2, int d);

void FreeUnipolyList (unipoly * UL, int n);

/*********  ifdef ANNA  *********/
#ifdef ANNA
void unipoly_sum_of_coeffs(UPCoeff *sum, unipoly p);
void unipoly_divide_by1minus_x(unipoly p);
#endif


int HasMPZCoeffs(unipoly p);

#endif




