#ifndef HilbAnnaUTILS_H
#define HilbAnnaUTILS_H

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

//#include <limits.h>
//#include <stdio.h>
#include <gmp.h>

/*  ============  auxtypes  ============  */

/* Miscellaneous OPTIONS */

#define NOT_COMPUTED -2
#define ANSWER_EQ 0
#define ANSWER_BOTH 0
#define ANSWER_FIRST 1
#define ANSWER_SECOND 2
#define ANSWER_NC 3

/* -----------[ utils ]-------------- */

typedef char coc_bool;

#define FALSE 0
#define TRUE 1

/* -----------[ ivec ]-------------- */
 
typedef int iv_elem;

typedef iv_elem ivec_elem;
typedef ivec_elem * ivec;
 
#define ivec_NULL NULL
 
/* -----------[ sbitsets ]-------------- */

typedef unsigned int bs_subset;
typedef bs_subset shortbitset;

/*  ============  coctypes  ============  */

/* -----------[ eterms ]-------------- */

typedef  int  * eterm;

/* -----------[ unipoly ]-------------- */
typedef union
{
  int i;
  mpz_t z;
} UPCoeff;

typedef UPCoeff * unipoly;

/*  ============  bgrobner  ============  */

/* -----------[ ints ]-------------- */

typedef int * ints;

/* -----------[ generic arrays ]-------------- */

typedef void ** flexarray;

/*  ============  general macros  ============  */

#define MAX(A, B)    ( (A) > (B) ? (A) : (B) )
#define MIN(A, B)    ( (A) > (B) ? (B) : (A) )

/*  ==========  memory management  ==========  */

#include <stdlib.h>
  
/* #ifdef COCOA_WITH_RUM */

#ifndef RUM_H
  #include "rum.h"
#endif /* endif RUM_H */
#define mymalloc(S,D) (((S)<=RUM_N_STACKS) ? rum_malloc(S) : malloc(S))
#define myfree(S,P,D) (((S)<=RUM_N_STACKS) ? rum_free(S,P) : free(P))
#define mycalloc(S,D,M) rum_calloc(S, D, M)
#define myrealloc(P,O,N,M) rum_realloc(P,O,N,M)

/* #else /\* ifndef COCOA_WITH_RUM *\/ */
/* #ifdef JDEBUG */

/* #ifndef JALLOC_H */
/*   #include "MEMORY/jalloc.h" */
/* #endif /\* endif JALLOC_H *\/ */
/* #define mymalloc(S,D) jmalloc(S) */
/* #define myfree(S,P,D) jfree(P) */
/* #define mycalloc(S,D,M) jcalloc(S, D) */
/* #define myrealloc(P,O,N,M) jrealloc(P,N) */

/* #else  /\* ifndef JDEBUG *\/ */
/*   #include "HilbTmpRum.h" */
/* #define mymalloc(S,D) malloc(S) */
/* #define myfree(S,P,D) free(P) */
/* #define mycalloc(S,D,M) rum_calloc(S, D, M) */
/* #define myrealloc(P,O,N,M) rum_realloc(P,O,N,M) */

/* #endif /\* endif JDEBUG *\/ */
/* #endif /\* endif COCOA_WITH_RUM *\/ */

#endif  /*  HilbAnnaUTILS_H  */
