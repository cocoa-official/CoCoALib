#ifndef IVECTOR_H
#define IVECTOR_H
 
//  Copyright (c)  2011  Anna Bigatti,  John Abbott

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

#ifndef AnnaUTILS_H
#include "AnnaUtils.h"
#endif

/* --------------------- *\
           IVEC
\* --------------------- */
 
#define malloc_ivec(n) (ivec)mymalloc(sizeof(ivec_elem)*(n+1),"ivec")
#define free_ivec(n,p) myfree(sizeof(ivec_elem)*(n+1),p,"ivec")

#define ivec_len(v) (v)[0]
#define ivec_set_len(v,n) (v)[0]=n
#define ivec_len_ck(v) (v==ivec_NULL ? ZERO : (v)[0])
#define ivec_nth(v,n) ((v)[n])
#define ivec_set_nth(v,n,x) (v)[n]=x
#define ivec_free(v) free_ivec(ivec_len(v),v)
#define ivec_sum_ck(v1,v2) (v1 == ivec_NULL ? v1 : ivec_sum(v1,v2))
 
#define ivec_dup_ck(v) (v==ivec_NULL ? ivec_NULL : ivec_dup(v))
#define ivec_dup_to_ck(v,n) (n==ZERO ? ivec_NULL : ivec_dup_to(v,n))
 
#define ivec_is_zero_ck(v) (v==ivec_NULL ? TRUE : ivec_is_zero(v))
 
ivec ivec_init(int n);
ivec ivec_dup(ivec v);
ivec ivec_sum(ivec v1, ivec v2);
coc_bool ivec_is_zero(ivec v);

/* --------------  ANNA  ------------- */
ivec ivec_sub(ivec v1, ivec v2);
/* ivec ivec_this_inv(ivec v); */  // private in toric.C
/* coc_bool ivec_neg_coprime(ivec v1, ivec v2); */

#endif  /* IVECTOR_H */


