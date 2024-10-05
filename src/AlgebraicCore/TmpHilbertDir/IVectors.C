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

/* --------------------- *\
           IVEC
\* --------------------- */

ivec ivec_init(int n)
{
  ivec aux_v;
  aux_v = malloc_ivec(n);
  ivec_set_len(aux_v,n);
  return aux_v;
}

ivec ivec_sum(ivec v1, ivec v2)
{
  ivec res_v;
  /*register*/ int i;
  res_v=malloc_ivec(ivec_len(v1));
  for (i=ivec_len(v1);i>0;--i)
    ivec_set_nth(res_v,i,ivec_nth(v1,i)+ivec_nth(v2,i));
  ivec_set_len(res_v,ivec_len(v1));
  return res_v;
}

ivec ivec_sub(ivec v1, ivec v2)
{
  ivec res_v;
  /*register*/ int i;
  res_v=malloc_ivec(ivec_len(v1));
  for (i=ivec_len(v1);i>0;--i)
    ivec_set_nth(res_v,i,ivec_nth(v1,i)-ivec_nth(v2,i));
  ivec_set_len(res_v,ivec_len(v1));
  return res_v;
}

ivec ivec_dup(ivec v)
{
  ivec res_v;
  /*register*/ int i;
  int n=ivec_len(v);
  res_v=malloc_ivec(ivec_len(v));
  for (i=n;i>0;--i) ivec_set_nth(res_v,i,ivec_nth(v,i));
  ivec_set_len(res_v,n);
  return res_v;
}

// /****  ANNA  ****/

// ivec ivec_this_inv(ivec v)
// {
//   /*register*/ int i;
//   for (i=ivec_len(v);i>0;--i) ivec_set_nth(v,i,-ivec_nth(v,i));
//   return v;
// }

// /* todo: remove */
// coc_bool ivec_neg_coprime(ivec v1, ivec v2)
// {
//   /*register*/ int i;
//   for (i=ivec_len(v1) ; i>0 ; --i ) 
//     if ( ivec_nth(v1,i)<0 && ivec_nth(v2,i)<0 ) return FALSE;
//   return TRUE;
// }
