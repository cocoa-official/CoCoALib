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

#ifndef ETERMS_H
# include "eterms.h"
#endif

#ifndef IVECTORS_H
# include "IVectors.h"
#endif


/*********  INT FLEXARRAY  *********/

ints ints_init(int n)
{
  ints res = ints_malloc(n);
  IntsSetSize (res, n);
  IntsSetLen (res, 0);
  return res;
}

ints ints_dup(ints Ints)
{
  ints res = ints_malloc(IntsGetSize(Ints));
  /*register*/ int i=IntsGetLen(Ints);

  IntsSetSize (res, IntsGetSize(Ints));
  IntsSetLen (res, i);
  for ( ; i>0 ; --i ) res[i]=Ints[i];

  return res;
}

/*********   *********/

eterm eterm_init(int n)
{
  eterm t = malloc_eterm(n);
  int i = n;
  
  SetSqFr(t, bs_init(n));
  eterm_set_degree(t,0);
  eterm_set_indetsNo(t,n);
  IntsSetLen(Indets(t),0);
  for ( ; i>0 ; --i ) t[i]=0;
  
  return t;
}

eterm eterm_dup(eterm t)
{
  int n = eterm_get_indetsNo(t);
  eterm aux_t = malloc_eterm(n);
  
  for ( n=eterm_last(n) ; n>0 ; --n )  aux_t[n] = t[n];
  SetSqFr(aux_t, SqFr(t));
  eterm_set_degree(aux_t, eterm_degree(t));
  eterm_set_indetsNo(aux_t, eterm_get_indetsNo(t));

  return aux_t;
}

coc_bool sp_BigMult(eterm t1, eterm t2)
{
  int i, *OccInd2 = Indets(t2);

  for ( i = IntsGetLen(OccInd2) ; i>0 ; --i)
    if ( t1[OccInd2[i]] <= t2[OccInd2[i]]) return FALSE;
  return TRUE;
}

coc_bool sp_SmallMult(eterm t1, eterm t2)
{
  int i, *OccInd2 = Indets(t2);

  for ( i = IntsGetLen(OccInd2) ; i>0 ; --i)
    if ( t2[OccInd2[i]] <= t1[OccInd2[i]]) return TRUE;
  return FALSE;
}

coc_bool sp_Mult(eterm t1, eterm t2)
{
  int i, *OccInd2 = Indets(t2);

  for ( i = IntsGetLen(OccInd2) ; i>0 ; --i)
    if ( t1[OccInd2[i]] < t2[OccInd2[i]]) return FALSE;
  return TRUE;
}

coc_bool sp_eq (eterm t1, eterm t2)
{
  int i = eterm_get_indetsNo(t2);

  for ( ; i>0 ; --i)  if ( t1[i] != t2[i] ) return FALSE;
  return TRUE;
}

coc_bool sp_eterm_coprime (eterm t1, eterm t2)
{
  int i, *OccInd2 = Indets(t2);

  for ( i =IntsGetLen(OccInd2) ; i>0 ; --i)
    if ( t1[OccInd2[i]] != 0 ) return FALSE;
  return TRUE;
}

eterm eterm_gcd(eterm t1, eterm t2)
{
  eterm aux_t =eterm_dup(t1);
  int *OccInd = Indets(aux_t), *OccInd1 = Indets(t1);
  size_t i =IntsGetLen(OccInd1), n, NewOccIndNo=0;

  for ( ; i>0 ; --i )
  {
    n =OccInd1[i];
    if ( t2[n] == 0 )
      //      eterm_put_nth (aux_t, n, 0);  // 2006
      eterm_put0_nth (aux_t, n);
    else
    {
      eterm_put_nth(aux_t, n, MIN(t1[n],t2[n]));
      OccInd[++NewOccIndNo] = n;
    }
  }
  IntsSetLen(OccInd, NewOccIndNo);
  return aux_t;
}

eterm eterm_lcm(eterm t1, eterm t2)
{
  ints OccInd2, OccInd;
  eterm aux_t;
  size_t i, n, aux_indNo;

  if ( IntsGetLen(Indets(t1)) > IntsGetLen(Indets(t2)) )
    aux_t = eterm_dup(t1);
  else
  {
    aux_t = eterm_dup(t2);
    t2 = t1;
  } 
  OccInd2 = Indets(t2);
  OccInd  = Indets(aux_t);

  aux_indNo = IntsGetLen(OccInd);
  for ( i = IntsGetLen(OccInd2) ; i>0 ; --i )
  {
    n = OccInd2[i];
    if ( t2[n] > aux_t[n] )
    {
      if ( aux_t[n] == 0 ) OccInd[++aux_indNo] = n;
      eterm_put_nth(aux_t, n, t2[n]);
    }
  }  
  IntsSetLen(OccInd, aux_indNo);

  return aux_t;
}

eterm eterm_colon(eterm t1, eterm t2)
{
  ints OccInd1 = Indets(t1);
  size_t i = IntsGetLen(OccInd1), n, exp1, exp2;
  
  for ( ; i>0 ; --i )
    if ( (exp2 = t2[(n = OccInd1[i])]) != 0 )    
    {
      if ( exp2 < (exp1 = t1[n]) )
	eterm_put_nth(t1, n, exp1-exp2);
      else
      {
        //	eterm_put_nth(t1, n, 0);  // 2006
	eterm_put0_nth(t1, n);
	IntsMoveLastToNth(OccInd1, i);
      }
    }
  return t1;
}

void eterm_mult_and_assign(eterm t1, eterm t2)
{
  ints OccInd2 = Indets(t2), OccInd1 = Indets(t1);
  size_t i =IntsGetLen(OccInd2), n;
  
  for ( ; i>0 ; --i )
    if ( t1[(n =OccInd2[i])] == 0 )
    {
      eterm_put_nth(t1, n, t2[n]);
      IntsPutLast(OccInd1, n);
    }
    else
      eterm_put_nth(t1, n, t1[n] + t2[n]);
}

void eterm_union_and_assign(eterm t1, eterm t2)
{
  ints OccInd2 = Indets(t2), OccInd1 = Indets(t1);
  size_t i, n;
 
  for ( i =IntsGetLen(OccInd2) ; i>0 ; --i )
    if ( t1[OccInd2[i]] == 0 )
    {
      n = OccInd2[i];
      eterm_put_nth (t1, n, 1);
      IntsPutLast(OccInd1, n);
    }
}


/************************  bgrobner  ************************/

void eterm_complete(eterm t)
{
  size_t i =eterm_get_indetsNo(t);
  ints OccInd = Indets(t);

  for ( ; i>0 ; --i )   if ( t[i]!=0 ) IntsPutLast(OccInd, i);
}

eterm special_eterm_colon_dup (eterm t1, eterm t2)
{
  size_t i =eterm_get_indetsNo(t1), sum=0;
  shortbitset res_sqfr =bs_init(i);
  eterm res =malloc_eterm(i);

  eterm_set_indetsNo(res,i);
  IntsSetLen(Indets(res),0);

  for ( ; i>0 ; --i )
    if ( (res[i]=t1[i]-t2[i])<=0 )
      res[i] = 0;
    else
    {
      sum += res[i];
      bs_put_n(res_sqfr, i);
    }
  eterm_set_degree(res, sum);
  SetSqFr (res, res_sqfr);

  return res;
}

eterm ivec_pos2eterm (ivec v)
{
  size_t i = ivec_len(v);  
  eterm res = eterm_init(i);  
  ints res_ind = Indets(res);

  for ( ; i>0 ; --i )
    if ( ivec_nth(v,i) > 0 )
    {      
      eterm_put_nth(res, i, ivec_nth(v,i));
      IntsPutLast(res_ind, i);      
    }
  /*
    else
    eterm_put0_nth (res, i);
    */
   
  return res;
}

eterm ivec_neg2eterm (ivec v)
{
  size_t i=ivec_len(v);  
  eterm res = eterm_init(i);
  ints res_ind = Indets(res);
  
  for ( ; i>0 ; --i )
    if ( ivec_nth(v,i) < 0 )
    {      
      eterm_put_nth (res, i, -ivec_nth(v,i));
      IntsPutLast(res_ind, i);
    }
    else
      eterm_put0_nth (res, i);
 
  return res;
}

#ifdef ANNA_ETERM_RUM

#define ETERMS_NO 2000
#define ETERMS_RUM_RELOAD 1000

static eterm_rum_stack ETERM_RUM;

void anna_eterm_rum_init(int n)
{
  ETERM_RUM.eterm_size = eterm_size(n);
  ETERM_RUM.slots = (eterm*)malloc((ETERMS_NO+1)*sizeof(eterm));
  ETERM_RUM.top   = 0;
}

eterm eterm_rum_malloc(int size)
{
  if ( ETERM_RUM.top == 0 ) 
  {
    int i = ETERMS_RUM_RELOAD;
    for  ( ; i>0 ; --i ) eterm_rum_free((eterm)malloc(size));
  }
  return ETERM_RUM.slots[ETERM_RUM.top--];
}

void eterm_rum_free(eterm p)
{
  if ( ETERM_RUM.top == ETERMS_NO ) free(p); 
  else ETERM_RUM.slots[++ETERM_RUM.top] = p;
}

#endif /* ETERM_RUM */
