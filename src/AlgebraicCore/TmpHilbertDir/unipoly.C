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

#include "gmp.h"
#include "unipoly.h"

/*********    *********/

unipoly *PowerList;

/*********    *********/

unipoly unipoly_dup(unipoly p)
{
  unipoly res_p;
  int d, size;

  if ( (size=UPSize(p)) <= UNIPOLY_RUM.MaxDeg ) size = UNIPOLY_RUM.MaxDeg;
  res_p = unipoly_malloc(size);
  UPSetSize(res_p, size);
  UPSetMax(res_p, UPMax(p));
  UPSetDeg(res_p, UPDeg(p));
  if ( UPMax(p)== 0 )
    for( d=UPSize(p); d>=0 ; --d )  mpz_init_set(res_p[d].z, p[d].z);
  else
    for( d=UPSize(p); d>=0 ; --d )  res_p[d] = p[d];

  return res_p;
}
 
unipoly unipoly_ext_dup(unipoly p, int size)
{
  unipoly res_p;
  int d;
 
  if ( size<=UNIPOLY_RUM.MaxDeg ) size = UNIPOLY_RUM.MaxDeg;
  res_p = unipoly_malloc(size);
  UPSetSize(res_p, size);
  UPSetMax(res_p, UPMax(p));
  UPSetDeg(res_p, UPDeg(p));
  if ( UPMax(p)== 0 )
  {
    for( d=size ; d>UPDeg(p) ; --d ) mpz_init(res_p[d].z);
    for( d=UPDeg(p) ; d>=0 ; --d )   mpz_init_set(res_p[d].z, p[d].z);
  }
  else
  {
    for( d=size; d>UPDeg(p) ; --d ) res_p[d].i = 0;
    for( d=UPDeg(p); d>=0 ; --d )   res_p[d] = p[d];
  }
  
  return res_p;
}

void unipoly_copy(unipoly From, unipoly To)
{
  int d;
 
  if ( UPMax(From)== 0 )
  {
    for( d=UPDeg(From); d>=0 ; --d )  mpz_init_set(To[d].z, From[d].z);
    for( d=-1 ; d>=UPFirst ; --d )    To[d] = From[d];
  }
  else
    for( d=UPDeg(From); d>=UPFirst ; --d )  To[d] = From[d];
}

unipoly NewUnipoly (int size)
{
  unipoly res_p;
  int d;
  
  if ( size<=UNIPOLY_RUM.MaxDeg ) size = UNIPOLY_RUM.MaxDeg;
  res_p = unipoly_malloc(size);
  UPSetSize(res_p,size);
  UPSetMax(res_p, 1);
  UPSetDeg(res_p, 0);
  for ( d=size ; d>=0 ; --d)  res_p[d].i = 0;
  
  return res_p;
}

unipoly NewUnipoly0 (int size)
{
  unipoly res_p;
  int d;
  
  if ( size<=UNIPOLY_RUM.MaxDeg ) size = UNIPOLY_RUM.MaxDeg;
  res_p = unipoly_malloc(size);
  UPSetSize(res_p,size);
  UPSetMax(res_p, 0);
  UPSetDeg(res_p, 0);
  for ( d=size ; d>=0 ; --d)  mpz_init(res_p[d].z);
  
  return res_p;
}

unipoly UnipolyOne (int size)
{
  unipoly res_p;
  int d;
  
  if ( size<=UNIPOLY_RUM.MaxDeg ) size = UNIPOLY_RUM.MaxDeg;
  res_p = unipoly_malloc(size);
  UPSetSize(res_p,size);
  UPSetMax(res_p, 1);
  UPSetDeg(res_p, 0);
  for ( d=size ; d>0 ; --d )  res_p[d].i = 0;
  res_p[0].i = 1;
 
  return res_p;
}

void FreeUnipoly (unipoly p)
{
  int d;
 
  if ( UPMax(p)== 0 )
    for ( d = UPSize(p) ; d>=0 ; --d ) if ( p[d].z!=NULL ) mpz_clear(p[d].z);
  unipoly_free(p);
}

void FreeUnipolyCheck (unipoly p)
{
  int d;

  if ( UPMax(p)== 0 )
    for ( d = UPSize(p) ; d>=0 ; --d ) if (p[d].z!=NULL) mpz_clear(p[d].z);
  unipoly_free(p);
} 

unipoly UnipolyChangeSize(unipoly p, int NewSize)
{  
  int d = UPSize(p);
  
  p = unipoly_realloc(p, NewSize);
  if ( p==(unipoly)UPFirst ) printf("Err UnipolyChangeSize, out of memory\n");
  UPSetSize(p, NewSize);
  if ( UPMax(p)== 0 )  for ( ++d ; d<=NewSize ; ++d) mpz_init(p[d].z);
  else                 for ( ++d ; d<=NewSize ; ++d) p[d].i = 0;

  return p;
}

unipoly UnipolyConvertCoeffs(unipoly p)
{  
  int d = UPSize(p);
  
  for ( ; d>UPDeg(p) ; --d ) mpz_init(p[d].z); // converts the 0s
  for ( ; d>=0 ; --d )       mpz_init_set_si(p[d].z, (long)p[d].i); // cast to long from int
  UPSetMax(p, 0);

  return p;
}

int UnipolyComputeMax(unipoly p)
{  /* assumed int coeffs */  
  int d = UPDeg(p), M = 0, m = 0;
  
  for ( ; d>=0 ; --d ) { if (p[d].i>M) M=p[d].i; if (p[d].i<m) m=p[d].i; }

  return MAX(M, -m);
}

unipoly UnipolyResetMax(unipoly p)
{  /* assumed int coeffs */  
  int d = UPDeg(p), M = 0, m = 0;
  
  for ( ; d>=0 ; --d ) { if (p[d].i>M) M=p[d].i; if (p[d].i<m) m=p[d].i; }
  UPSetMax(p, MAX(M, -m));
  if ( UPMax(p)>INT_MAX/2 )
  {
    for ( d=UPSize(p) ; d>UPDeg(p) ; --d ) mpz_init(p[d].z);
    for ( ; d>=0 ; --d )  mpz_init_set_si(p[d].z, (long)p[d].i); // cast to long from int
    UPSetMax(p, 0);
  }

  return p;
}

/********************  + *  ********************/

void MultByOneMinusXExp (unipoly p, int exp)
{
  int    d =UPDeg(p), de;

  de = d+exp;
  if ( de>UPSize(p) ) 
  { printf("Err:MultByOneMinusXExp\n");p=UnipolyChangeSize(p,de);}
  if ( UPMax(p)!=0 )
  {
    if (d >= exp) for ( ; de>UPDeg(p) ; --de,--d )  p[de].i = -p[d].i;
    for ( ; d>=0 ; --de,--d)   p[de].i = p[de].i - p[d].i;
    UPSetDeg(p, UPDeg(p)+exp);
    UPSetMax(p, 2*UPMax(p));
    if ( UPMax(p)>INT_MAX/2 ) UnipolyResetMax(p);
  }
  else
  {
    if (d >= exp)
      for ( ; de>UPDeg(p) ; --de,--d )	mpz_neg(p[de].z, p[d].z);
    for ( ; d>=0 ; --de,--d) mpz_sub(p[de].z, p[de].z, p[d].z);
    UPSetDeg(p, UPDeg(p)+exp);
  }
}

unipoly P1PlusXExpP2 (unipoly p1, int exp, unipoly p2)
{
  int deg1 =UPDeg(p1), d1, d2 =UPDeg(p2);
  
  d1 = d2+exp;

  if ( UPMax(p1)!=0 && UPMax(p2)!=0 )
  {
    if (deg1 >= d1)      /*    ///XXX///    */
    {
      for ( ; d2>=0 ; --d1,--d2)   p1[d1].i = p1[d1].i + p2[d2].i;
      if ( p1[deg1].i==0 && deg1!=0 )  
      {
	for ( --deg1 ; p1[deg1].i==0 ; --deg1) {}
	UPSetDeg(p1, deg1);
      }
      UPSetMax(p1, UPMax(p1)+UPMax(p2));
      if ( UPMax(p1)>INT_MAX/2 ) UnipolyResetMax(p1);
      unipoly_free(p2); 
      return p1;
    }
    
    if ( UPSize(p1) < d1 )  p1 = UnipolyChangeSize(p1, d1);
    UPSetDeg (p1, d1);
    if (deg1 < exp)      /*    ////  \\\\   */
    {
      for ( ; d2>=0 ; --d1,--d2)  p1[d1] = p2[d2];
      UPSetMax(p1, MAX(UPMax(p1),UPMax(p2)));
      unipoly_free(p2); 
      return p1;
    }
                         /*    ////XX\\\\   */
    for ( ; d1>deg1 ; --d1,--d2)  p1[d1] = p2[d2];
    for ( ; d2>=0 ; --d1,--d2)    p1[d1].i = p1[d1].i + p2[d2].i; 
    UPSetMax(p1, UPMax(p1)+UPMax(p2));
    if ( UPMax(p1)>INT_MAX/2 ) UnipolyResetMax(p1);
    unipoly_free(p2); 
  }
  else
  {
    if ( UPMax(p1)!=0 ) UnipolyConvertCoeffs(p1);
    if ( UPMax(p2)!=0 ) UnipolyConvertCoeffs(p2);
    if (deg1 >= d1)      /*    ///XXX///    */
    {
      for ( ; d2>=0 ; --d1,--d2)  mpz_add(p1[d1].z, p1[d1].z, p2[d2].z);
      if ( mpz_sgn(p1[deg1].z)==0 && ( deg1!=0 ) )  
      {
	for ( --deg1 ; mpz_sgn(p1[deg1].z)==0 ; --deg1) {}
	UPSetDeg (p1, deg1);
      }
      FreeUnipoly(p2); 
      return p1;
    }
    
    if ( UPSize(p1) < d1 )   p1 = UnipolyChangeSize(p1, d1);
    UPSetDeg (p1, d1);
    if (deg1 < exp)      /*    ////  \\\\   */
    {
      for ( ; d2>=0 ; --d1,--d2) mpz_set(p1[d1].z, p2[d2].z);
      FreeUnipolyCheck (p2); 
      return p1;
    }
                         /*    ////XX\\\\   */
    for ( ; d1>deg1 ; --d1,--d2)  mpz_set(p1[d1].z, p2[d2].z);
    for ( ; d2>=0 ; --d1,--d2)    mpz_add(p1[d1].z, p1[d1].z, p2[d2].z);
    FreeUnipolyCheck (p2); 
  }

  return p1;
}

unipoly P1MinusXExpP2 (unipoly p1, int exp, unipoly p2)
{
  int deg1 =UPDeg(p1), d1, d2 =UPDeg(p2);
  
  d1 =d2+exp;
  if ( UPMax(p1)!=0 && UPMax(p2)!=0 )
  {
    if (deg1 >= d1)      /*    ///XXX///    */
    {
      for ( ; d2>=0 ; --d1,--d2)   p1[d1].i = p1[d1].i - p2[d2].i;
      if ( p1[deg1].i==0 && deg1!=0 )  
      {
	for ( --deg1 ; p1[deg1].i==0 ; --deg1) {}
	UPSetDeg(p1, deg1);
      }
      UPSetMax(p1, UPMax(p1)+UPMax(p2));
      if ( UPMax(p1)>INT_MAX/2 ) UnipolyResetMax(p1);
      unipoly_free(p2); 
      return p1;
    }
    
    if ( UPSize(p1) < d1 )   p1 = UnipolyChangeSize(p1, d1);
    UPSetDeg (p1, d1);
    if (deg1 < exp)      /*    ////  \\\\   */
    {
      for ( ; d2>=0 ; --d1,--d2) p1[d1].i = -p2[d2].i;
      UPSetMax(p1, MAX(UPMax(p1),UPMax(p2)));
      unipoly_free(p2); 
      return p1;
    }
                         /*    ////XX\\\\   */
    for ( ; d1>deg1 ; --d1,--d2) p1[d1].i = -p2[d2].i;
    for ( ; d2>=0 ; --d1,--d2)   p1[d1].i = p1[d1].i - p2[d2].i; 
    UPSetMax(p1, UPMax(p1)+UPMax(p2));
    if ( UPMax(p1)>INT_MAX/2 ) UnipolyResetMax(p1);
    unipoly_free(p2); 
  }
  else
  {
    if ( UPMax(p1)!=0 ) UnipolyConvertCoeffs(p1);
    if ( UPMax(p2)!=0 ) UnipolyConvertCoeffs(p2);
    if (deg1 >= d1)      /*    ///XXX///    */
    {
      for ( ; d2>=0 ; --d1,--d2) mpz_sub(p1[d1].z, p1[d1].z, p2[d2].z);
      if ( mpz_sgn(p1[deg1].z)==0 )  
      {
	for ( --deg1 ; mpz_sgn(p1[deg1].z)==0 ; --deg1 ) {}
	UPSetDeg (p1, deg1);
      }
      FreeUnipoly(p2); 
      return p1;
    }
    
    if ( UPSize(p1) < d1 )   p1 = UnipolyChangeSize(p1, d1);
    UPSetDeg (p1, d1);
    if (deg1 < exp)      /*    ////  \\\\   */
    {
      for ( ; d2>=0 ; --d1,--d2) mpz_neg(p1[d1].z, p2[d2].z);
      FreeUnipolyCheck (p2); 
      return p1;
    }
                         /*    ////XX\\\\   */
    for ( ; d1>deg1 ; --d1,--d2) mpz_neg(p1[d1].z, p2[d2].z);
    for ( ; d2>=0 ; --d1,--d2)   mpz_sub(p1[d1].z, p1[d1].z, p2[d2].z);
    FreeUnipolyCheck (p2); 
  }
  return p1;
}

unipoly P1TimesP2 (unipoly p1, unipoly p2)
{
  unipoly res;
  int d2 =UPDeg(p2), d1 =UPDeg(p1), j;
  UPCoeff p2d2;
  double bound = ((double)UPMax(p1))*((double)UPMax(p2));

  res = NewUnipoly(d2+d1);
  UPSetDeg (res, d2+d1);
  //  if ( bound!=0 && bound*MIN(UPDeg(p1),UPDeg(p2))<INT_MAX ) // 2006
  if ( bound!=0 && bound*MIN(d1,d2)<INT_MAX )
  {
    for ( ; d2>=0 ; --d2)   
      if ( p2[d2].i!=0 )
      {
	p2d2 = p2[d2];
	for ( j=d1 ; j>=0 ; --j)  res[d2+j].i += p2d2.i*p1[j].i;
      }
    UPSetMax(res, (int)bound);
    if ( UPMax(res)>INT_MAX/2 ) UnipolyResetMax(res);
    unipoly_free(p1);
    unipoly_free(p2);
  }
  else
  {
    mpz_t z;
    int d;
    mpz_init(z);
    for ( d=UPSize(res) ; d>=0 ; --d ) mpz_init(res[d].z);
    if ( UPMax(p1)!=0 ) UnipolyConvertCoeffs(p1);
    if ( UPMax(p2)!=0 ) UnipolyConvertCoeffs(p2);
    for ( ; d2>=0 ; --d2)   
      if ( mpz_sgn(p2[d2].z)!=0 )
      {
	p2d2 = p2[d2];    
	for ( j=d1; j>=0 ; --j) 
	{
	  mpz_mul(z, p2d2.z, p1[j].z);
	  mpz_add(res[d2+j].z, res[d2+j].z, z);  
	}
      }
    UPSetMax(res, 0);
    FreeUnipoly(p1);
    FreeUnipoly(p2);
    mpz_clear(z);
  }

  return res;
}


unipoly* MakePowerList(int MaxDeg)
{
  int d =0;
  unipoly aux_p =UnipolyOne(MaxDeg), *resPowerList;

  resPowerList = (unipoly*)mycalloc((MaxDeg+1),sizeof(unipoly),"*unipoly");
  
  for ( ; d<MaxDeg ; ++d )
  {
    resPowerList[d] = unipoly_dup(aux_p);
    MultByOneMinusXExp(aux_p, 1);
  }
  resPowerList[MaxDeg] = aux_p;

  return resPowerList;  
}

unipoly PowerExtDup (int d, int size)
{
  int i;
  unipoly res;
  
  if ( d <= PoincareMaxPower ) return unipoly_ext_dup(PowerList[d], size);
  else
  {
    res = unipoly_ext_dup(PowerList[PoincareMaxPower], size);
    for ( i=d ; i>PoincareMaxPower ; --i )  MultByOneMinusXExp(res, 1);
  }
  
  return res;
}
 
/*********  UNIPOLY-RUM  *********/

#define UNIPOLY_NO 150
#define UNIPOLY_RUM_RELOAD 100

unipoly_rum_stack UNIPOLY_RUM;

void unipoly_rum_init(int n)
{
  UNIPOLY_RUM.MaxDeg = n;
  UNIPOLY_RUM.slots = (unipoly*)malloc((UNIPOLY_NO+1)*sizeof(unipoly));
  UNIPOLY_RUM.top   = 0;
}

void unipoly_rum_reset(int n)
{
  if ( n!=UNIPOLY_RUM.MaxDeg )
  {
    int i = UNIPOLY_RUM.top;
    for ( ; i>0 ; --i )
      myfree(UPRealSize(UNIPOLY_RUM.MaxDeg), UNIPOLY_RUM.slots[i] ,"UP");
    UNIPOLY_RUM.MaxDeg = n;
    UNIPOLY_RUM.top   = 0;
  }
}

unipoly unipoly_rum_malloc()
{
  if ( UNIPOLY_RUM.top == 0 ) 
  { 
    int sz = UPRealSize(UNIPOLY_RUM.MaxDeg), i = UNIPOLY_RUM_RELOAD;
    for  ( ; i>0 ; --i ) unipoly_rum_free((unipoly)malloc(sz));
  }
  return UNIPOLY_RUM.slots[UNIPOLY_RUM.top--];
}

void unipoly_rum_free(unipoly p)
{
  if ( UNIPOLY_RUM.top == UNIPOLY_NO )  free(p);
  else UNIPOLY_RUM.slots[++UNIPOLY_RUM.top] = p;
}

int UnipolySubCoeffsOfXd(unipoly p1, unipoly p2, int d)
{
  int i;
  mpz_t z1, z2;

  if ( UPDeg(p1)>=d )
    if ( UPMax(p1)!=0 )  mpz_init_set_si(z1, p1[d].i);
    else mpz_init_set(z1, p1[d].z);
  else mpz_init(z1);
  
  if ( UPDeg(p2)>=d )
  {
    if ( UPMax(p2)!=0 )  mpz_init_set_si(z2, p2[d].i);
    else mpz_init_set(z2, p2[d].z);
    mpz_sub(z1, z1, z2);
    mpz_clear(z2);
  }
  i = mpz_get_si(z1);
  mpz_clear(z1);
  return i;
}

int UnipolyAddCoeffsOfXd(unipoly p1, unipoly p2, int d)
{
  int i;
  mpz_t z1, z2;

  if ( UPDeg(p1)>=d )
    if ( UPMax(p1)!=0 )  mpz_init_set_si(z1, p1[d].i);
    else mpz_init_set(z1, p1[d].z);
  else mpz_init(z1);
  
  if ( UPDeg(p2)>=d )
  {
    if ( UPMax(p2)!=0 )  mpz_init_set_si(z2, p2[d].i);
    else mpz_init_set(z2, p2[d].z);
    mpz_add(z1, z1, z2);
    mpz_clear(z2);
  }
  i = mpz_get_si(z1);
  mpz_clear(z1);
  return i;
}


void FreeUnipolyList (unipoly * UL, int n)
{
  int d = n;

  for ( ; d>=0 ; --d ) FreeUnipoly(UL[d]);
  myfree( (n+1)*sizeof(unipoly), UL, "*unipoly");  
}

/*********  ifdef ANNA  *********/

#ifdef ANNA

void unipoly_sum_of_coeffs(UPCoeff *sum, unipoly p)
{
  int d;
  
  if ( UPMax(p)== 0 )  
  {
    mpz_set_si((*sum).z, 0);
    for (d=UPDeg(p) ; d>=0 ; --d) mpz_add((*sum).z, (*sum).z, p[d].z);
  }
  else
  {
    (*sum).i = 0;
    for (d=UPDeg(p) ; d>=0 ; --d) (*sum).i += p[d].i;
  }  
}

void unipoly_divide_by1minus_x(unipoly p)
{
  int d;

  if ( UPMax(p)== 0 )
    for ( d=1 ; d<UPDeg(p) ; ++d)  mpz_add(p[d].z, p[d].z, p[d-1].z);
  else
    for ( d=1 ; d<UPDeg(p) ; ++d)  p[d].i += p[d-1].i;
  UPSetDeg(p,UPDeg(p)-1);
}

#endif



/*********  new interface functions  *********/

int HasMPZCoeffs(unipoly p)
{
  return ( UPMax(p)== 0 );
}


/*********  Unipoly.c: END  *********/

