//   Copyright (c)  1997-2006  John Abbott,  Anna M Bigatti

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

/* need stddef for NULL */
#include <stddef.h>

#include "FF.h"
#include "jaaerror.h"
#include "jalloc.h"
#include "primes.h"
#include "FindPrimRoot.h"

/* Global in which the "current" finite field is specified */
struct FFstruct CurrentFF;


FF FFctor(FFelem p)
{
  FF ans;
  FFelem r, i, s, p1;
  FFelem *Exp;
  unsigned int *Log;

  if (p > MAX_PRIME || !isprime(p))
  {
    JERROR(JERROR_FF_BAD_P);
    return (FF)NULL;
  }
  ans = (FF)MALLOC(sizeof(struct FFstruct));
  ans->prime = p;
  ans->LogTable = NULL;
  ans->ExpTable = NULL;
  ans->k = MAX_FFelem/p/p/2;
  ans->shift = p*(MAX_FFelem/p/2);
  if (p == 2) return ans;      /* Don't create log/exp tables for p=2.       */
  if (p > 65536) return ans;   /* Don't create log/exp tables for big primes */

  Log = (unsigned int*)MALLOC(p*sizeof(unsigned int));
  Exp = (FFelem*)MALLOC((2*p-2)*sizeof(FFelem));
  if ((Log == NULL) || (Exp == NULL) || (ans == NULL))
  {
    if (Log) FREE(Log);
    if (Exp) FREE(Exp);
    if (ans) FREE(ans);
    JERROR(JERROR_NO_MEM);
  }
  r = FindPrimRoot(p);
  p1 = p-1; /* constant */

  s = 1;
  for (i=0; i<p1/2; i++)
  {
    Log[s] = i;
    Log[p-s] = i+p1/2;
    Exp[i] = s;
    Exp[i+p1/2] = p-s;
    Exp[i+p1] = s;
    Exp[i+3*p1/2] = p-s;
    s = (s*r)%p;
  }
  ans->LogTable = Log;
  ans->ExpTable = Exp;
  return ans;
}


void FFdtor(FF Fp)
{
  if (Fp->LogTable) FREE(Fp->LogTable);
  if (Fp->ExpTable) FREE(Fp->ExpTable);
  FREE(Fp);
}


void FFselect(const FF Fp)
{
  CurrentFF.prime = Fp->prime;
  CurrentFF.LogTable = Fp->LogTable;
  CurrentFF.ExpTable = Fp->ExpTable;
  CurrentFF.k = Fp->k;
  CurrentFF.shift = Fp->shift;
}


FFelem FFchar(const FF Fp)
{
  return(Fp->prime);
}


FFelem FFadd(const FFelem x, const FFelem y)
{
  FFelem sum = x+y;
  FFelem p = CurrentFF.prime;

  if (sum >= p) return sum-p;
  return sum;
}


FFelem FFsub(const FFelem x, const FFelem y)
{
  FFelem p = CurrentFF.prime;

  if (x < y) return x+(p-y);
  return x-y;
}


FFelem FFmul(const FFelem x, const FFelem y)
{
  unsigned int *Log = CurrentFF.LogTable;
  FFelem *Exp = CurrentFF.ExpTable;

  if (Log == NULL) return (x*y)%CurrentFF.prime;
  if ((x==0) || (y==0)) return 0;
  return Exp[Log[x]+Log[y]];
}


/* These ought to be inline... (not possible in standard C) */
/* Will crash if called when the Log/Exp tables are NULL.   */
unsigned int FFlog(FFelem x)
{
  return CurrentFF.LogTable[x];
}


FFelem FFexp(unsigned int x)
{
  return CurrentFF.ExpTable[x];
}


/* Compute num/den modulo p using extended Euclidean algorithm.        */
/* Note that cofacA is minus what you might reasonably guess it to be. */

static FFelem FFquot(FFelem num, FFelem den, FFelem p)
{
  FFelem cofacA, cofacB;
  FFelem A, B, q;
  
  cofacA = 0;
  cofacB = 1;
  A = p;
  B = den;
  while (B > 1)
  {
    q = A/B;
    A -= q*B;
    cofacA += q*cofacB;
    q = B/A;
    B -= q*A;
    cofacB += q*cofacA;
  }
  if (B) return (num*cofacB)%p;
  return (num*(p-cofacA))%p;
}


FFelem FFdiv(const FFelem x, const FFelem y)
{
  unsigned int *Log = CurrentFF.LogTable;
  FFelem *Exp = CurrentFF.ExpTable;
  int p = CurrentFF.prime;

  if (y==0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return 0;
  }
  if (x==0) return 0;
  if (Log == NULL)
  {
    if (p == 2) return 1;
    else return FFquot(x, y, CurrentFF.prime);
  }
  return Exp[Log[x]-Log[y]+p-1];
}


FFelem FFpower(const FFelem x, int e)
{
  if (e < 0) return FFpower(FFdiv(1, x), -e);
  if (e == 0) return 1;
  if (x == 0 || x == 1 || e == 1) return x;

  if (CurrentFF.LogTable != NULL)
  {
    int p = CurrentFF.prime;
    unsigned int *Log = CurrentFF.LogTable;
    FFelem *Exp = CurrentFF.ExpTable;
    
    return Exp[((e%(p-1))*Log[x])%(p-1)];
  }
  if (e&1) return FFmul(x, FFpower(FFmul(x, x), (e-1)/2));
  return FFpower(FFmul(x, x), e/2);
}
