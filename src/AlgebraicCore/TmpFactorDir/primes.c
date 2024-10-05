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

#include "primes.h"
#include "jalloc.h"
#include "jaaerror.h"

/* mask[n%30] == 0 iff n is divisible by at least one of 2, 3 and 5.  */
/* n+skip[n%30] is the next number after n not divisible by 2, 3 or 5 */
int mask[30] = {0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,0,1};
int skip[30] = {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2};
int fall[30] = {1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6};

/* Function to compute a^n mod p */
static FFelem pwr(FFelem a, FFelem n, FFelem p)
{
  FFelem b;

  if (n == 0) return 1;
  if (n == 1) return a;
  b = pwr((a*a)%p, n/2, p);
  if ((n&1)) b = (a*b)%p;
  return b;
}


static int spsp(FFelem b, FFelem n)
{
  int i, r;
  FFelem q, n1, power;

  if (n < 2 || b%n == 0) return 0;
  q = n-1;
  r = 0;
  while ((q&1) == 0) { r++; q /= 2; }
  n1 = n-1;
  power = pwr(b, q, n);
  if (power == n1 || power == 1) return 1;
  for (i=0; i < r; i++)
  {
    power = (power*power)%n;
    if (power == n1) return 1;
    if (power == 1) return 0;
  }
  return 0;
}


/* This definition is for 64-bit machines, but remains valid for 32-bitters */
int isprime(FFelem n)
{
  if (n ==  746331041UL ||
      n == 2840871041UL ||
      n == 3014101261UL ||
      n == 3215031751UL)
    return 0;
  if (n == 2 || n == 3 || n == 5 || n == 7) return 1;
  if (n < 11) return 0;
  return spsp(2, n) && spsp(5, n) && spsp(7, n);
}


FFelem nextprime(FFelem n)
{
  int n30, delta;

  if (n < 2) return 2;
  if (n == 2) return 3;
  if (n < 5) return 5;

  n30 = n%30;
  do
  {
    delta = skip[n30];
    if (n30 == 29) n30=1;
    else n30 += delta;
    n += delta;
    if (n > MAX_PRIME) return 0; /* to indicate overflow */
  } while (!isprime(n));
  return n;
}


FFelem prevprime(FFelem n)
{
  int n30, delta;

  if (n < 3) return 0;
  if (n == 3) return 2;
  if (n < 6) return 3;
  if (n < 8) return 5;

  n30 = n%30;
  do
  {
    delta = fall[n30];
    n30 -= delta;
    if (n30 < 0) n30 += 30;
    n -= delta;
  } while (!isprime(n));
  return n;
}


PrimeSource PrimeSourceCtor(void)
{
  PrimeSource ans;

  ans = (PrimeSource)MALLOC(sizeof(struct PrimeSourceStruct));
  ans->start = nextprime(MAX_PRIME/3); /* heuristic/guess */
  ans->p = ans->start;
  return ans;
}

void PrimeSourceDtor(PrimeSource PS)
{
  FREE(PS);
}

FFelem NextPrime(PrimeSource PS)
{
  if (PS->start > 0)
  {
    PS->p = nextprime(PS->p);
    if (PS->p != 0) return PS->p;
    PS->p = PS->start;
    PS->start = 0;
  }
  PS->p = prevprime(PS->p);
  if (PS->p == 0) JERROR(JERROR_PRIMES_EXHAUSTED);
  return PS->p;
}
