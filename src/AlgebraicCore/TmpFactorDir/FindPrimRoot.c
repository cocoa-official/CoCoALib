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

#include "FindPrimRoot.h"
#include "primes.h"

/* No number less than 2^32 can have more than 10 distinct prime factors.*/
static FFelem prime_factor[10];
static FFelem prime_power[10];

static void Factorize32bits(FFelem n)
/* Factors are placed in prime_factor[], and exponents in prime_power[]. */
/* The end of the factors is marked by prime_factor[k] = 1.              */
{
  FFelem p;
  int nfactors = 0;

  for (p = 2;; p = nextprime(p))
  {
    if (p*p > n) break;
    if (n%p != 0) continue;

    prime_factor[nfactors] = p;
    prime_power[nfactors] = 1;
    n /= p;
    while (n%p == 0)
    {
      n /= p;
      prime_power[nfactors]++;
    }
    nfactors++;
  }
  prime_factor[nfactors] = n;
  if (n==1) return;
  prime_power[nfactors] = 1;
  prime_factor[nfactors+1] = 1;
}


static FFelem BinaryPowerModP(FFelem n, FFelem e, FFelem p) /* p<65536 */
{
  if (e==0) return 1;
  if (e==1) return n;
  if ((e&1)==0) return BinaryPowerModP((n*n)%p, e/2, p);
  return (BinaryPowerModP((n*n)%p, e/2, p) * n)%p;
}


static int FindPrimRootTest(FFelem r, FFelem p)
{
  int i, q;

  i = 0;
  for (q=prime_factor[i]; q != 1; q=prime_factor[++i])
    if (BinaryPowerModP(r, p/q, p) == 1) return 0;
  return 1;
}

FFelem FindPrimRoot(FFelem p)
{
  int q;

  if (p == 2) return 1;
  if (p > MAX_PRIME || !isprime(p)) return 0; /* 0 indicates failure */
  Factorize32bits(p-1);
  for (q=2;;q++)
    if (FindPrimRootTest(q, p)) return q;
}
