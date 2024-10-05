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

#include "mpz_from_padic.h"
#include "jalloc.h"

/* This one uses memory but is pretty quick.  A recursive one might  */
/* use a little less memory and be a little quicker, but this is far */
/* simpler! (and probably has better malloc/free characteristics).   */

void mpz_from_padic(mpz_t ans, FFelem p, int k, FFelem *digits)
{
  int i, n;
  mpz_t *table, q;

  mpz_init_set_ui(q, (unsigned long)p);
  n = k;
  table = (mpz_t*)MALLOC(n*sizeof(mpz_t));
  for(i=0; i < n; i++) mpz_init_set_ui(table[i], (unsigned long)digits[i]);
  while (n > 1)
  {
    for (i=0; i < n-1; i+=2)
    {
      mpz_mul(table[i+1], table[i+1], q);
      mpz_add(table[i/2], table[i], table[i+1]);
    }
    if (n&1) mpz_set(table[n/2], table[n-1]);
    n = (n+1)/2;
    if (n > 1) mpz_mul(q, q, q);
  }
  mpz_set(ans, table[0]);
  for(i=0; i < k; i++) mpz_clear(table[i]);
  FREE(table);
  mpz_clear(q);
}

#if 0
#include "mpz_from_padic.h"

/* This one is short and sweet.  Unfortunately it is a bit slow when  */
/* there are many digits.  The longer routine above is 4 times faster */
/* when there are about 16000 digits.                                 */

void mpz_from_padic(mpz_t ans, FFelem p, int k, FFelem *digits)
{
  int i;

  mpz_set_ui(ans, 0);
  for (i=k-1; i >= 0; i--)
  {
    mpz_mul_ui(ans, ans, (unsigned long)p);
    mpz_add_ui(ans, ans, (unsigned long)digits[i]);
  }
}
#endif
