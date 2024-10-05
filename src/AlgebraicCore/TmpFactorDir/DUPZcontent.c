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

#include "DUPZcontent.h"


void DUPZcontent(mpz_t ans, const DUPZ f)
{
  int i, min, minsz, sz, df;

  df = DUPZdeg(f);
  /* Search for the "smallest" non-zero coefficient as our initial value */
  min = 0; minsz = mpz_sizeinbase(f->coeffs[0], 2);
  for (i=1; i <= df; i++)
  {
    if (mpz_sgn(f->coeffs[i]) == 0) continue;
    sz = mpz_sizeinbase(f->coeffs[i], 2);
    if (sz < minsz) { minsz = sz; min = i; }
  }

  mpz_set(ans, f->coeffs[min]);
  for (i=0; i <= df; i++)
  {
    if (i == min) continue; /* we started with this value so can skip it now */
    if (mpz_cmp_ui(ans, 1) == 0) break;  /* answer is 1 */
    if (mpz_sgn(f->coeffs[i]) == 0) continue; /* ignore zero coeffs */
    mpz_gcd(ans, ans, f->coeffs[i]);
  }
}
