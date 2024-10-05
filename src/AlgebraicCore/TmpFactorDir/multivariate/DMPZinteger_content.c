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

#include "DMPZinteger_content.h"
#include <stddef.h>

/* Compute integer content of a DMPZ polynomial and put result in content. */
/* The algorithm is not terribly clever.                                   */
/* If f = 0 result is zero, otherwise it is positive.                      */

void DMPZinteger_content(mpz_t content, DMPZ f)
{
  DMPZ iter, min_ptr;
  size_t min;

  if (f == NULL) { mpz_set_ui(content, 0); return; } /* silly case, f==0 */
  /* Start content with value of the "smallest" coefficient. */
  min_ptr = f;
  min = mpz_sizeinbase(f->coeff, 2);
  for (iter=f; iter; iter = iter->next)
    if (mpz_sizeinbase(iter->coeff, 2) < min)
    {
      min = mpz_sizeinbase(iter->coeff, 2);
      min_ptr = iter;
    }
  mpz_abs(content, min_ptr->coeff);

  /* Now gcd content with all the other coefficients. */
  for (iter=f; (mpz_cmp_ui(content, 1) != 0) && iter; iter = iter->next)
    if (iter != min_ptr)
      mpz_gcd(content, content, iter->coeff);
}
