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

/* Email sent to GMP people hoping that this will be implemented more */
/* sensibly/efficiently in some later release.  28 Oct 1997.          */

#include "mpz_cmp_abs.h"

int mpz_cmp_abs(mpz_t A, mpz_t B)
{
  mpz_t tmp1, tmp2;
  int result;

  mpz_init(tmp1);
  mpz_init(tmp2);
  mpz_abs(tmp1, A);
  mpz_abs(tmp2, B);
  result = mpz_cmp(tmp1, tmp2);
  mpz_clear(tmp2);
  mpz_clear(tmp1);
  return result;
}


int mpz_cmp_abs_ui(mpz_t A, long B)
{
  mpz_t tmp1;
  int result;

  mpz_init(tmp1);
  mpz_abs(tmp1, A);
  result = mpz_cmp_ui(tmp1, B);
  mpz_clear(tmp1);
  return result;
}
