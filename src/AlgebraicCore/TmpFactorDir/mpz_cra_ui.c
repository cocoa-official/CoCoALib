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

#include "mpz_cra_ui.h"
#include "FF.h"

int mpz_cra_ui_raw(mpz_t r1, const mpz_t m1, FFelem r2, FFelem m2,
		   mpz_t tmp, FFelem m1modm2)
{
  FFelem r1modm2, k;
  int sign;

  if (mpz_cmp_ui(m1, 1) == 0) mpz_set_ui(r1, 0); /* if m1=1, make sure r1=0 */
  r1modm2 = mpz_fdiv_ui(r1, (unsigned long)m2);
  k = FFdiv(FFsub(r2, r1modm2), m1modm2);
  sign = 1;
  if (k > m2/2) { k = m2-k; sign = -1; }
  mpz_mul_ui(tmp, m1, (unsigned long)k);
  if (sign == -1) mpz_neg(tmp, tmp);
  mpz_add(r1, r1, tmp);

  return (r1modm2 != r2);
}


int mpz_cra_ui(mpz_t r1, const mpz_t m1, FFelem r2, FFelem m2)
{
  mpz_t tmp;
  int changed;

  mpz_init(tmp);
  changed = mpz_cra_ui_raw(r1, m1, r2, m2, tmp, mpz_fdiv_ui(m1, (unsigned long)m2));
  mpz_clear(tmp);
  return changed;
}
