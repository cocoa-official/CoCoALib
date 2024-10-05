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

#include "DUPZcra.h"
#include "jaaerror.h"
#include "mpz_cra_ui.h"

/**********************************************************************/
/* NOTE the CRA construction is symmetric (-modulus/2,..., modulus/2) */
/**********************************************************************/

int DUPZcra(DUPZ f1, const mpz_t m1, const DUPFF f2, FFelem m2)
{
  int i, changed, df1, df2;
  mpz_t tmp; /* used as workspace in mpz_cra_ui_raw */
  FFelem m1modm2;

  df1 = DUPZdeg(f1);
  df2 = DUPFFdeg(f2);
  if (df2 > f1->maxdeg) { JERROR(JERROR_DEG_TOO_LOW); return 1; }

  mpz_init(tmp);
  m1modm2 = mpz_fdiv_ui(m1, (unsigned long)m2);
  changed = 0;
  for (i=1+df1; i <= df2; i++) mpz_set_ui(f1->coeffs[i], 0);
  for (i=0; i <= df2; i++)
    changed |= mpz_cra_ui_raw(f1->coeffs[i], m1, f2->coeffs[i], m2, tmp, m1modm2);
  for (i=1+df2; i <= df1; i++)
    changed |= mpz_cra_ui_raw(f1->coeffs[i], m1, 0, m2, tmp, m1modm2);

  if (df2 > df1) f1->deg = df2;
  mpz_clear(tmp);
  return changed;
}
