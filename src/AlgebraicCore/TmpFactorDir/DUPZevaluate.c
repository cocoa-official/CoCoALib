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

#include "DUPZevaluate.h"

void DUPZevaluate(mpz_t value, const DUPZ f, int x)
{
  int i, df, negative;

  df = DUPZdeg(f);
  mpz_set_ui(value, 0);
  negative = x < 0;
  if (negative) x = -x;
  for(i=df; i >= 0; i--)
  {
    mpz_mul_ui(value, value, x);
    if (negative) mpz_neg(value, value);
    mpz_add(value, value, f->coeffs[i]);
  }
}

