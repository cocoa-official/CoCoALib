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

#include "DUPZformal_deriv.h"

#include "gmp.h"
#include "jaaerror.h"

DUPZ DUPZformal_deriv(const DUPZ f)
{
  DUPZ ans;
  int df = DUPZdeg(f);
  
  if (df < 1) return DUPZnew(-1);
  ans = DUPZnew(df-1);
  DUPZformal_deriv2(ans, f);
  return ans;
}


/* f and f prime may be identical */
void DUPZformal_deriv2(DUPZ fprime, const DUPZ f)
{
  int i;
  const int df = DUPZdeg(f);

  if (df < 1) { fprime->deg = -1; return; }
  if (df - 1 > fprime->maxdeg) { JERROR(JERROR_DEG_TOO_LOW); return; }
  fprime->deg = df-1;
  for (i=1; i <= df; i++)
    mpz_mul_ui(fprime->coeffs[i-1], f->coeffs[i], i);
}
