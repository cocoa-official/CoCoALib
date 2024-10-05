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

#include "DUPFFderiv.h"
#include "jaaerror.h"

void DUPFFformal_deriv2(DUPFF fprime, const DUPFF f)
/* fprime and f may be identical */
{
  int df, dans, j, p;

  p = CurrentFF.prime;
  df = DUPFFdeg(f);
  /* find degree of the answer */
  for (dans=df; dans >= 1; dans--)
    if (f->coeffs[dans] != 0 && (dans%p != 0)) break;
  dans--;
  if (dans < 0)
  {
    fprime->deg = -1;
    return;
  }
  if (fprime->maxdeg < dans)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    return;
  }
  fprime->deg = dans;
  for (j = 1; j <= dans+1; j++)
    fprime->coeffs[j-1] = FFmul(j%p, f->coeffs[j]);
  return;
}
