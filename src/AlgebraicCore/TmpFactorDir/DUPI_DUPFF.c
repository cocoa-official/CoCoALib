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

#include "DUPI_DUPFF.h"
#include "jaaerror.h"


DUPI DUPFF_to_DUPI(const DUPFF f)
{
  DUPFF tmp;

  tmp = DUPFFcopy(f);
  return (DUPI)tmp;
}


void DUPFF_to_DUPI2(DUPI f, const DUPFF g)
{
  int i;
  const int dg = DUPFFdeg(g);

  if (f->maxdeg < dg) { JERROR(JERROR_DEG_TOO_LOW); return; }
  f->deg = dg;
  for (i=0; i <= dg; i++) f->coeffs[i] = g->coeffs[i];
}


DUPFF DUPI_to_DUPFF(const DUPI f)
{
  DUPFF ans;
  int i, df;
  const int p = CurrentFF.prime;

  df = DUPIdeg(f);
  while (df >= 0 && f->coeffs[df] % p == 0) df--;
  ans = DUPFFnew(df);
  ans->deg = df;
  for (i=0; i <= df; i++)
  {
    if (f->coeffs[i] >= 0) ans->coeffs[i] = f->coeffs[i]%p;
    else
    {
      ans->coeffs[i] = p - ((- f->coeffs[i])%p);
      if (ans->coeffs[i] == p) ans->coeffs[i] = 0;
    }
  }
 
  return ans;
}


void DUPI_to_DUPFF2(DUPFF f, const DUPI g)
{
  int i, dg;
  const int p = CurrentFF.prime;

  dg = DUPIdeg(g);
  while (dg >= 0 && g->coeffs[dg] % p == 0) dg--;
  if (f->maxdeg < dg) { JERROR(JERROR_DEG_TOO_LOW); return; }
  f->deg = dg;
  for (i=0; i <= dg; i++)
  {
    if (g->coeffs[i] >= 0) f->coeffs[i] = g->coeffs[i]%p;
    else
    {
      f->coeffs[i] = p - ((- g->coeffs[i])%p);
      if (f->coeffs[i] == p) f->coeffs[i] = 0;
    }
  }
}
