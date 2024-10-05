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

#include "DUPZ_DUPFF.h"
#include "jaaerror.h"

DUPZ DUPFF_to_DUPZ(const DUPFF f)
{
  DUPZ ans;
  int df, i;

  df = DUPFFdeg(f);
  ans =  DUPZnew(df);
  for (i=0; i <= df; i++) mpz_set_ui(ans->coeffs[i], f->coeffs[i]);
  ans->deg = df;

  return ans;
}


void DUPFF_to_DUPZ2(DUPZ dest, const DUPFF f)
{
  int df, i;

  df = DUPFFdeg(f);
  if (df > dest->maxdeg) JERROR(JERROR_DEG_TOO_LOW);
  for (i=0; i <= df; i++) mpz_set_ui(dest->coeffs[i], f->coeffs[i]);
  dest->deg = df;
}


DUPFF DUPZ_to_DUPFF(const DUPZ f)
{
  DUPFF ans;
  int i, p, df;

  p = CurrentFF.prime;
  df = DUPZdeg(f);
  /* find the degree of ans */
  for (i=df; (i >= 0) && (mpz_fdiv_ui(f->coeffs[i], p) == 0); i--) ;
  ans = DUPFFnew(i);
  ans->deg = i;
  for (; i >= 0; i--)
    ans->coeffs[i] = mpz_fdiv_ui(f->coeffs[i], p);

  return ans;
}


void DUPZ_to_DUPFF2(DUPFF fbar, const DUPZ f)
{
  int i, p, df;

  p = CurrentFF.prime;
  df = DUPZdeg(f);
  /* check the degree of fbar */
  for (i=df; (i >= 0) && (mpz_fdiv_ui(f->coeffs[i], p) == 0); i--) ;
  if (fbar->maxdeg < i) { JERROR(JERROR_DEG_TOO_LOW); return; }
  fbar->deg = i;
  for (; i >= 0; i--)
    fbar->coeffs[i] = mpz_fdiv_ui(f->coeffs[i], p);
}
