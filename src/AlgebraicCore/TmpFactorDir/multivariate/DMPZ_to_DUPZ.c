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

#include "DMPZ_to_DUPZ.h"
#include "jalloc.h"


DUPZ DMPZ_to_DUPZ(const DMPZ f, int var)
{
  int df, e;
  DUPZ ans;
  DMPZ term;

  df = DMPZdeg(f, var);
  ans = DUPZnew(df);
  /* DMPZs are not canonical, so same exponent may appear several times */
  for (term = f; term; term = term->next)
  {
    e = term->exps[var];
    mpz_add(ans->coeffs[e], ans->coeffs[e], term->coeff);
  }
  while (df >= 0 && mpz_sgn(ans->coeffs[df]) == 0) df--;
  ans->deg = df;
  return ans;
}


DMPZ DUPZ_to_DMPZ(const DUPZ f, int var)
{
  DMPZ ans;
  int i, j, *exps;

  ans = NULL;
  for (i=0; i <= DUPZdeg(f); i++)
  {
    if (mpz_sgn(f->coeffs[i]) == 0) continue;
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (j=0; j < NVARS; j++) if (j == var) exps[j] = i; else exps[j] = 0;
    ans = DMPZprepend(f->coeffs[i], exps, ans);
  }
  return ans;
}

