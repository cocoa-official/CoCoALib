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

#include <stdio.h>
#include "DUPZprint.h"

void DUPZprint(const DUPZ f)
{
  int i, df;

  df = DUPZdeg(f);
  if (df < 0) { printf("0;\n"); return; }
  if (df > 99) df = 99; /* avoid printing out too much */
  printf(" ");
  for (i=0; i<df ; i++)
  {
    if (mpz_sgn(f->coeffs[i]) == 0) continue;
    mpz_out_str(stdout, 10, f->coeffs[i]);
    printf("*x^%d+ ", i);
  }
  if (DUPZdeg(f) > 99) printf(".... +");
  mpz_out_str(stdout, 10, f->coeffs[DUPZdeg(f)]);
  printf("*x^%d;\n", DUPZdeg(f));
}
