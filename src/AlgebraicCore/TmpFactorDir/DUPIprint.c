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

#include "DUPIprint.h"
#include <stdio.h>

void DUPIprint(const DUPI f)
{
  int i, df;

  df = DUPIdeg(f);
  if (df < 0) { printf("0\n"); return; }
  if (df > 99) df = 99; /* avoid printing out too much */
  printf("[");
  for (i=0; i<df ; i++) printf("%d, ", f->coeffs[i]);
  if (DUPIdeg(f) > 99) printf(".... ,");
  printf("%d]\n", f->coeffs[DUPIdeg(f)]);
}
