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

#ifndef DUPZPADIC_EXPAND_H
#define DUPZPADIC_EXPAND_H

#include "DUPFF.h"
#include "DUPZ.h"
#include "DUPI.h"
#include <gmp.h>


/* Using a struct avoid polluting the top level namespace.                 */
struct DUPZpadic_struct
{
  int p, k, kmax;
  mpz_t lcf_copy;
  mpz_t tmp;
  int* lcf_cmpts;
  DUPZ f_copy;
  DUPFF *f_cmpts, *monic_cmpts;
  DUPFF deg0;   /* used to view constants as polynomials */
  DUPI E, fk;
};

typedef struct DUPZpadic_struct *DUPZpadic_expansion;


DUPZpadic_expansion DUPZpadic_expand_init(const DUPZ f, int p, int kmax);
int DUPZpadic_expand_step(DUPZpadic_expansion z);
void DUPZpadic_expand_end(DUPZpadic_expansion z);

#endif
