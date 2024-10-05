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

#ifndef DUPZFACTORS_H
#define DUPZFACTORS_H

#include "gmp.h"
#include "DUPZ.h"
#include "DUPZlist.h"

struct DUPZfactors_struct
{
  mpz_t content;
  DUPZlist list;
  int multiplicity;
  int reversed;  /* flag to say whether to reverse factors or not */
};

typedef struct DUPZfactors_struct *DUPZfactors;


DUPZfactors DUPZfactors_ctor();
void DUPZfactors_dtor(DUPZfactors A);
void DUPZfactors_content(DUPZfactors A, mpz_t n);
void DUPZfactors_negate(DUPZfactors A);
void DUPZfactors_add(DUPZfactors A, const DUPZ f);
void DUPZfactors_multiplicity(DUPZfactors A, int m);
void DUPZfactors_reversed(DUPZfactors A, int yes_no);
int DUPZfactors_nfactors(const DUPZfactors A);


#endif
