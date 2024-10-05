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

#ifndef DMPZfactors_H
#define DMPZfactors_H

#include "DMPZ.h"
#include "DMPZlist.h"

/* To make sure that the product of the factors in the list multiplied    */
/* by the content is equal to the original polynomial (with correct sign) */
/* the sign of the leading coefficient (wrt to deglex) of the original    */
/* polynomial must be included in the integer content; all polynomial     */
/* factors are forced to have positive leading coefficient (wrt deglex).  */

struct DMPZfactors_struct
{
  mpz_t content;    /* signed integer content */
  DMPZlist list;    /* list of factors, each tagged with its multiplicity */
  int multiplicity; /* used by DMPZfactors_add to know multiplicity of factor */
};

typedef struct DMPZfactors_struct *DMPZfactors;


DMPZfactors DMPZfactors_ctor();
void DMPZfactors_dtor(DMPZfactors A);
void DMPZfactors_content(DMPZfactors A, mpz_t n);
void DMPZfactors_negate(DMPZfactors A);
void DMPZfactors_add(DMPZfactors A, const DMPZ f);
void DMPZfactors_multiplicity(DMPZfactors A, int m);
//void DMPZfactors_reversed(DMPZfactors A, int yes_no);
int DMPZfactors_nfactors(const DMPZfactors A);


#endif
