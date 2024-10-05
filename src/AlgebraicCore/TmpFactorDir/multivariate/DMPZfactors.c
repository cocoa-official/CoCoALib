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

#include "DMPZfactors.h"
#include "jalloc.h"


DMPZfactors DMPZfactors_ctor()
{
  DMPZfactors THIS;

  THIS = (DMPZfactors)MALLOC(sizeof(struct DMPZfactors_struct));
  THIS->list = NULL;
  THIS->multiplicity = 1;
  mpz_init_set_ui(THIS->content, 1);
  return THIS;
}


void DMPZfactors_dtor(DMPZfactors THIS)
{
  mpz_clear(THIS->content);
  DMPZlist_dtor(THIS->list);
  FREE(THIS);
}

void DMPZfactors_content(DMPZfactors THIS, mpz_t n)
{
  mpz_mul(THIS->content, THIS->content, n);
}


void DMPZfactors_negate(DMPZfactors THIS)
{
  if (THIS->multiplicity & 1) mpz_neg(THIS->content, THIS->content);
}


void DMPZfactors_add(DMPZfactors THIS, const DMPZ f)
{
  DMPZ fcopy;
  
  fcopy = DMPZcopy(f);
  /* force LC (wrt to deglex) > 0; order deglex imposed by DMPZfactor(..) */
  fcopy = DMPZsort(fcopy, DMPZorder_deglex);
  if (mpz_sgn(fcopy->coeff) < 0) DMPZnegate(fcopy);
  THIS->list = DMPZlist_append(THIS->list,
                               DMPZlist_ctor(fcopy, THIS->multiplicity));
}

void DMPZfactors_multiplicity(DMPZfactors THIS, int m)
{
  THIS->multiplicity = m;
}

int DMPZfactors_nfactors(const DMPZfactors THIS)
{
  DMPZlist iter;
  int nfactors;

  nfactors = 0;
  for (iter=THIS->list; iter; iter = iter->next)
    nfactors += iter->deg;
  return nfactors;
}
