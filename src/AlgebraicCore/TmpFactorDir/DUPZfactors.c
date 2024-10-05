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

#include <stddef.h>
#include "DUPZfactors.h"
#include "jalloc.h"
#include "jaaerror.h"


DUPZfactors DUPZfactors_ctor()
{
  DUPZfactors THIS;

  THIS = (DUPZfactors)MALLOC(sizeof(struct DUPZfactors_struct));
  THIS->list = NULL;
  THIS->multiplicity = 1;
  THIS->reversed = 0;
  mpz_init_set_ui(THIS->content, 1);
  return THIS;
}


void DUPZfactors_dtor(DUPZfactors THIS)
{
  mpz_clear(THIS->content);
  DUPZlist_dtor(THIS->list);
  FREE(THIS);
}

void DUPZfactors_content(DUPZfactors THIS, mpz_t n)
{
  mpz_mul(THIS->content, THIS->content, n);
}


void DUPZfactors_negate(DUPZfactors THIS)
{
  if (THIS->multiplicity & 1) mpz_neg(THIS->content, THIS->content);
}


void DUPZfactors_add(DUPZfactors THIS, const DUPZ f)
{
  DUPZ fcopy;
  
#ifdef FACTOR_DEBUG
  printf("%d ", DUPZdeg(f));
#endif
  fcopy = DUPZcopy(f);
  if (THIS->reversed) DUPZreverse(fcopy);
  /* The next lines just make sure that the lc of each factor is >0       */
  /* Also negate the "content" if lc < 0 and the multiplicity is odd.     */
  if (mpz_sgn(DUPZlc(fcopy)) < 0)
  {
    DUPZneg1(fcopy);
    if (THIS->multiplicity&1) mpz_neg(THIS->content, THIS->content);
  }
  THIS->list = DUPZlist_append(THIS->list,
                               DUPZlist_ctor(fcopy, THIS->multiplicity));
}

void DUPZfactors_multiplicity(DUPZfactors THIS, int m)
{
  THIS->multiplicity = m;
}

void DUPZfactors_reversed(DUPZfactors THIS, int yes_no)
{
  THIS->reversed = yes_no;
}

int DUPZfactors_nfactors(const DUPZfactors THIS)
{
  DUPZlist iter;
  int nfactors;

  nfactors = 0;
  for (iter=THIS->list; iter; iter = iter->next)
    nfactors += iter->deg;
  return nfactors;
}
