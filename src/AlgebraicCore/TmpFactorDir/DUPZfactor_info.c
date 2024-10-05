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
#include "DUPZfactor_info.h"
#include "logi.h"
#include "jalloc.h"
#include "DUPZfactor_liftq.h"

DUPZfactor_info DUPZfactor_info_ctor(DUPZ f, DUPZfactors output)
{
  int i, df;
  DUPZfactor_info ans;

  ans = (DUPZfactor_info)MALLOC(sizeof(struct DUPZfactor_info_struct));
  ans->f = f; /* deliberately do not copy */
  ans->irreds = output;
  ans->nprimes = 0;
  ans->max_nprimes = 25;
  ans->FFq = (FF*)MALLOC(ans->max_nprimes*sizeof(struct FFstruct));
  ans->qfactors = (DUPFFlist*)MALLOC(ans->max_nprimes*sizeof(DUPFFlist));
  df = DUPZdeg(f);
  ans->fds = (int*)MALLOC((1+df)*sizeof(int));
  for (i=0; i <= df; i++) ans->fds[i] = 1;

  mpz_init(ans->Q);
  mpz_init(ans->recip_lcf);
  /* Set these fields to NULL to mean "unused" */
  ans->bounds = NULL;
  ans->lifter = NULL;

  return ans;
}


void DUPZfactor_info_dtor(DUPZfactor_info THIS)
{
  int i;

  mpz_clear(THIS->recip_lcf);
  mpz_clear(THIS->Q);
  DUPZfree(THIS->f);
  for (i=0; i < THIS->nprimes; i++)
  {
    FFdtor(THIS->FFq[i]);
    DUPFFlist_dtor(THIS->qfactors[i]);
  }
  FREE(THIS->FFq);
  FREE(THIS->qfactors);
  if (THIS->fds) FREE(THIS->fds);
  /* The two fields below may be "unused", if so don't call their dtors */
  if (THIS->bounds) DUPZfactor_bound_dtor(THIS->bounds);
  if (THIS->lifter) DUPZfactor_lift_dtor(THIS->lifter);
  FREE(THIS);
}

/***************************************************************************/
/* Set the target_height based on information in bounds and the prime.     */

void DUPZfactor_info_set_target_height(DUPZfactor_info THIS)
{
  THIS->target_height = 1 + (int)((logi(2) + THIS->bounds->leading_coeff + THIS->bounds->lift_bound)/logi(THIS->p));
}
