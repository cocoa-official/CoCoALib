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

#ifndef DUPZfactor_info_h
#define DUPZfactor_info_h

#include "DUPZ.h"
#include "DUPZfactors.h"
#include "DUPFFlist.h"
#include "FF.h"
#include "DUPZfactor_bound.h"
#include "DUPZfactor_liftq.h"

#include "DUPZfactor_info1.h"
/* typedef struct DUPZfactor_info_struct *DUPZfactor_info; */

struct DUPZfactor_info_struct
{
  DUPZ f;                           /* poly to be factorized           */
  DUPZfactors irreds;               /* where to put factors            */
  int p;                            /* chosen prime modulus            */
  int p_index;                      /* index into FFq and qfactors     */
  DUPFFlist pfactors;               /* factors modulo p                */
  int nprimes;                      /* number of primes used           */
  int max_nprimes;                  /* max number of primes to be tried*/
  FF *FFq;                          /* finite fields for those primes  */
  DUPFFlist *qfactors;              /* several modular factorizations  */
  int *fds;                         /* factor degree set               */
  int dmax;                         /* maximum degree of any factor    */

  bound_info bounds;                /* bound information for poly      */
  int target_height;                /* a priori lift height needed     */
  int current_height;               /* height of factors in lifter     */
  mpz_t recip_lcf;                  /* 1/lc(f) modulo Q                */
  mpz_t Q;                          /* p^height                        */
  DUPZfactor_lifter lifter;         /* tree for Hensel lifting         */
};


DUPZfactor_info DUPZfactor_info_ctor(DUPZ f, DUPZfactors output);
void DUPZfactor_info_dtor(DUPZfactor_info THIS);
void DUPZfactor_info_set_target_height(DUPZfactor_info THIS);

#endif
