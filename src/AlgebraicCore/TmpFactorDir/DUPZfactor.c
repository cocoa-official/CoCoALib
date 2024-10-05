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
#include <stdlib.h>

#include "DUPZfactor.h"

#include "logi.h"
#include "mpz_lift_recip.h"
#include "primes.h"
#include "FF.h"
#include "DUPZfactor_liftq.h"
#include "DUPZfactor_combine.h"
#include "DUPZ.h"
#include "DUPZfactor_info.h"
#include "DUPZfactor_prime.h"
#include "DUPZfactor_bound.h"
#include "DUPZlist.h"
#include "DUPZsqfrd.h"
#include "DUPZcontent.h"
#include "DUPZevaluate.h"


/***************************************************************************/

/* advance declarations */
static void DUPZfactor_sqfr(DUPZfactors ans, const DUPZ fin);
void DUPZfactor_finish(DUPZfactor_info THIS);



DUPZfactors DUPZfactor(const DUPZ fin)
{
  mpz_t content;
  DUPZfactors ans;
  DUPZlist sqfr_cmpts, iter;
  DUPZ f;
  int df;

  df = DUPZdeg(fin);
  ans = DUPZfactors_ctor();
  /* Zero polynomial -- 8 lines of code for a stupid case, groan! */
  if (df < 0)
  {
    mpz_t zero;
    mpz_init_set_ui(zero, 0);
    DUPZfactors_content(ans, fin->coeffs[0]);
    mpz_clear(zero);
    return ans;
  }
  /* Non-zero constant polynomial... */
  if (df == 0) { DUPZfactors_content(ans, fin->coeffs[0]); return ans; }

  /* remove content, and include it in the output if it is not 1 */
  f = DUPZcopy(fin);
  mpz_init(content);
  DUPZcontent(content, f);
  /* Fiddle the sign of the content so that DUPZsqfrd works as we want. */
  if (mpz_sgn(DUPZlc(f)) < 0) mpz_neg(content, content);
  DUPZdiv2z(f, content);
  DUPZfactors_content(ans, content);
  mpz_clear(content);

  /* square-free decompose f, then factorize each component */
  sqfr_cmpts = DUPZsqfrd(f);
  for (iter=sqfr_cmpts; iter; iter = iter->next)
  {
    DUPZfactors_multiplicity(ans, iter->deg);
    DUPZfactor_sqfr(ans, iter->poly);
  }

  DUPZlist_dtor(sqfr_cmpts);
  DUPZfree(f);
  return ans;
}

/***************************************************************************/
/* This function is to avoid problems of zero appearing when we evaluate   */
/* at a small random integer as a ruse to discard "bad combinations" fast. */
/* This code is simple rather than efficient; it should be fast enough.    */
/* Input f is assumed square-free.                                         */

static void DUPZfactor_remove_small_linears(DUPZfactors ans, DUPZ f)
{
  const int range=10; /* = random eval point range in DUPZfactor_combine */
  DUPZ linear, zero;
  int n;
  mpz_t y;

  mpz_init(y);
  linear = DUPZnew(1); mpz_set_ui(linear->coeffs[1], 1); linear->deg = 1;
  zero = DUPZnew(0);
  for (n=-range; n <=range; n++)
  {
    DUPZevaluate(y, f, n);
    if (mpz_sgn(y) != 0) continue;
    mpz_set_si(linear->coeffs[0], -n);
    DUPZfactors_add(ans, linear);
    DUPZdiv4(f, zero, f, linear);
  }
  
  DUPZfree(linear);
  DUPZfree(zero);
  mpz_clear(y);
}


/***************************************************************************/
/* In this function the input polynomial is content-free, and square-free. */
/* All irreducible factors of fin are added to ans.                        */
/* It is assumed that fin has degree at least 1.                           */
/* We remove small linear factors -- necessary for correct searching       */
/* through tuples (we need value of f at small integers to be non-zero).   */
/* Then decide whether to use forward or reverse polynomial by comparing   */
/* leading and trailing coeffs.                                            */

static void DUPZfactor_sqfr(DUPZfactors ans, const DUPZ fin)
{
  int irred;
  DUPZ f;
  DUPZfactor_info info;

  DUPZfactors_reversed(ans, 0); /* currently the normal way round */
  /* if degree == 1 then poly is irred */
  if (DUPZdeg(fin) == 1) { DUPZfactors_add(ans, fin); return; }

  /* Decide between forward and reversed polynomial.                 */
  /* This is a little messy as we also remove "small" linear factors */
  /* from both the forward and reversed polynomials.                 */
  f = DUPZcopy(fin);

  DUPZfactor_remove_small_linears(ans, f);
  if (DUPZdeg(f) == 0) /* if nothing left, just return. */
  {
    if (mpz_sgn(f->coeffs[0]) < 0) DUPZfactors_negate(ans);
    DUPZfree(f);
    return;
  }
  DUPZreverse(f);
  DUPZfactors_reversed(ans, 1);
  DUPZfactor_remove_small_linears(ans, f);
  if (DUPZdeg(f) == 0) /* if nothing left, just return. */
  {
    if (mpz_sgn(f->coeffs[0]) < 0) DUPZfactors_negate(ans);
    DUPZfree(f);
    return;
  }
  /* Here is where we really make the decision. */
  if (mpz_sizeinbase(DUPZlc(f), 2) > mpz_sizeinbase(f->coeffs[0], 2))
  {
    DUPZfactors_reversed(ans, 0);
    DUPZreverse(f);
  }

  /* Despite earlier juggling f may now have a negative leading coeff... */
  info = DUPZfactor_info_ctor(f, ans);          /* info now owns f */
  irred = DUPZfactor_pick_prime(info);
  if (irred) { DUPZfactors_add(ans, f); goto tidy_mem;} /* fds proves f to be irred */
  info->bounds = DUPZfactor_bound_ctor(f, info->dmax);
  DUPZfactor_info_set_target_height(info);
  DUPZfactor_lift_init(info);
  DUPZfactor_combine_init();
  DUPZfactor_finish(info);
  DUPZfactor_combine_done();
tidy_mem:
  DUPZfactor_info_dtor(info);
}

/***************************************************************************/
/* This function does both Hensel lifting and factor recombination.        */
/* The two are intertwined as we are using "early detection".              */


void DUPZfactor_finish(DUPZfactor_info THIS)
{
  int next_height;
  mpz_t nextQ;
  DUPZfactor_combiner combiner;
  DUPZ f_monic;

//  printf("Hensel lifting with early detection from height %d to height %d.\n", THIS->current_height, THIS->target_height);

  mpz_init(nextQ);
  f_monic = DUPZnew(DUPZdeg(THIS->f));
  while (THIS->current_height < THIS->target_height)
  {
    if (THIS->target_height > 16*THIS->current_height)
      next_height = 2*THIS->current_height;
    else
    {
      for (next_height = THIS->target_height; next_height > 2*THIS->current_height; next_height = (next_height+1)/2) {}
    }
    mpz_ui_pow_ui(THIS->Q, THIS->p, (next_height+1)/2);
    mpz_mul(nextQ, THIS->Q, THIS->Q);
    mpz_lift_recip(THIS->recip_lcf, DUPZlc(THIS->f), nextQ);
    DUPZmmul3z(f_monic, THIS->f, THIS->recip_lcf, nextQ);
    DUPZfactor_lift_step(THIS->lifter, f_monic, THIS->Q);
    mpz_set(THIS->Q, nextQ);
    THIS->current_height = next_height;

    /* Perform a speculative early search unless lifting is complete */
    if (THIS->current_height >= THIS->target_height) break;
continue; /* disable early detection for the moment */
//if (THIS->current_height < 3 + logi(DUPZdeg(THIS->f)) + THIS->bounds->root_bound) continue;
//???    if (THIS->current_height < 4) continue;
#ifdef FACTOR_DEBUG
    printf("\nStarting speculative early search at height %d.\n", THIS->current_height);
#endif
    combiner = DUPZfactor_combiner_ctor(THIS);
    DUPZfactor_combine_early(combiner);
    if (DUPZdeg(THIS->f) > 0 &&
	combiner->nfactors != combiner->nfactors_orig)
    {
      DUPZfactor_lift_revise(THIS);
#ifdef FACTOR_DEBUG
      printf("Revised lifting bound is %d\n", THIS->target_height);
#endif
    }
    DUPZfactor_combiner_dtor(combiner);
    if (DUPZdeg(THIS->f) <= 0) break;

  }
  mpz_clear(nextQ);
  DUPZfree(f_monic);

  if (DUPZdeg(THIS->f) <= 0) return;
#ifdef FACTOR_DEBUG
  printf("\nSome factors remain; starting final search\n");
#endif

  /* try combinations */

  combiner = DUPZfactor_combiner_ctor(THIS);
  DUPZfactor_combine_final(combiner);
  DUPZfactor_combiner_dtor(combiner);
}

