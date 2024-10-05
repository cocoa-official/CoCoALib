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


#include "DUPZfactor_combine.h"
#include "mpz_log.h"
/* We use a random number in DUPZfactor_combine_init, so include stdlib.h */
#include <stdlib.h>
#include <stddef.h>
/* Use float.h just to get DBL_EPSILON and DBL_MANT_DIG, used in DUPZfactor_combiner_set_n1coeffd_bounds */
#include <float.h>

#include <math.h>
#include "jalloc.h"
#include "logi.h"
#include "WGD.h"
#include "DUPZevaluate.h"
#include "DUPZfactor_bound.h"
#include "DUPFF.h"
#include "DUPZ_DUPFF.h"
#include "DUPZfactor_refine_fds.h"

/* This is defined in DUPZfactor.c. */
extern void DUPZfactor_finish(DUPZfactor_info THIS);

/***************************************************************************/
/* These static variables are initialised by DUPZfactor_combine_init.      */

static mpz_t num, den, tmp, lcm;

/* for accounting/debugging */
static int divisions;
static int count1, count2, count3, count4, count5, count6;

/***************************************************************************/
/* This function sets values for P and Q the numerator and denominator     */
/* bounds to be used during rational reconstruction.  During early         */
/* searching the value for Q is "guessed" (approx 1/(2*sqrt(lc(f)))).      */

static void DUPZfactor_combiner_set_PQ(DUPZfactor_combiner THIS)
{
  /* numerator/denominator bounds for rational reconstruction */
  /* Updated whenever non-monic factors are found             */
  if (mpz_log(THIS->modulus) > logi(2) + 2*mpz_log(DUPZlc(THIS->info->f)))
    mpz_abs(THIS->Q, DUPZlc(THIS->info->f));
  else
  {
    mpz_set_ui(THIS->Q, 1);
    mpz_mul_2exp(THIS->Q, THIS->Q, mpz_sizeinbase(THIS->modulus, 4)-1);
  }
  
  mpz_fdiv_q(THIS->P, THIS->modulus, THIS->Q);
  mpz_fdiv_q_ui(THIS->P, THIS->P, 2);
}

#if 0
static void DUPZfactor_combiner_set_n1coeff_bounds(DUPZfactor_combiner THIS)
{
  mpz_t tmp;
  int i;
  double log_bound;

  mpz_init(tmp);
  for (i=0; i < nfactors; i++)
  {
    fi = *THIS->factors[i];
    mpz_mul(tmp, fi->coeff[degfi - delta], DUPZlc(THIS->info->f));
    mpz_fdiv_r(tmp, tmp, modulus);
    mpz_mul_2exp(tmp, tmp, 32);    /* assumes unsigned is 32 bits */
    mpz_fdiv_q(tmp, tmp, modulus);
    n1coeff[i] = mpz_get_ui(tmp);
  }

  log_bound = THIS->info->bounds->root_bound + THIS->info->bounds->leading_coeff - THIS->height * logi(THIS->p);
  if (log_bound < -32*logi(2)) n1bound = 1;
  else n1bound = 1 + exp(32*logi(2) + log_bound);
  mpz_clear(tmp);
}
#endif

static void DUPZfactor_combiner_set_n1coeffd_bounds(DUPZfactor_combiner THIS)
{
/*  const int nfactors = DUPFFlist_length(THIS->info->pfactors);*/
  int log2_modulus, shift;
  double tmpd, log_modulus;
  const int float_bits = DBL_MANT_DIG;
  const double eps = DBL_EPSILON;
  DUPZ fi;
  mpz_t n1coeff, half_modulus; /* workspace */
  int i, delta;

  log_modulus = mpz_log(THIS->modulus);
  mpz_init(n1coeff);
  mpz_init_set(half_modulus, THIS->modulus);
  mpz_fdiv_q_ui(half_modulus, half_modulus, 2);
  log2_modulus = mpz_sizeinbase(THIS->modulus, 2);
  shift = log2_modulus - float_bits - 10; /* 10 is rather arbitrary */
  if (shift < 0) shift = 0;
delta=1;
  for (i=0; i < THIS->nfactors_orig; i++)
  {
    if (THIS->used[i]) continue;
    fi = *THIS->factors[i];
    mpz_set(n1coeff, fi->coeffs[DUPZdeg(fi)-delta]);
/*printf("d-1 coeff of f[%d]=",i);mpz_out_str(stdout,10,n1coeff);printf("\n");*/
    mpz_mul(n1coeff, n1coeff, DUPZlc(THIS->info->f));
    mpz_fdiv_r(n1coeff, n1coeff, THIS->modulus);
    if (mpz_cmp(n1coeff, half_modulus) > 0)
      mpz_sub(n1coeff, n1coeff, THIS->modulus);
    mpz_tdiv_q_2exp(n1coeff, n1coeff, shift);
    THIS->n1coeffd[i] = mpz_get_d(n1coeff);
    mpz_tdiv_q_2exp(n1coeff, THIS->modulus, shift);
    THIS->n1coeffd[i] /= mpz_get_d(n1coeff);
/*DUPZprint(fi);*/
/*printf("THIS->n1coeffd[%d]=%g\n",i,THIS->n1coeffd[i]);*/
  }
  mpz_clear(n1coeff);
  mpz_clear(half_modulus);
  THIS->n1coeffd_bound = exp(THIS->info->bounds->leading_coeff + THIS->info->bounds->root_bound - mpz_log(THIS->modulus));
  THIS->n1coeffd_bound += 2*eps;
  THIS->n1coeffd_sum[0] = 0;
/*
 *printf("mpz_log(THIS->modulus)=%g\n",mpz_log(THIS->modulus));
 *printf("THIS->info->bounds->root_bound=%g\n",THIS->info->bounds->root_bound);
 *printf("THIS->info->bounds->leading_coeff=%g\n",THIS->info->bounds->leading_coeff);
 *printf("THIS->n1coeffd_bound=%g\n",THIS->n1coeffd_bound);
 */
  THIS->n1coeffd_flag = 0;
  tmpd = THIS->info->bounds->leading_coeff + THIS->info->bounds->root_bound;
  if (log_modulus < logi(10) + logi(DUPZdeg(THIS->info->f)) + tmpd) return;
  THIS->n1coeffd_flag = 1;
}



DUPZfactor_combiner DUPZfactor_combiner_ctor(DUPZfactor_info info)
{
  DUPZ f = info->f;  /* aliasing for convenience */
  DUPZfactor_combiner THIS;
  int i;
  double rt_bd;
  const int nfactors = DUPFFlist_length(info->pfactors);
  const int range = 10; /* should be the same as range in DUPZfactor_sqfr */

  THIS = (DUPZfactor_combiner)MALLOC(sizeof(struct DUPZfactor_combine_struct));
  THIS->info = info;
  THIS->tuple_deg = 0;
  
  mpz_init_set(THIS->modulus, info->Q);
  THIS->nfactors_orig = nfactors; /* this never changes */
  THIS->nfactors = nfactors;      /* this counts down as factors are used */
  THIS->factors = (DUPZ**)MALLOC(nfactors*sizeof(DUPZ*));
  DUPZfactor_lift_output(THIS->factors, info->lifter, 0);
/*
 *printf("There are %d modular factors:\n", nfactors);
 *for(i=0;i<nfactors;i++)DUPZprint(*THIS->factors[i]);
 */
  THIS->combination = (int*)MALLOC(nfactors*sizeof(int));
  THIS->used = (int*)MALLOC(nfactors*sizeof(int));
  for (i=0; i < nfactors; i++) THIS->used[i] = 0;

  THIS->n1coeffd = (double*)MALLOC(nfactors*sizeof(double));
  THIS->n1coeffd_sum = (double*)MALLOC(nfactors*sizeof(double));
  DUPZfactor_combiner_set_n1coeffd_bounds(THIS);

  /* pick random value from [-range, range] avoiding 0 */
  THIS->evaluation_point = rand()%(2*range) - range;
  if (THIS->evaluation_point == 0) THIS->evaluation_point = range;

  /* the two lines below are odd to avoid risk of overflow */
  rt_bd = info->bounds->root_bound;
  THIS->log_ev_pt = rt_bd + log(1+exp(logi(abs(THIS->evaluation_point)) - rt_bd));
  mpz_init(THIS->f_value);
  DUPZevaluate(THIS->f_value, f, THIS->evaluation_point);
  /* we have removed small linear factors from f, so f_value is non-zero */
  THIS->factor_value = (mpz_t*)MALLOC(nfactors*sizeof(mpz_t));
  for (i=0; i < nfactors; i++)
  {
    mpz_init(THIS->factor_value[i]);
    DUPZevaluate(THIS->factor_value[i], *THIS->factors[i], THIS->evaluation_point);
  }

  mpz_init(THIS->P);
  mpz_init(THIS->Q);
  DUPZfactor_combiner_set_PQ(THIS);
  THIS->factor = DUPZnew(DUPZdeg(f)); /* Must be big enough to contain either the */
  THIS->quot = DUPZnew(DUPZdeg(f));   /* putative factor or its complement.       */
  THIS->rem = DUPZnew(DUPZdeg(f));
  mpz_init(THIS->value_at_1);
  DUPZevaluate(THIS->value_at_1, f, 1);
  mpz_init(THIS->value_at_m1);
  DUPZevaluate(THIS->value_at_m1, f, -1);

  THIS->tc_combination = (int*)MALLOC(nfactors*sizeof(int));
  THIS->tc_combination[0] = -1;
  THIS->tc_product = (mpz_t*)MALLOC(nfactors*sizeof(mpz_t));
  for (i=0; i < nfactors; i++) mpz_init(THIS->tc_product[i]);

  return THIS;
}


void DUPZfactor_combiner_dtor(DUPZfactor_combiner THIS)
{
  int i;

  DUPZfree(THIS->factor);
  DUPZfree(THIS->quot);
  DUPZfree(THIS->rem);
/*****for (i=0; i < THIS->nfactors_orig; i++) DUPZfree(THIS->factors[i]);*****/
  FREE(THIS->factors);
  mpz_clear(THIS->modulus);
  mpz_clear(THIS->f_value);
  mpz_clear(THIS->value_at_1);
  mpz_clear(THIS->value_at_m1);
  mpz_clear(THIS->P);
  mpz_clear(THIS->Q);
  for(i=0; i < THIS->nfactors_orig; i++) mpz_clear(THIS->factor_value[i]);
  FREE(THIS->factor_value);
  FREE(THIS->used);
  FREE(THIS->n1coeffd);
  FREE(THIS->n1coeffd_sum);

  FREE(THIS->combination);
  FREE(THIS->tc_combination);
  for (i=0; i < THIS->nfactors_orig; i++) mpz_clear(THIS->tc_product[i]);
  FREE(THIS->tc_product);

  FREE(THIS);
}



/***************************************************************************/
/* These are not constructor/destructor; "_init" must be called before     */
/* searching combinations, and "_done" must be called after all searching  */
/* has been finished.  These functions deal with globals.                  */

void DUPZfactor_combine_init()
{
  divisions=count1=count2=count3=count4=count5=count6=0;
  mpz_init(num);
  mpz_init(den);
  mpz_init(tmp);
  mpz_init(lcm);
}

void DUPZfactor_combine_done()
{
  mpz_clear(num); mpz_clear(den);  mpz_clear(tmp); mpz_clear(lcm);
#ifdef FACTOR_DEBUG
  printf("\nSearch complete.  count1=%d, count2=%d, count3=%d, count4=%d, count5=%d, count6=%d\n", count1, count2, count3, count4, count5, count6);
  printf("Total number of divisions: %d\n", divisions);
#endif
}


/***************************************************************************/

static int DUPZfactor_combine1(DUPZfactor_combiner THIS, int n, int r, int i);

static void DUPZfactor_combine(DUPZfactor_combiner THIS)
{
  int tuple_size, max_tuple_size;
  double log_modulus, log_lc, factor_bound;

  THIS->single_factor = 0;
  log_modulus = mpz_log(THIS->modulus);
  if (THIS->early)
  {
    max_tuple_size = (int)(10/logi(THIS->nfactors));
    if (max_tuple_size > THIS->nfactors/2) max_tuple_size = THIS->nfactors/2;
  }
  else
  {
    log_lc = THIS->info->bounds->leading_coeff;
    factor_bound = THIS->info->bounds->all_factors;
    /* Search through only half the combinations if modulus is big enough */
    if (log_modulus > 2*log_lc + factor_bound + 0.7) max_tuple_size = THIS->nfactors/2;
    else { THIS->single_factor = 1; max_tuple_size = THIS->nfactors - 1; }
  }
#ifdef FACTOR_DEBUG
  printf("max_tuple_size=%d\n", max_tuple_size);
#endif
  for (tuple_size=1; tuple_size <= max_tuple_size; tuple_size++)
  {
    if (tuple_size >= THIS->nfactors) break; /* only relevant during early searching */
#ifdef FACTOR_DEBUG
    printf("\nInterim.  count1=%d, count2=%d, count3=%d, count4=%d, count5=%d, count6=%d", count1, count2, count3, count4, count5, count6);
    printf("\n%d-tuples: ", tuple_size);
#endif
    THIS->tuple_size = tuple_size;
    THIS->tc_combination[tuple_size-1] = -1; /* so first call to tctest works properly */
    DUPZfactor_combine1(THIS, THIS->nfactors, 0, 0);
    if (THIS->early) continue;
    /* Check if we can reduce to a search through half the combinations. */
    log_lc = THIS->info->bounds->leading_coeff;
    factor_bound = THIS->info->bounds->all_factors;
    THIS->single_factor = 0;
//   2020-04-25 See redmine 1449.  This clever trick seems to be "bogus"  :-(
//    if (log_modulus > 2*log_lc + factor_bound + 0.7)
//    { max_tuple_size = THIS->nfactors/2; }
//    else { THIS->single_factor = 1; max_tuple_size = THIS->nfactors - 1; }
    THIS->single_factor = 1; max_tuple_size = THIS->nfactors - 1;
  }
  if (THIS-> early || DUPZdeg(THIS->info->f) <= 0) return;
  DUPZfactors_add(THIS->info->irreds, THIS->info->f);
  THIS->info->f->deg = 0;
}


void DUPZfactor_combine_early(DUPZfactor_combiner THIS)
{
  THIS->early = 1; /* flag that we're doing early searching */
  DUPZfactor_combine(THIS);
}


void DUPZfactor_combine_final(DUPZfactor_combiner THIS)
{
  THIS->early = 0; /* flag that this is the final search */
  DUPZfactor_combine(THIS);
}



/***************************************************************************/
/***************************************************************************/
/* n-1 coeff test */

static int DUPZfactor_combine_n1test(DUPZfactor_combiner THIS)
{
  int i, deg, df, D;
  const DUPZ f = THIS->info->f;
  DUPZ fi; /* just a convenient alias */
  double n1_coeff_bound;

  deg = THIS->tuple_deg;

  df = DUPZdeg(f);

  /* Compute in tmp the n-1 coeff of the current tuple (or its complement).    */
  /* Cannot just subtract from n-1 coeff in f because that value isn't to hand. */
  mpz_set_ui(tmp, 0);
  if (THIS->single_factor || deg <= df/2) /* tuple has low degree, so compute n-1 coeff directly */
  {
    D = deg;
    for (i=0; i < THIS->tuple_size; i++)
    {
      fi = *THIS->factors[THIS->combination[i]];
      mpz_add(tmp, tmp, fi->coeffs[DUPZdeg(fi)-1]);
    }
  }
  else /* tuple has high degree, so compute n-1 coeff of its complement */
  {
    D = df - deg;
    for (i=0; i < THIS->nfactors_orig; i++)
    {
      if (THIS->used[i]) continue;
      fi = *THIS->factors[i];
      mpz_add(tmp, tmp, fi->coeffs[DUPZdeg(fi)-1]);
    }
  }
  mpz_fdiv_r(tmp, tmp, THIS->modulus);
  /* Now compute the rational corresponding to the modular value. */
  modular_to_rational(num, den, THIS->Q, tmp, THIS->modulus);

  /* check magnitude of numerator */
  if (mpz_sgn(num) == 0) return 1;
  n1_coeff_bound = THIS->info->bounds->root_bound + logi(D) + 0.001;
  if (mpz_log(num) - mpz_log(den) > n1_coeff_bound) return 0;

  /* check denominator divides leading coeff */
  mpz_fdiv_r(num, DUPZlc(f), den);
  if (mpz_sgn(num) != 0) return 0;
  return 1;
}


/***************************************************************************/
/* Trailing coeff test                                                     */
/* SIDE EFFECT: sets THIS->log_tc (used in coefftest below)                */

static int DUPZfactor_combine_tctest(DUPZfactor_combiner THIS)
{
  DUPZ f = THIS->info->f;  /* convenient alias */
  int i, deg, df;
  DUPZ fi; /* just a convenient alias */

  deg = THIS->tuple_deg;
  df = DUPZdeg(f);
  /* compute in tmp the trailing coeff of the current combination */
  mpz_set_ui(tmp, 1);
  if (THIS->single_factor || deg <= df/2)
  {
    for (i=0; i < THIS->tuple_size; i++)
    {
      if (THIS->combination[i] == THIS->tc_combination[i]) continue;
      THIS->tc_combination[i] = THIS->combination[i];
      THIS->tc_combination[i+1] = -1;
      fi = *THIS->factors[THIS->combination[i]];
      if (i == 0) mpz_set(THIS->tc_product[i], fi->coeffs[0]);
      else
      {
        mpz_mul(THIS->tc_product[i], THIS->tc_product[i-1], fi->coeffs[0]);
        mpz_fdiv_r(THIS->tc_product[i], THIS->tc_product[i], THIS->modulus);
      }
    }
    mpz_set(tmp, THIS->tc_product[i-1]);
  }
  else
  {
    for (i=0; i < THIS->nfactors_orig; i++)
    {
      if (THIS->used[i]) continue;
      fi = *THIS->factors[i];
      mpz_mul(tmp, tmp, fi->coeffs[0]);
      mpz_fdiv_r(tmp, tmp, THIS->modulus);
    }
  }

  /* NB tmp != 0 because we excluded primes dividing the trailing coeff */
  modular_to_rational(num, den, THIS->Q, tmp, THIS->modulus);

  /* check numerator divides trailing coeff of original polynomial */
  mpz_fdiv_r(tmp, f->coeffs[0], num);
  if (mpz_sgn(tmp) != 0) return 0;

  /* check denominator divides leading coeff */
  mpz_fdiv_r(tmp, DUPZlc(f), den);
  if (mpz_sgn(tmp) != 0) return 0;

  THIS->log_tc = mpz_log(num) - mpz_log(den); /* used in coefftest below */
  return 1;
}


static int DUPZfactor_combine_evaltest(DUPZfactor_combiner THIS)
{
  int i, deg, df, D;
  DUPZ f = THIS->info->f; /* convenient alias */

  df = DUPZdeg(f);
  deg = THIS->tuple_deg;
  if (THIS->single_factor || deg <= df/2) D = deg; else D = df - deg;
  /* Do eval at r test only if P/Q > (r + root_bound)^D */
  if (mpz_log(THIS->P) - mpz_log(THIS->Q) < D*THIS->log_ev_pt + 0.001) return 1;
  /* compute in tmp the value at r of the current combination (or its complement) */
  mpz_set_ui(tmp, 1);
  if (THIS->single_factor || deg <= df/2) /* tuple has low degree, so compute its value directly */
  {
    for (i=0; i < THIS->tuple_size; i++)
    {
      mpz_mul(tmp, tmp, THIS->factor_value[THIS->combination[i]]);
      mpz_fdiv_r(tmp, tmp, THIS->modulus);
    }
  }
  else /* tuple has high degree, so compute value of its complement instead */
  {
    for (i=0; i < THIS->nfactors_orig; i++)
    {
      if (THIS->used[i]) continue;
      mpz_mul(tmp, tmp, THIS->factor_value[i]);
      mpz_fdiv_r(tmp, tmp, THIS->modulus);
    }
  }

  modular_to_rational(num, den, THIS->Q, tmp, THIS->modulus);

  if (mpz_sgn(num) == 0) return 0;

  /* check numerator divides the value of f */
  mpz_fdiv_r(num, THIS->f_value, num);
  if (mpz_sgn(num) != 0) return 0;

  /* check denominator divides lc(f) */
  mpz_fdiv_r(den, DUPZlc(f), den);
  if (mpz_sgn(den) != 0) return 0;

  return 1;
}


/***************************************************************************/
/* Coefficient test: the modular factors are multiplied, and the coeffs of */
/* the product converted to rationals checking numerator and denominator   */
/* have possible values.  This test is relatively expensive.               */

static int DUPZfactor_combine_coefftest(DUPZfactor_combiner THIS)
{
  int i, j, deg, df, D;
  DUPZ f = THIS->info->f;             /* alias */
  DUPZ g = THIS->factor;              /* alias */
  double log_coeff, forward, reverse; /* forward and reverse bounds */

  deg = THIS->tuple_deg;
  df = DUPZdeg(f);
  /* if deg <= df/2 then we multiply out the combination,     */
  /* otherwise we multiply out the complementary combination. */
  if (THIS->single_factor || deg <= df/2)
  {
    D = deg;
    DUPZcopy2(g, *THIS->factors[THIS->combination[0]]);
    for (i=1; i < THIS->tuple_size; i++)
    {
      DUPZmul3(g, g, *THIS->factors[THIS->combination[i]]);
      for (j=DUPZdeg(g); j >= 0; j--)
	mpz_fdiv_r(g->coeffs[j], g->coeffs[j], THIS->modulus);
    }
  }
  else
  {
    D = df - deg;
    mpz_set_ui(g->coeffs[0], 1);
    g->deg = 0;
    for (i=0; i < THIS->nfactors_orig; i++)
    {
      if (THIS->used[i]) continue;
      DUPZmul3(g, g, *THIS->factors[i]);
      for (j=DUPZdeg(g); j >= 0; j--)
	mpz_fdiv_r(g->coeffs[j], g->coeffs[j], THIS->modulus);
    }
  }


  forward = 0.001;          /* add 0.001 to allow for rounding error */
  reverse = THIS->log_tc + D*THIS->info->bounds->reverse_root_bound + 0.001;
  mpz_set_ui(lcm, 1);
  for (i=D-1; i >= 0; i--)
  {
    forward += logi(i+1) - logi(D-i) + THIS->info->bounds->root_bound;
    reverse += logi(i+1) - logi(D-i) - THIS->info->bounds->reverse_root_bound;
    if (mpz_sgn(g->coeffs[i]) == 0) continue;
    modular_to_rational(num, den, THIS->Q, g->coeffs[i], THIS->modulus);
    mpz_fdiv_r(tmp, DUPZlc(f), den);          /* check denom divides lc(f) */
    if (mpz_sgn(tmp) != 0) return 0;          /* if not, factor is false */

    /* check coeff is within permitted range */
    log_coeff = mpz_log(num) - mpz_log(den);
    if (log_coeff > forward || log_coeff > reverse) return 0;

    /* store rational in THIS->factor and THIS->quot */
    mpz_set(g->coeffs[i], num);
    mpz_set(THIS->quot->coeffs[i], den);

    /* update the common denominator */
    mpz_gcd(tmp, lcm, den);
    mpz_divexact(den, den, tmp);
    mpz_mul(lcm, lcm, den);
  }
  
  /* make THIS->factor have a common denominator */
  for (i=D-1; i >= 0; i--)
  {
    if (mpz_sgn(g->coeffs[i]) == 0) continue;
    if (mpz_cmp(lcm, THIS->quot->coeffs[i]) == 0) continue;
    mpz_divexact(den, lcm, THIS->quot->coeffs[i]);
    mpz_mul(g->coeffs[i], g->coeffs[i], den);
  }
  mpz_set(g->coeffs[D], lcm);

  return 1;
}


/***************************************************************************/
/* Division test: we divide the putative factor into the polynomial.       */
/* We initially check that the values at 1 and -1 divide those of f.       */
/* Then we just do long division.  This test is potentially very expensive */
/* if the putative factor is not a true one -- should hardly ever happen.  */

static int DUPZfactor_combine_dividetest(DUPZfactor_combiner THIS)
{
  DUPZ f = THIS->info->f;  /* convenient alias */

  /* These two preliminary checks are most relevant during early detection */
  /* check value of factor at 1 divides the value of f at 1 */
  DUPZevaluate(tmp, THIS->factor, 1);
  if (mpz_sgn(tmp) == 0) return 0;
  mpz_fdiv_r(tmp, THIS->value_at_1, tmp);
  if (mpz_sgn(tmp) != 0) return 0;

  /* check value of factor at -1 divides the value of f at -1 */
  DUPZevaluate(tmp, THIS->factor, -1);
  if (mpz_sgn(tmp) == 0) return 0;
  mpz_fdiv_r(tmp, THIS->value_at_m1, tmp);
  if (mpz_sgn(tmp) != 0) return 0;

  /* we expect this division almost always to be exact */
  DUPZdiv4(THIS->quot, THIS->rem, f, THIS->factor);
  divisions++;
  if (DUPZdeg(THIS->rem) >= 0) return 0;

  /* Swap factor and quot if we have built the complement in factor. */
  if (!THIS->single_factor && THIS->tuple_deg > DUPZdeg(THIS->info->f)/2)
    DUPZswap(THIS->quot, THIS->factor);
  return 1;
}


/***************************************************************************/
/* This function is called if we cannot immediately show that the factor   */
/* is irreducible.  A more stringent irreducibility test is conducted, and */
/* if necessary, a subsidiary factorization is performed.                  */

static void DUPZfactor_combine_fork(DUPZfactor_combiner old, DUPZ factor)
{
  int i, irred, dmax;
  DUPZfactor_info fork;
  DUPFF fp, rem;
  DUPFFlist iter;
  bound_info bd;
  double lift_height = mpz_log(old->modulus) - 0.001; /* allow for rounding */

  /* initial quick tests for "obvious" irreducibility */
  if (old->tuple_size == 1 || !old->early)
  {
    DUPZfactors_add(old->info->irreds, factor);
    return;
  }
  dmax = DUPZdeg(factor)/2;
  while (!old->info->fds[dmax]) dmax--;
  bd = DUPZcopy_refine_bound(old->info->bounds, factor, dmax);
  if (bd->lift_bound < lift_height) 
  {
    FREE(bd);
    DUPZfactors_add(old->info->irreds, factor);
    return;
  }
  fork = (DUPZfactor_info)MALLOC(sizeof(struct DUPZfactor_info_struct));
  fork->f = DUPZcopy(factor);
  fork->irreds = old->info->irreds;
  fork->p = old->info->p;

  fork->nprimes = old->info->nprimes;
  fork->FFq = (FF*)MALLOC(fork->nprimes*sizeof(FF));
  fork->qfactors = (DUPFFlist*)MALLOC(fork->nprimes*sizeof(DUPFFlist));
  fork->fds = (int*)MALLOC((1+DUPZdeg(fork->f))*sizeof(int));
  for (i=0; i <= DUPZdeg(fork->f); i++) fork->fds[i] = old->info->fds[i];
  fork->bounds = NULL;
  mpz_init_set(fork->recip_lcf, old->info->recip_lcf);
  mpz_mul(fork->recip_lcf, fork->recip_lcf, DUPZlc(old->info->f));
  mpz_divexact(fork->recip_lcf, fork->recip_lcf, DUPZlc(factor));
  mpz_init_set(fork->Q, old->info->Q);
  fork->lifter = NULL;

  /* determine new factorizations modulo the various primes */
  fp = DUPFFnew(DUPZdeg(fork->f));
  rem = DUPFFnew(DUPZdeg(fork->f));
  irred = 0;
  for (i=0; i < old->info->nprimes; i++)
  {
    FFselect(old->info->FFq[i]);
    fork->FFq[i] = FFctor(CurrentFF.prime); /* wasteful copy */
    fork->qfactors[i] = NULL;
    DUPZ_to_DUPFF2(fp, fork->f);
    for (iter = old->info->qfactors[i]; iter; iter = iter->next)
    {
      DUPFFcopy2(rem, fp);
      DUPFFrem2(rem, iter->poly);
      if (DUPFFdeg(rem) >= 0) continue;
      fork->qfactors[i] = DUPFFlist_append(DUPFFlist_ctor(DUPFFcopy(iter->poly), 1), fork->qfactors[i]);
    }
    irred = DUPZfactor_refine_fds(fork->fds, fork->qfactors[i]);
/*    if (irred) break;*/ /* including this line makes deallocation trickier */
  }
  /* select correct prime for rest of computation */
  fork->p_index = old->info->p_index;
  fork->pfactors = fork->qfactors[fork->p_index];
  FFselect(fork->FFq[fork->p_index]);
  DUPFFfree(fp);
  DUPFFfree(rem);
  if (irred) 
  {
    FREE(bd);
    DUPZfactors_add(fork->irreds, fork->f);
    DUPZfactor_info_dtor(fork);
    return;
  }

  /* update bounds */
  fork->dmax = old->info->dmax;
  if (fork->dmax > DUPZdeg(fork->f)/2) fork->dmax = DUPZdeg(fork->f)/2;
  
  while (!fork->fds[dmax]) dmax--;
  fork->dmax = dmax;
/*  fork->dmin = old->info->dmin;                */
/*  while (!fork->fds[fork->dmin]) fork->dmin++; */
  DUPZrefine_bound(bd, fork->f, dmax);
  fork->bounds = bd;
  if (bd->lift_bound < lift_height)
  {
    DUPZfactors_add(fork->irreds, fork->f);
    DUPZfactor_info_dtor(fork);
    return;
  }
  DUPZfactor_info_set_target_height(fork);


/*************************************************************  DUPZfactor_lift_revise1(fork, old->info->lifter->r); **********************/
#ifdef FACTOR_DEBUG
  printf("\nFactor may be reducible: starting subsidiary factorization\n");
#endif
  DUPZfactor_finish(fork);
#ifdef FACTOR_DEBUG
  printf("Subsidiary factorization complete; resuming previous factorization\n");
#endif
  DUPZfactor_info_dtor(fork);
}

/***************************************************************************/
/* Upon entry to this function THIS->factor contains the factor built from */
/* the current tuple, while THIS->quot contains its complement (i.e. the   */
/* quotient of f by the factor found.                                      */

static void DUPZfactor_combine_found(DUPZfactor_combiner THIS)
{
  const DUPZ f = THIS->info->f;  /* convenient alias */
DUPZ fi;
  double rt_bd;
  int i, irred, df;
  DUPFF rem, fp;
  DUPFFlist iter, next, *prev;

  DUPZfactor_combine_fork(THIS, THIS->factor); /* recurse on factor, needed only if using early detection */
for(i=0;i<THIS->tuple_size;i++){fi=*THIS->factors[THIS->combination[i]];fi->deg=-fi->deg;}
  DUPZcopy2(f, THIS->quot);
  mpz_mul(THIS->info->recip_lcf, THIS->info->recip_lcf, DUPZlc(THIS->factor));
  df = DUPZdeg(f);
  THIS->nfactors -= THIS->tuple_size;

  /* ---- we may be able to prove that what's left is irreducible ----- */
  /* determine new factorizations modulo the various primes */
#ifdef FACTOR_DEBUG
  printf("Checking if quotient is irred\n");
#endif
  fp = DUPFFnew(df);
  rem = DUPFFnew(df);
  irred = 0;
  for (i=0; i < THIS->info->nprimes; i++)
  {
    FFselect(THIS->info->FFq[i]);
    prev = &THIS->info->qfactors[i];
    DUPZ_to_DUPFF2(fp, f);
    for (iter = THIS->info->qfactors[i]; iter; iter = next)
    {
      next = iter->next;
      DUPFFcopy2(rem, fp);
      DUPFFrem2(rem, iter->poly);
      if (DUPFFdeg(rem) >= 0) DUPFFlist_elem_dtor(iter);
      else { *prev = iter; prev = &iter->next; }
    }
    *prev = NULL;
    irred = DUPZfactor_refine_fds(THIS->info->fds, THIS->info->qfactors[i]);
#ifdef FACTOR_DEBUG
    printf("irred (via FDS) for remaining factor: %d\n",irred);
#endif
    if (irred) break;
  }

  THIS->info->pfactors = THIS->info->qfactors[THIS->info->p_index];
  FFselect(THIS->info->FFq[THIS->info->p_index]);
  DUPFFfree(fp);
  DUPFFfree(rem);
  if (irred) 
  {
    DUPZfactors_add(THIS->info->irreds, f);
    for (i=0; i < THIS->nfactors_orig; i++) THIS->used[i] = 1;
for (i=0; i < THIS->nfactors_orig; i++) {fi=*THIS->factors[i];if (DUPZdeg(fi)>0) fi->deg = -fi->deg; }
    THIS->nfactors = 0;
    f->deg = 0;
    return;
  }
  if (THIS->info->dmax > df/2) THIS->info->dmax = df/2;
  while (!THIS->info->fds[THIS->info->dmax]) THIS->info->dmax--;
/*  while (!THIS->info->fds[THIS->info->dmin]) THIS->info->dmin++; */

  /* update bounds */
#ifdef FACTOR_DEBUG
  printf("Updating factor bounds...\n");
#endif
  DUPZrefine_bound(THIS->info->bounds, f, THIS->info->dmax);
  DUPZfactor_info_set_target_height(THIS->info);
/* #ifdef FACTOR_DEBUG */
/* printf("bound info\n----------\n"); */
/* printf("deg\t\t%d\n", THIS->info->bounds->df); */
/* printf("lc\t\t%f\n", THIS->info->bounds->leading_coeff); */
/* printf("tc\t\t%f\n", THIS->info->bounds->trailing_coeff); */
/* printf("dmax\t\t%d\n", THIS->info->bounds->dmax); */
/* printf("Bombieri2\t%f\n", THIS->info->bounds->Bombieri2); */
/* printf("root_bound\t%f\n", THIS->info->bounds->root_bound); */
/* printf("rev_root_bound\t%f\n", THIS->info->bounds->reverse_root_bound); */
/* printf("all_factors\t%f\n", THIS->info->bounds->all_factors); */
/* printf("lift_bound\t%f\n", THIS->info->bounds->lift_bound); */
/* #endif */

  if (THIS->early && THIS->info->target_height <= THIS->info->current_height)
  {
    THIS->early = 0;
#ifdef FACTOR_DEBUG
    printf("[early search has now become final search]");
#endif
  }

  /* update log_ev_pt */
  rt_bd = THIS->info->bounds->root_bound;
  THIS->log_ev_pt = rt_bd + log(1+exp(logi(abs(THIS->evaluation_point)) - rt_bd));

  /* update f_value */
  DUPZevaluate(THIS->f_value, f, THIS->evaluation_point);

  /* update P, Q */
  DUPZfactor_combiner_set_PQ(THIS);

  /* update n1coeff_bounds */
  DUPZfactor_combiner_set_n1coeffd_bounds(THIS);
}



static int DUPZfactor_combine1(DUPZfactor_combiner THIS, int n, int r, int i)
/* any true factors found will be added to ans.   */
/* n is number of factors still to be considered, */
/* r is number of factors already picked,         */
/* i is index of next factor to consider          */
{
  int deg, found;

  if (n < THIS->tuple_size - r) return 0;
  if (r < THIS->tuple_size)
  {
    /*    if (deg > dmax) return;*/
    if (THIS->used[i]) return DUPZfactor_combine1(THIS, n, r, i+1);
    THIS->combination[r] = i;
    THIS->used[i] = 2;
    deg = DUPZdeg(*THIS->factors[i]);
    THIS->tuple_deg += deg;
    THIS->n1coeffd_sum[r+1] = THIS->n1coeffd_sum[r] + THIS->n1coeffd[i];
     /* first try with this factor */
    found = DUPZfactor_combine1(THIS, n-1, r+1, i+1);
    THIS->tuple_deg -= deg;
    if (!found) THIS->used[i] = 0;
    /* special case if n = nfactors/2 */
    if (!found && !THIS->single_factor && (r == 0) && (THIS->nfactors == 2*THIS->tuple_size)) return 0;
    if (!found) return DUPZfactor_combine1(THIS, n-1, r, i+1);
    if (r > 0) return 1; /* unwind to "top level" */
    /* now we must be at "top level" */
    if (!THIS->single_factor && THIS->nfactors < 2*THIS->tuple_size) return 0;
    if (THIS->nfactors <= THIS->tuple_size) return 0;

    /* now try without this factor */
    return DUPZfactor_combine1(THIS, n - THIS->tuple_size, r, i+1);
  }

#ifdef FACTOR_DEBUG
  printf("Trying tuple: \t"); for (i=0; i < THIS->tuple_size; i++) printf("%d ", THIS->combination[i]); printf("\n");
#endif
count1++;
  
  /* check degree is possible */
  if (THIS->info->fds[THIS->tuple_deg] == 0) return 0;
count2++;


  if (THIS->n1coeffd_flag)
  {
    /* Compute in deg the nearest integer to THIS->n1coeffd_sum */
    if (THIS->n1coeffd_sum[THIS->tuple_size] >= 0)
      deg = (int)(THIS->n1coeffd_sum[THIS->tuple_size] + 0.5);
    else deg = (int)(THIS->n1coeffd_sum[THIS->tuple_size] - 0.5);

    if (fabs(THIS->n1coeffd_sum[THIS->tuple_size] - deg) >= THIS->tuple_deg * THIS->n1coeffd_bound) return 0;
  }
  else
    if (!DUPZfactor_combine_n1test(THIS)) return 0;

/*
 *  if (!DUPZfactor_combine_1test(THIS)) {printf("!");fflush(stdout);return 0;}
 */
count3++;

  if (!DUPZfactor_combine_tctest(THIS)) return 0;
count4++;

  if (!DUPZfactor_combine_evaltest(THIS)) return 0;
count5++;

  if (!DUPZfactor_combine_coefftest(THIS)) return 0;
count6++;

  if (!DUPZfactor_combine_dividetest(THIS)) return 0;

  DUPZfactor_combine_found(THIS);

  return 1;
}
