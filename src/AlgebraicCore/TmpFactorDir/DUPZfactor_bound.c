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


#include "DUPZfactor_bound.h"
#include "DUPZroot_bound.h"
#include "mpz_log.h"
#include "logi.h"
#include "jalloc.h"
#include "add_logs.h"


static double log_choose(int n, int r)
{
  double ans;
  int i;

  ans = 0.0;
  if (2*r > n) r = n-r;
  for (i=1; i <= r; i++)
    ans += logi(n+1-i) - logi(i);
  return ans;
}

/***************************************************************************/
/* This returns an upper bound for the log of the Mahler measure of f.       */
/* M(f) is the product of max(1, abs(alpha_i)) where alpha_i are roots of f. */
/* An easy upper bound is the standard 2-norm of f.                          */

static double Mahler(const DUPZ f)
{
  double lcf, norm2;
  int i;

  if (DUPZdeg(f) < 0) return -1;
  lcf = mpz_log(DUPZlc(f));
  norm2 = 2*lcf;
  for (i=0; i < DUPZdeg(f); i++)
    if (mpz_sgn(f->coeffs[i]) != 0)
      norm2 = add_logs(norm2, 2*mpz_log(f->coeffs[i]));

  return norm2/2 - lcf;
}

static double Mignotte_bound(const bound_info bd)
{
  return bd->dmax*logi(2) + bd->Mahler + bd->leading_coeff;
}

/***************************************************************************/

static double Bombieri_norm(const DUPZ f)
{
  int i, df;
  double Bombieri, log_binomial;

  df = DUPZdeg(f);
  /* we have to work with logs to avoid overflow problems */
  Bombieri = 2*mpz_log(f->coeffs[0]);
  log_binomial = 0.0;
  for (i=1; i <= df; i++)
  {
    log_binomial += logi(df+1-i) - logi(i);
    if (mpz_sgn(f->coeffs[i]) != 0)
      Bombieri = add_logs(Bombieri, 2*mpz_log(f->coeffs[i])-log_binomial);
  }
  return Bombieri/2;
}


static double Beauzamy_bound(const bound_info bd)
{
  double bit_inside_sqrt;

  bit_inside_sqrt = log_choose(bd->df, bd->dmax) + log_choose(bd->dmax, bd->dmax/2) - logi(2);
  return bit_inside_sqrt/2 + bd->Bombieri2;
}


static double JAA_bound(const bound_info bd)
{
  int i;
  double coeff_bound, *coeffs, *revcoeffs;
  const int dmax = bd->dmax;

  /* the next section is a pretty "stupid" implementation but should work */
  coeffs = (double*)MALLOC((dmax+1)*sizeof(double));
  revcoeffs = (double*)MALLOC((dmax+1)*sizeof(double));
  coeffs[dmax] = bd->leading_coeff;
  for (i=dmax-1; i >= 0; i--)
  {
    coeffs[i] = coeffs[i+1] + bd->root_bound + logi(i+1) - logi(dmax-i);
  }
  revcoeffs[0] = bd->trailing_coeff;
  for (i=1; i <= dmax; i++)
  {
    revcoeffs[i] = revcoeffs[i-1] + bd->reverse_root_bound + logi(dmax-i+1) - logi(i);
  }
  /* The bound is max{ min(coeffs[i], revcoeffs[i]): i=0,1,..,dmax } */
  coeff_bound = 0.0;
  for (i=0; i <= dmax; i++)
  {
    if (coeffs[i] > revcoeffs[i])
    {
      if (coeff_bound < revcoeffs[i]) coeff_bound = revcoeffs[i];
    }
    else
      if (coeff_bound < coeffs[i]) coeff_bound = coeffs[i];
  }
  FREE(coeffs);
  FREE(revcoeffs);
  return coeff_bound;
}


static double single_factor_bound(const bound_info bd)
{
  double single_factor;

  single_factor = bd->Bombieri2/2 + 0.3466*bd->df + 0.25/bd->df - 0.375*logi(bd->df) + 0.004;
  if (bd->df == 2) single_factor = bd->Bombieri2/2 + logi(2)/4;
  return single_factor;
}


static void set_lift_bound(bound_info bd)
{
  double Beauzamy, JAA, single_factor, Mignotte;

  Beauzamy = Beauzamy_bound(bd);
  JAA = JAA_bound(bd);
  Mignotte = Mignotte_bound(bd);
  bd->all_factors = Beauzamy;
  if (JAA < bd->all_factors) bd->all_factors = JAA;
  if (Mignotte < bd->all_factors) bd->all_factors = Mignotte;
  single_factor = single_factor_bound(bd);
  /* Use single factor bound only if it is "much" better; 1.4 is a guess. */
  if (bd->all_factors < 1.4 * single_factor) bd->lift_bound = bd->all_factors;
  else bd->lift_bound = single_factor;
bd->lift_bound = bd->all_factors; /* DISABLE SINGLE FACTOR BOUND COMPLETELY */
#ifdef FACTOR_DEBUG
printf("Beauzamy bound is %f\n", Beauzamy);
printf("JAA bound is %f\n", JAA);
printf("Mignotte bound is %f\n", Mignotte);
printf("Single factor bound is %f\n", single_factor);
#endif
}



bound_info DUPZfactor_bound_ctor(const DUPZ f, int dmax)
{
  DUPZ revf;
  bound_info bd = (bound_info)MALLOC(sizeof(struct DUPZfactor_bound_struct));

  bd->df = DUPZdeg(f);
  bd->Bombieri2 = Bombieri_norm(f);
  bd->Mahler = Mahler(f);

  bd->leading_coeff = mpz_log(DUPZlc(f));
  bd->trailing_coeff = mpz_log(f->coeffs[0]);
  bd->dmax = dmax;
  bd->root_bound = DUPZroot_bound(f);
  revf = DUPZcopy(f);
  DUPZreverse(revf);
  bd->reverse_root_bound = DUPZroot_bound(revf);
  DUPZfree(revf);

  set_lift_bound(bd);
#ifdef FACTOR_DEBUG
printf("bound info\n----------\n");
printf("deg\t\t%d\n", bd->df);
printf("lc\t\t%f\n", bd->leading_coeff);
printf("tc\t\t%f\n", bd->trailing_coeff);
printf("dmax\t\t%d\n", bd->dmax);
printf("Bombieri2\t%f\n", bd->Bombieri2);
printf("root_bound\t%f\n", bd->root_bound);
printf("rev_root_bound\t%f\n", bd->reverse_root_bound);
printf("all_factors\t%f\n", bd->all_factors);
printf("lift_bound\t%f\n", bd->lift_bound);
#endif

  return bd;
}


void DUPZfactor_bound_dtor(bound_info bd)
{
  FREE(bd);
}


/***************************************************************************/
/* It is assumed that f here is a factor of the polynomial used originally */
/* to create bd; if not, the result could be junk data.                    */

void DUPZrefine_bound(bound_info bd, const DUPZ f, int dmax)
{
  double tmp;
  DUPZ revf;

  bd->df = DUPZdeg(f);
  bd->leading_coeff = mpz_log(DUPZlc(f));
  bd->trailing_coeff = mpz_log(f->coeffs[0]);
  bd->Bombieri2 = Bombieri_norm(f);
  bd->Mahler = Mahler(f);

  if (dmax < bd->dmax) bd->dmax = dmax;
  tmp = DUPZroot_bound(f);
  if (tmp < bd->root_bound) bd->root_bound = tmp;
  revf = DUPZcopy(f);  DUPZreverse(revf);
  tmp = DUPZroot_bound(revf);
  DUPZfree(revf);
  if (tmp < bd->reverse_root_bound) bd->reverse_root_bound = tmp;
  set_lift_bound(bd);
}


bound_info DUPZcopy_refine_bound(bound_info bd, const DUPZ f, int dmax)
{
  bound_info ans;

  ans = (bound_info)MALLOC(sizeof(struct DUPZfactor_bound_struct));
  ans->root_bound = bd->root_bound;
  ans->reverse_root_bound = bd->reverse_root_bound;
  ans->dmax = bd->dmax;
  DUPZrefine_bound(ans, f, dmax);
  return ans;
}
