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
//#include <math.h>
#include "DUPZbivariate_lift.h"
#include "jalloc.h"
#include "jaaerror.h"
#include "DUPZcontent.h"
#include "DUPZexgcd.h"
#include "DUPZfactors.h"
#include "primes.h"
#include "DUPFF.h"
#include "DUPZ_DUPFF.h"
#include "DUPZcra.h"

static int DUPFFbivariate_lift_single_factor(DUPFF *factor, const DUPFF *f, int degy)
{
  DUPFF cofac, other_cofac, *other_factor, tmp, tmp2, gcd, corr, other_corr, error, zero;
  int degx, i, level, j;
/*  int lcfactor;*/ /* Not used any more*/

printf("*");
//#ifdef FACTOR_DEBUG
//printf("modular univariate factor is ");DUPFFprint(factor[0]);
//#endif
  zero = DUPFFnew(-1);
  degx = DUPFFdeg(f[0]);
  other_factor = (DUPFF*)MALLOC((1+degy)*sizeof(DUPFF));
//  for (i=0; i <= degy; i++) other_factor[i] = DUPFFnew(degx - DUPFFdeg(factor[0]));
  other_factor[0] = DUPFFnew(degx); /* too big, should be degx - deg(factor[0])*/
  DUPFFdiv4(other_factor[0], zero, f[0], factor[0]);
//DUPFFmake_monic(factor[0]);
//DUPFFmake_monic(other_factor[0]);

  cofac = DUPFFnew(DUPFFdeg(other_factor[0])-1);
  other_cofac = DUPFFnew(DUPFFdeg(factor[0])-1);
  gcd = DUPFFexgcd(&cofac, &other_cofac, factor[0], other_factor[0]);
  DUPFFdiv2ff(cofac, DUPFFlc(gcd));
  DUPFFdiv2ff(other_cofac, DUPFFlc(gcd));
  if (DUPFFdeg(gcd) != 0) return 0; /* tidy memory */

  corr = DUPFFnew(DUPFFdeg(f[0]) + DUPFFdeg(other_cofac));
  other_corr = DUPFFnew(DUPFFdeg(f[0]) + DUPFFdeg(cofac));
  tmp = DUPFFnew(DUPFFdeg(f[0]));
  tmp2 = DUPFFnew(DUPFFdeg(f[0]));
  for (level=1; level <= degy; level++)
  {
    tmp->deg = -1;
    for (j=1; j < level; j++)
    {
      DUPFFmul3(tmp2, factor[j], other_factor[level-j]);
      DUPFFadd3(tmp, tmp, tmp2);
    }
    error = DUPFFcopy(f[level]);
//    DUPFFdiv2ff(error, DUPFFlc(f[0])); /* simulate f being monic */
    error = DUPFFsub(error, tmp);
    DUPFFmul3(corr, error, other_cofac);
    factor[level] = DUPFFrem(corr, factor[0]);
//    DUPFFdiv2ff(factor[level], gcd->coeffs[0]);
    DUPFFmul3(other_corr, error, cofac);
    other_factor[level] = DUPFFrem(other_corr, other_factor[0]);
//    DUPFFdiv2ff(other_factor[level], gcd->coeffs[0]);
    DUPFFfree(error);
  }

  for (i=0; i <= degy; i++) DUPFFfree(other_factor[i]);
//  for (i=0; i <= degy; i++) DUPFFmul2ff(factor[i], DUPFFlc(factor[0])); /* force lc */
  FREE(other_factor);
  DUPFFfree(cofac); DUPFFfree(other_cofac);
  DUPFFfree(corr); DUPFFfree(other_corr);
  DUPFFfree(zero);
  DUPFFfree(gcd); DUPFFfree(tmp); DUPFFfree(tmp2);
  return 1;
}

void DUPZbivariate_lift_single_factor(DUPZ *factor, const DUPZ *f, int degy)
{
  int p, i, changed, OK;
  FF Fp;
  mpz_t modulus, content, content2;
  DUPFF *factorp, *fp;

//#ifdef FACTOR_DEBUG
//printf("DUPZbivariate_lift_single_factor called on:\n");
//printf("f is:\n");
//for (i=0; i<=degy; i++)
//{printf("component in y^%d = ", i); DUPZprint(f[i]);}
//printf("Univariate factor is: ");DUPZprint(factor[0]);
//#endif
  factorp = (DUPFF*)MALLOC((1+degy)*sizeof(DUPFF));
  fp = (DUPFF*)MALLOC((1+degy)*sizeof(DUPFF));
  for (i=0; i <= degy; i++)
  {
    factorp[i] = DUPFFnew(DUPZdeg(factor[0]));
    fp[i] = DUPFFnew(DUPZdeg(f[0]));
    if (i > 0) factor[i] = DUPZnew(DUPZdeg(factor[0]));
  }
  p = 31;
  mpz_init_set_ui(modulus, 1);
  changed = 1;
  while (changed)
  {
    do
    {
      p = nextprime(p);
    }
    while (mpz_fdiv_ui(DUPZlc(f[0]), p) == 0);

    Fp = FFctor(p);
    FFselect(Fp);
    for (i=0; i <= degy; i++)
      DUPZ_to_DUPFF2(fp[i], f[i]);
    DUPZ_to_DUPFF2(factorp[0], factor[0]);
    OK = DUPFFbivariate_lift_single_factor(factorp, fp, degy);
    if (!OK) continue; /* memory leak */
//#ifdef FACTOR_DEBUG
//printf("modular bivariate factor is:\n");
//for (i=0; i <= degy; i++)
//{printf("component in y^%d: ", i); DUPFFprint(factorp[i]);}
//#endif
    changed = 0;
    for (i=1; i <= degy; i++)
      changed |= DUPZcra(factor[i], modulus, factorp[i], p);
    FFdtor(Fp); /* finished with Fp now */
    mpz_mul_ui(modulus, modulus, p);
//#ifdef FACTOR_DEBUG
//printf("modulus=");mpz_out_str(stdout,10,modulus);printf("\n");
//printf("lifted factor so far is:\n");
//for (i=0; i <= degy; i++)
//{printf("component in y^%d: ", i); DUPZprint(factor[i]);}
//#endif
  }
  /* free up factorp and fp */
  mpz_clear(modulus);

  mpz_init(content);
  mpz_init(content2);
  for (i=0; i <= degy; i++)
  {
    DUPZcontent(content2, factor[i]);
    mpz_gcd(content, content, content2);
  }
  for (i=0; i <= degy; i++)
    DUPZdiv2z(factor[i], content);
  mpz_clear(content);
  mpz_clear(content2);
}

#if 0
void DUPZbivariate_lift_single_factor(DUPZ *factor, const DUPZ *f, int degy)
{
  DUPZ cofac, other_cofac, *other_factor, tmp, tmp2, gcd, corr, other_corr, error, zero;
  int degx, i, level, j;

  zero = DUPZnew(-1);
  degx = DUPZdeg(f[0]);
  other_factor = (DUPZ*)MALLOC((1+degy)*sizeof(DUPZ));
//  for (i=0; i <= degy; i++) other_factor[i] = DUPZnew(degx - DUPZdeg(factor[0]));
  other_factor[0] = DUPZnew(degx); /* too big, should be degx - deg(factor[0])*/
  DUPZdiv4(other_factor[0], zero, f[0], factor[0]);
  
  cofac = DUPZnew(DUPZdeg(other_factor[0])-1);
  other_cofac = DUPZnew(DUPZdeg(factor[0])-1);
  gcd = DUPZexgcd(cofac, other_cofac, factor[0], other_factor[0]);
  if (DUPZdeg(gcd) != 0) JERROR(JERROR_HENSEL);

  corr = DUPZnew(DUPZdeg(f[0]) + DUPZdeg(other_cofac));
  other_corr = DUPZnew(DUPZdeg(f[0]) + DUPZdeg(cofac));
  tmp = DUPZnew(DUPZdeg(f[0]));
  tmp2 = DUPZnew(DUPZdeg(f[0]));
  for (level=1; level <= degy; level++)
  {
    tmp->deg = -1;
    for (j=1; j < level; j++)
    {
      DUPZmul3(tmp2, factor[j], other_factor[level-j]);
      DUPZadd3(tmp, tmp, tmp2);
    }
    error = DUPZsub(f[level], tmp);
    DUPZmul3(corr, error, other_cofac);
    factor[level] = DUPZrem(corr, factor[0]);
    DUPZdiv2z(factor[level], gcd->coeffs[0]);
    DUPZmul3(other_corr, error, cofac);
    other_factor[level] = DUPZrem(other_corr, other_factor[0]);
    DUPZdiv2z(other_factor[level], gcd->coeffs[0]);
    DUPZfree(error);
  }

  for (i=0; i <= degy; i++) DUPZfree(other_factor[i]);
  FREE(other_factor);
  DUPZfree(cofac); DUPZfree(other_cofac);
  DUPZfree(corr); DUPZfree(other_corr);
  DUPZfree(zero);
  DUPZfree(gcd); DUPZfree(tmp); DUPZfree(tmp2);
}
#endif

DUPZ** DUPZbivariate_lift(const DUPZ* f, DUPZfactors factors)
{
  DUPZ **ans;
  int multiplicity, finished, i, j;
  int degx, degy, nfactors;
  mpz_t j_plus_1;
  DUPZlist iter;

  degx = DUPZdeg(f[0]);
  degy = degx;
  nfactors = DUPZlist_length(factors->list);
  ans = (DUPZ**)MALLOC(nfactors*sizeof(DUPZ*));
  for (i=0, iter = factors->list; i < nfactors; i++, iter = iter->next)
  {
    ans[i] = (DUPZ*)MALLOC(degy*sizeof(DUPZ));
    ans[i][0] = DUPZcopy(iter->poly);
    for (j=1; j <= degy; j++) ans[i][j] = DUPZnew(degx);
  }
  multiplicity = 1;
  finished = 0;
  while (!finished)
  {
    finished = 1;
    for (i=0, iter = factors->list; iter; i++, iter = iter->next)
    {
      if (iter->deg > multiplicity) finished = 0;
      if (iter->deg != multiplicity) continue;
#ifdef FACTOR_DEBUG
      printf("Lifting a factor of multiplicity %d\n", multiplicity);
#endif
      DUPZbivariate_lift_single_factor(ans[i], f, degy);
#ifdef FACTOR_DEBUG
printf("factor lifted:\n");
for(j=degy;j>=0;j--){printf("coeff of y^%d: ",j);DUPZprint(ans[i][j]);}
#endif
    }
    if (finished) break;
#ifdef FACTOR_DEBUG
    printf("differentiating\n");
#endif
    /* Formally differentiate f wrt y */
    for (j=0; j < degy; j++)
    {
      DUPZcopy2(f[j], f[j+1]);
      mpz_set_ui(j_plus_1, j+1);
      DUPZmul2z(f[j], j_plus_1);
    }
    f[degy]->deg = -1; /* set to zero */
    multiplicity++;
  }
  return ans;
}

