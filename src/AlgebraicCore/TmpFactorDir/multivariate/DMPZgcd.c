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

#include <stdlib.h>
#include "jalloc.h"
#include "DMPZgcd.h"
#include "DMPZ_to_DUPFF.h"
#include "DMPZ_to_DUPZ.h"
#include "primes.h"
#include "DUPZgcd.h"
#include "DMPZlift.h"
#include "DMPZmap_to_univariate.h"

int DMPZgcd_level = 0;

/***************************************************************************/
/* Compute content of a (multivariate) polynomial wrt a given variable.    */

DMPZ DMPZcontent_var(const DMPZ f, int var)
{
  int d;
  DMPZ coeff, ans, next;

  ans = NULL;
  for (d=DMPZdeg(f, var); d >= 0; d--)
  {
    coeff = DMPZcoeff(f, var, d);
    next = DMPZgcd(ans, coeff);
    DMPZdtor(coeff);
    DMPZdtor(ans);
    ans = next;
  }

  return ans;
}


/***************************************************************************/
/* GCD of two DMPZs which both happen to be univariate (in the same var).  */
/* At least one of f and g has positive degree.                            */

static DMPZ DMPZunivariate_gcd(const DMPZ f, const DMPZ g)
{
  DUPZ f1, g1, gcd1;
  DMPZ ans;
  int var;
  for (var=0; var < NVARS; var++)
    if (DMPZdeg(f, var) > 0 || DMPZdeg(g, var) > 0) break;
  /*assert(var < NVARS); */

  f1 = DMPZ_to_DUPZ(f, var);
  g1 = DMPZ_to_DUPZ(g, var);
  gcd1 = DUPZgcd(f1, g1);
  DUPZfree(f1); DUPZfree(g1);
  ans = DUPZ_to_DMPZ(gcd1, var);
  DUPZfree(gcd1);
  return ans;
}


/***************************************************************************/
/* Prepare for a multivariate lift.  The catch is that the univariate      */
/* images may not be coprime -- a prerequisite for lifting.  We solve this */
/* by adding some integer multiple of g to f.  Sooner or later we are      */
/* guaranteed to have two coprime univariate factors which we can lift.    */
/* If lifting fails for some reason, we return NULL (equivalent to 0).     */
/* ASSUMES that f, g are primitive wrt var, and have no integer content.   */

static DMPZ DMPZgcd_lift(const DMPZ f, const DMPZ g, const DUPZ f1_, const DUPZ g1, int var, int *a)
{
  DMPZ F;     /* will be a copy of f once trivial cases have been eliminated */
  DMPZ G = g; /* Just an alias */
  DUPZ gcd1;
  DUPZ fcofac, gcofac;
  DUPZ f1;
  DMPZlifter lifter;

  /* We first check to see if the GCD is equal to f or g.                 */
  /* This is partly for efficiency, but mostly because the lifter crashes */
  /* if this case is allowed to reach the general code below.             */
  gcd1 = DUPZgcd(f1_, g1);
  if (DUPZdeg(gcd1) == 0) /* Trivial case: polys are coprime */
  {
    DUPZfree(gcd1);
    return DMPZone();
  }
  if (DUPZdeg(gcd1) == DUPZdeg(f1_)) /* GCD = f if we have good reduction */
  {
    DMPZ quot;
    DUPZfree(gcd1);
    quot = DMPZdiv_exact(g, f);
    if (quot == NULL) return NULL;
    DMPZdtor(quot);
    return DMPZcopy(f);
  }
  if (DUPZdeg(gcd1) == DUPZdeg(g1)) /* GCD = g if we have good reduction */
  {
    DMPZ quot;
    DUPZfree(gcd1);
    quot = DMPZdiv_exact(f, g);
    if (quot == NULL) return NULL;
    DMPZdtor(quot);
    return DMPZcopy(g);
  }
  /* These 5 lines allow f1 to contain a poly of degree up to max(deg(f1_), deg(g1)). */
  /* Then initialize f1 to value of f1_.  Similarly allow fcofac to have degree up to */
  /* max(deg(f1), deg(g1)) - deg(gcd1).  The extra space is needed if the loop below  */
  /* has to add an integer multiple of g1 to f1 to make fcofac and gcd1 coprime.      */
  if (DUPZdeg(g1) > DUPZdeg(f1_)) f1 = DUPZnew(DUPZdeg(g1));
  else f1 = DUPZnew(DUPZdeg(f1_));
  DUPZcopy2(f1, f1_);
  if (DUPZdeg(g1) > DUPZdeg(f1_)) fcofac = DUPZnew(DUPZdeg(g1)-DUPZdeg(gcd1));
  else fcofac = DUPZnew(DUPZdeg(f1)-DUPZdeg(gcd1));

  /*  The block below effectively computes fcofac = DUPZdiv(f1, gcd1); */
  {
    DUPZ zero = DUPZnew(-1);
    
    DUPZdiv4(fcofac, zero, f1, gcd1); /* This will error if the division is not exact. */
    DUPZfree(zero);
  }
  gcofac = DUPZdiv(g1, gcd1);

  /* Must have fcofac and gcd1 coprime for lifting to work. */
  /* Achieve this by adding an integer multiple of G to F.  */
  F = DMPZcopy(f);
  while (1)
  {
    int d;
    {
      DUPZ tmp = DUPZgcd(fcofac, gcd1);
      d = DUPZdeg(tmp);
      DUPZfree(tmp);
    }
    if (d == 0) break;
    d = DUPZdeg(f1);
    once_more:
    DUPZadd3(f1, f1, g1);
    DUPZadd3(fcofac, fcofac, gcofac);
    DMPZadd3(&F, F, G);
    if (DUPZdeg(f1) < d) goto once_more; /* cond is true only if deg(f1)=deg(g1) && lc(f1)+lc(g1)=0. */
                                         /* Can happen at most once. */
  }

  lifter = DMPZlifter_ctor(F, gcd1, fcofac, a, var);
  DMPZlift(lifter);
  DMPZdtor(F);
  DUPZfree(f1);
  DUPZfree(fcofac);
  DUPZfree(gcofac);
  DUPZfree(gcd1);
  if (lifter->FAILED) { DMPZlifter_dtor(lifter); return NULL; }
  {
    DMPZ ans, gcofac;
    ans = lifter->g;
    lifter->g = NULL;
    DMPZlifter_dtor(lifter);
    /* Must check that the result also divides g -- it could be too big */
    gcofac = DMPZdiv_exact(g, ans);
    if (gcofac == NULL) { DMPZdtor(ans); return NULL; }
    DMPZdtor(gcofac);
    return ans;
  }
}


/***************************************************************************/
/* Auxiliary function, may destroy its args.                               */
/* Any variable appearing in f appears in g, and vice versa.               */
/* No integer or monomial content.                                         */

static DMPZ DMPZgcd_aux(const DMPZ fin, const DMPZ gin)
{
  int var, i, niters;
  int *df, *dg, *dans, *substitution;
  DMPZ f, g, content, ans;
  
#ifdef FACTOR_DEBUG
printf("DMPZgcd_aux: fin = ");DMPZprint(fin);  
printf("DMPZgcd_aux: gin = ");DMPZprint(gin);  
#endif
  df = (int*)MALLOC(NVARS*sizeof(int));
  DMPZdegs(df, fin);

  /* The next block looks for two simple special cases (i.e. nvars <= 1). */
  {
    int nvars = 0;
    for (var = 0; var < NVARS; var++)
      if (df[var]) nvars++;
    if (nvars == 0) { FREE(df); return DMPZone(); }
    if (nvars == 1) { FREE(df); return DMPZunivariate_gcd(fin, gin); }
  }

#ifdef FACTOR_DEBUG
printf("DMPZgcd_aux: not a trivial case\n");
#endif
  /* We have at least two variables. */
  dg = (int*)MALLOC(NVARS*sizeof(int));
  DMPZdegs(dg, gin);
  f = DMPZcopy(fin);
  g = DMPZcopy(gin);
#ifdef FACTOR_DEBUG
printf("Really inside DMPZgcd; degrees are as follows:\n");
printf("f:");for(i=0;i<NVARS;i++)printf("%3d ",df[i]);printf("\n");
printf("g:");for(i=0;i<NVARS;i++)printf("%3d ",dg[i]);printf("\n");
#endif

  /* Predict degrees of the result in each variable.  We may guess too high. */
  dans = (int*)MALLOC(NVARS*sizeof(int));

  substitution = (int*)MALLOC(NVARS*sizeof(int));
  for (var = 0; var < NVARS; var++)
  {
    int p;
    if (df[var] == 0) { dans[var] = 0; continue; }
/*    printf("loop: doing var x(%d)\n",var);                                */
/*    printf("loop: orig degs are deg(f)=%d  deg(g)=%d\n",df[var],dg[var]); */
    for (p = 7; ; p = nextprime(p))
    {
      DUPFF fx, gx;
      FF Fp = FFctor(p);
      FFselect(Fp);
/*      printf("inner loop: p=%d\n",p); */
      for (i=0; i < NVARS; i++)
        substitution[i] = rand()%p; /* need values for only SOME of the variables     */
/*      printf("inner loop: subst=[");for(i=0;i<NVARS;++i)printf("%d ",substitution[i]);printf("]\n");*/
      fx = DMPZ_to_DUPFF(f, var, substitution);
      gx = DMPZ_to_DUPFF(g, var, substitution);
/*      printf("inner loop: deg(fx)=%d  deg(gx)=%d\n",DUPFFdeg(fx),DUPFFdeg(gx));  */
/*      printf("inner loop: fx="); DUPFFprint(fx);                                 */
/*      printf("inner loop: gx="); DUPFFprint(gx);                                 */
      if (DUPFFdeg(fx) < df[var] || DUPFFdeg(gx) < dg[var]) { DUPFFfree(gx); DUPFFfree(fx); FFdtor(Fp); continue; }
      DUPFFgcd2(fx, gx);
      dans[var] = DUPFFdeg(fx);
      DUPFFfree(gx);
      DUPFFfree(fx);
      FFdtor(Fp);
      break;
    }
    if (dans[var] > 0) continue;
    /* if var does not appear in the answer, remove it from f and g */
#ifdef FACTOR_DEBUG
    printf("DMPZgcd_aux: var[%d] does not appear in answer, so remove it by taking contents\n",var);
#endif
    {
      /* This block simply replaces f and g by their contents w.r.t. var */
      /* Cannot use assignment as it is not defined for DMPZs -- naive assignment leaks! */
      DMPZ tmp = DMPZcontent_var(f, var);
      DMPZdtor(f); f = tmp;
      tmp = DMPZcontent_var(g, var);
      DMPZdtor(g); g = tmp;
    }
#ifdef FACTOR_DEBUG
    printf("DMPZgcd_aux: after removal of var[%d] f and g are:\n",var);
    printf("DMPZgcd_aux: f=");DMPZprint(f);
    printf("DMPZgcd_aux: g=");DMPZprint(g);
#endif
    DMPZdegs(df, f);
    DMPZdegs(dg, g);
    var = -1;
  }
  FREE(substitution);
  /* At this point the predicted degrees are probably correct. */

#ifdef FACTOR_DEBUG
printf("Predicted degrees are: ");for(i=0;i<NVARS;i++)printf("%3d ",dans[i]);printf("\n");
#endif

  /* We pick a variable to be retained in univariate image. */
  /* FOR THE MOMENT we just take the first one we find. */
  for (var = 0; var < NVARS; var++)
    if (dans[var] != 0) break;
  if (var == NVARS)  /* This means we have eliminated all the variables! */
  {
    FREE(dans);
    FREE(dg);
    FREE(df);
    DMPZdtor(g);
    DMPZdtor(f);
    return DMPZone();
  }
#ifdef FACTOR_DEBUG
printf("variable to be kept in univariate image: x(%d)\n", var);
#endif
  {
    DMPZ contentf, contentg;
    contentf = DMPZcontent_var(f, var);
    contentg = DMPZcontent_var(g, var);
    {
      /* This block simply replaces f and g by their primitive parts w.r.t. var */
      /* Cannot use assignment as it is not defined for DMPZs -- naive assignment leaks! */
      DMPZ tmp = DMPZdiv_exact(f, contentf);
      DMPZdtor(f); f = tmp;
      tmp = DMPZdiv_exact(g, contentg);
      DMPZdtor(g); g = tmp;
    }
    content = DMPZgcd(contentf, contentg);
    DMPZdtor(contentf);
    DMPZdtor(contentg);
  }
#ifdef FACTOR_DEBUG
printf("Removed content wrt x[%d]\n", var);  
printf("f=");DMPZprint(f);
printf("g=");DMPZprint(g);
#endif
/*
 *  lcf = DMPZcoeff(f, var, df[var]);
 *  lcg = DMPZcoeff(g, var, dg[var]);
 *  lcans = DMPZgcd(lcf, lcg);
 */
  
  for (niters=1; ; ++niters)
  {
    DUPZ f1, g1;
    int *substitution = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++)
      substitution[i] = (rand()%(9+10*niters)); /* only need values for SOME of the variables */
                                                /* Modulus is a heuristic/guess. */

#ifdef FACTOR_DEBUG
printf("Map is : ");
for (i=0; i < NVARS; i++)
  printf("%5d",substitution[i]);
printf("\n");
#endif
    f1 = DMPZmap_to_univariate(f, var, substitution);
    g1 = DMPZmap_to_univariate(g, var, substitution);
#ifdef FACTOR_DEBUG
printf("DMPZgcd_aux: univariate image of f is "); DUPZprint(f1);
printf("DMPZgcd_aux: univariate image of g is "); DUPZprint(g1);
#endif
    if (DUPZdeg(f1) < df[var] || DUPZdeg(g1) < dg[var])
    {
      DUPZfree(f1); DUPZfree(g1);
      FREE(substitution);
      continue;
    }
    ans = DMPZgcd_lift(f, g, f1, g1, var, substitution);
    DUPZfree(f1); DUPZfree(g1);
    FREE(substitution);
    if (ans != NULL) break;
  }

  DMPZdtor(f);
  DMPZdtor(g);
  FREE(df);
  FREE(dg);
  FREE(dans);
  {
    DMPZ tmp;
    tmp = DMPZmul(content, ans);
    DMPZdtor(content);
    DMPZdtor(ans);
    ans = tmp;
  }
  return ans;
}

/***************************************************************************/
/* Multivariate GCD -- general entry point.                                */
/* fin and gin satisfy no special criteria.  Essentially this function     */
/* removes content repeatedly to reduce to the case of computing the GCD   */
/* of two polynomials which involve the same variables, and which are      */
/* primitive wrt a given one of those variables.                           */

DMPZ DMPZgcd(const DMPZ fin, const DMPZ gin)
{
  int var;
  int *df, *dg, *monomial_content;
  DMPZ f, g, ans;
  mpz_t integer_content;

  /* Two trivial cases */
  if (fin == NULL) { g = DMPZcopy(gin); if (mpz_sgn(g->coeff) < 0) { DMPZnegate(g); } return g; }

  if (gin == NULL) { f = DMPZcopy(fin); if (mpz_sgn(f->coeff) < 0) { DMPZnegate(f); } return f; }
#ifdef FACTOR_DEBUG
  ++DMPZgcd_level;
printf("DMPZgcd(%d): ENTERED WITH\n", DMPZgcd_level);
printf("DMPZgcd(%d): fin=",DMPZgcd_level);DMPZprint(fin);
printf("DMPZgcd(%d): gin=",DMPZgcd_level);DMPZprint(gin);
#endif
  mpz_init(integer_content);
  {
    mpz_t cf, cg;
    /* Remove integer contents. */
    f = DMPZcopy(fin);
    g = DMPZcopy(gin);
    mpz_init(cf);
    mpz_init(cg);
    DMPZcontent(cf, f);
    DMPZcontent(cg, g);
    DMPZdiv2z(f, cf);
    DMPZdiv2z(g, cg);
    mpz_gcd(integer_content, cf, cg);
    mpz_clear(cf);
    mpz_clear(cg);
  }
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): after removing integer contents we have...\n",DMPZgcd_level);
printf("DMPZgcd(%d): f=",DMPZgcd_level);DMPZprint(f);
printf("DMPZgcd(%d): g=",DMPZgcd_level);DMPZprint(g);
#endif
  /* Determine any monomial content. */
  df = (int*)MALLOC(NVARS*sizeof(int));
  dg = (int*)MALLOC(NVARS*sizeof(int));
  monomial_content = (int*)MALLOC(NVARS*sizeof(int));
  DMPZmindegs(df, f);
  DMPZmindegs(dg, g);
  for (var=0; var < NVARS; var++)
  {
    monomial_content[var] = (df[var] < dg[var]) ? df[var] : dg[var];
    df[var] = -df[var];
    dg[var] = -dg[var];
  }
  /* Divide out the monomial contents */
  DMPZshift(f, df);
  DMPZshift(g, dg);
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): after removing monomial contents we have...\n",DMPZgcd_level);
printf("DMPZgcd(%d): f=",DMPZgcd_level);DMPZprint(f);
printf("DMPZgcd(%d): g=",DMPZgcd_level);DMPZprint(g);
#endif

  /* Find out which variables really occur in both polynomials. */
  DMPZdegs(df, f);
  DMPZdegs(dg, g);
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): degrees in the variables are\n",DMPZgcd_level);
printf("DMPZgcd(%d): f: ",DMPZgcd_level);for(var=0;var<NVARS;var++)printf("%3d ",df[var]);printf("\n");
printf("DMPZgcd(%d): g: ",DMPZgcd_level);for(var=0;var<NVARS;var++)printf("%3d ",dg[var]);printf("\n");
#endif
restart_loop:
  for (var=0; var < NVARS; var++)
  {
    if (df[var] == 0 && dg[var] != 0)
    {
      DMPZ tmp;
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): eliminating variable x[%d] from g by content computation.\n",DMPZgcd_level,var);
#endif
      tmp = DMPZcontent_var(g, var);
      DMPZdtor(g);
      g = tmp;
      DMPZdegs(dg, g);
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): eliminated variable x[%d] from g by content computation.\n",DMPZgcd_level,var);
printf("DMPZgcd(%d): g: ",DMPZgcd_level);for(var=0;var<NVARS;var++)printf("%3d ",dg[var]);printf("\n");
#endif
      goto restart_loop;
    }
    if (dg[var] == 0 && df[var] != 0)
    {
      DMPZ tmp;
#ifdef FACTOR_DEBUG
printf("DMPZgcd: eliminating variable x[%d] from f by content computation.\n",var);
#endif
      tmp = DMPZcontent_var(f, var);
      DMPZdtor(f);
      f = tmp;
      DMPZdegs(df, f);
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): eliminated variable x[%d] from f by content computation.\n",DMPZgcd_level,var);
printf("DMPZgcd(%d): f: ",DMPZgcd_level);for(var=0;var<NVARS;var++)printf("%3d ",df[var]);printf("\n");
#endif
      goto restart_loop;
    }
  }
  FREE(df);
  FREE(dg);
  /* Look for any homogeneity ... NYI */
#ifdef FACTOR_DEBUG
printf("DMPZgcd(%d): after removing unnecessary variables we have...\n",DMPZgcd_level);
printf("DMPZgcd(%d): f=",DMPZgcd_level);DMPZprint(f);
printf("DMPZgcd(%d): g=",DMPZgcd_level);DMPZprint(g);
#endif

  /* Pick variable order? */
  ans = DMPZgcd_aux(f, g);
  if (mpz_sgn(ans->coeff) < 0) { DMPZnegate(ans); } /* force LC > 0 */
  DMPZdtor(f);
  DMPZdtor(g);
  {
    DMPZ content, tmp;
    content = DMPZctor(integer_content, monomial_content);
    tmp = DMPZmul(ans, content);
    DMPZdtor(ans);
    DMPZdtor(content); /* WARNING! This effectively includes   FREE(monomial_content) */
    ans = tmp;
  }
  mpz_clear(integer_content);
#ifdef FACTOR_DEBUG
  printf("DMPZgcd(%d): RETURNING: ",DMPZgcd_level);DMPZprint(ans);
  --DMPZgcd_level;
#endif
  return ans;
}

