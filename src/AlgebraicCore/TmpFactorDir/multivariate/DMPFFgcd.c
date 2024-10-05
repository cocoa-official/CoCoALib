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

#include "DMPFFgcd.h"
#include "primes.h"

/***************************************************************************/
/* GCD of two DMPFFs both of which happen to be univariate.                */

DMPFF DMPFFunivariate_gcd(DMPFF f, DMPFF g, int var)
{
  DUPFF f1, g1;
  DMPFF ans;

  f1 = DMPFF_to_DUPFF(f, var);
  g1 = DMPFF_to_DUPFF(g, var);
  DUPFFgcd2(f1, g1);  /* overwrites f1 with the answer */
  ans = DUPFF_to_DMPFF(f1, var);
  DUPFFdtor(f1);
  DUPFFdtor(g1);
  return ans;
}

/***************************************************************************/
/* Content of DMPFF with respect to a given variable.                      */
/* Result is defined only up to multiplication by a scalar.                */

DMPFF DMPFFcontent(DMPFF f, int var)
{
  DMPFF content, fcopy, top, prev;
  int d;

  content = NULL;
  fcopy = DMPFFcopy(f);
  DMPFForder_main_var(var);
  fcopy = DMPFFsort(fcopy, DMPFForder_main_var_lex);
  while (fcopy != NULL)
  {
    /* Put in top that part of f having degree d in var */
    d = fcopy->exps[var]; /* = DMPFFdeg(fcopy, var) */
    top = fcopy;
    for (; fcopy != NULL && fcopy->exps[var] == d; fcopy = fcopy->next)
    {
      fcopy->exps[var] = 0;
      prev = fcopy;
    }
    prev->next = NULL;
    content = DMPFFgcd(content, top);
    DMPFFdtor(top);
  }
  return content;
}

/***************************************************************************/
/* Multivariate GCD.                                                       */
/* Computes GCD of two DMPFFs.  Result is defined only up to multiplication*/
/* by a scalar.                                                            */

DMPFF DMPFFgcd(DMPFF f, DMPFF g)
{
  int *df, *dg, *monomial_content;

  /* Two trivial cases */
  if (f == NULL) return DMPFFcopy(g);
  if (g == NULL) return DMPFFcopy(f);

  /* Determine any monomial content. */
  df = (int*)MALLOC(NVARS*sizeof(int));
  dg = (int*)MALLOC(NVARS*sizeof(int));
  DUPFFmindegs(df, f);
  DUPFFmindegs(dg, g);

  monomial_content = (int*)MALLOC(NVARS*sizeof(int));
  for (var=0; var < NVARS; var++)
  {
    monomial_content[var] = (df[var] < dg[var]) ? df[var] : dg[var];
    df[var] = -df[var];
    dg[var] = -dg[var];
  }
  /* Divide out the monomial contents */
  DMPFFshift(f, df);
  DMPFFshift(g, dg);

  /* Find out which variables really occur in both polynomials. */
  DUPFFdegs(df, f);
  DUPFFdegs(dg, g);
  for (var=0; var < NVARS; var++)
  {
    if (df[var] == 0 && dg[var] != 0)
    {
      tmp = DMPFFcontent(g, var);
      DMPFFdtor(g);
      g = tmp;
    }
    if (dg[var] == 0 && df[var] != 0)
    {
      tmp = DMPFFcontent(f, var);
      DMPFFdtor(f);
      f = tmp;
    }
  }

  /* Look for any homogeneity ... NYI */

  /* Pick variable order? */
  return DMPFFshift(DMPFFgcd_aux(f, g), monomial_content);
}

/***************************************************************************/
/* Auxiliary function, may destroy its args.  No trivial cases.            */
/* Any variable appearing in f appears in g, and vice versa.               */

DMPFF DMPFFgcd_aux(DMPFF f, DMPFF g)
{
  int var, i, p, nvars;
  int *df, *dg, *dans, *substitution;
  FF Fp;
  DUPFF fx, gx;
  
  DMPFFdegs(df, f);
  nvars = 0;
  for (var = 0; var < NVARS; var++) if (df[var]) nvars++;
  if (nvars == 0) return one;
  if (nvars == 1) return DMPFFunivariate_gcd(f, g, var);

  /* We have at least two variables . */
  DMPFFdegs(dg, g);

  /* Predict degrees of the result in each variable.  We may guess too high. */
  /* For each var map to univariate by random substitution and compute gcd */
  p = CurrentFF.prime;
  for (var = 0; var < NVARS; var++)
  {
    if (df[var] == 0) continue;
    do
    {
      for (i=0; i < NVARS; i++)
        substitution[i] = rand()%p; /* only need values for SOME variables  */
      fx = DMPFF_to_DUPFF(f, var, substitution);
      gx = DMPFF_to_DUPFF(g, var, substitution);
    } while (DUPFFdeg(fx) < df[var] || DUPFFdeg(gx) < dg[var]);
    DUPFFgcd2(fx, gx);
    dans[var] = DUPFFdeg(fx);
    DUPFFdtor(fx);
    DUPFFdtor(gx);
    if (dans[var] > 0) continue;
    /* if var does not appear in the answer, remove it from f and g */
    nvars--;
    if (nvars == 0) return one;
    f = DMPFFcontent(f, var);
    g = DMPFFcontent(g, var);
    if (nvars == 1) return DMPFFunivariate_gcd(f, g, var);
  }
  /* At this point the predicted degrees are probably correct. */

printf("Predicted degrees are: ");for(i=0;i<NVARS;i++)printf("%3d ",dans[vars]);printf("\n");

  /* Pick the 'main' variable. */
  for (var=0; dans[var] == 0; var++);
  fc = DMPFFcontent(f, var);
  gc = DMPFFcontent(g, var);
  f = DMPFFdivexact(f, fc);
  g = DMPFFdivexact(g, gc);
  hc = DMPFFgcd(fc, gc);

  fx = DMPFF_to_DUPFF(f, var, substitution);
  gx = DMPFF_to_DUPFF(g, var, substitution);
  hx = DUPFFgcd2(fx, gx);
  fcofac = DUPFFdiv(fx, hx);
  gcofac = DUPFFdiv(gx, hx);
  tmp = DUPFFgcd2(fcofac, hx);
  if (DUPFFdeg(tmp) == 0)
  {
    DMPFFlift(f, hx, fcofac, substitution, dans);
  }
  tmp = DUPFFgcd2(gcofac, hx);
  if (DUPFFdeg(tmp) == 0)
  {
    DMPFFlift(g, hx, gcofac, substitution, dans);
  }
  for (i=1; i < p; i++)
  {
    DUPFFadd3(fcofac, fcofac, gcofac);
    DUPFFadd3(fx, fx, gx);
    DMPFFadd3(f, f, g);
    tmp = DUPFFgcd(fcofac, hx);
    if (DUPFFdeg(tmp) == 0) break
  }
  DMPFFlift(f, hx, fcofac, substitution, dans);
}


#if 0
/***************************************************************************/
/* Compute full content decomposition f.                                   */
/* For each subset of the variables compute the product of all irreducible */
/* factors of f each containing exactly those variables.  Return a list of */
/* the ones which are non-trivial.                                         */

DMPFFlist DMPFFcontent_decomp(DMPFF f)
{
  
}
#endif
