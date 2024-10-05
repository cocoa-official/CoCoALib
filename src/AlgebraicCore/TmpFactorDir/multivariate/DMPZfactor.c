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

#include "DMPZfactor.h"

#include "DUPZformal_deriv.h"
#include "DUPZgcd.h"
#include "DMPZgcd.h"
#include "DMPZlift.h"
#include "DMPZ_to_DUPZ.h"
#include "DUPZcontent.h"
#include "DUPZfactor.h"
#include "jalloc.h"
#include "DMPZeval.h"



/***************************************************************************/
/* The following two functions deal specially with the case that the       */
/* polynomial is really univariate in the "variable" var.  Note that var   */
/* is in fact a power product.                                             */

/* fin MUST be a polynomial in var (which is itself really a power product) */

static DMPZfactors DMPZfactor_univariate(const DMPZ fin)
{
  DMPZfactors ans;
  DMPZ term, multivariate_copy;
  DUPZ f;
  DUPZfactors facs;
  DUPZlist iter;
  int df, e, var;

  for (var=0; var < NVARS; var++)
    if (DMPZdeg(fin, var) > 0) break;
  /*assert(var < NVARS); */

  df = DMPZdeg(fin, var);
  f = DUPZnew(df);
  /* DMPZs are not canonical, so same exponent may appear several times */
  for (term = fin; term; term = term->next)
  {
    e = term->exps[var];
    mpz_add(f->coeffs[e], f->coeffs[e], term->coeff);
  }
  while (df >= 0 && mpz_sgn(f->coeffs[df]) == 0) df--;
  f->deg = df;
  ans = DMPZfactors_ctor();
  if (df < 1) { DMPZfactors_content(ans, f->coeffs[0]); return ans; }
  facs = DUPZfactor(f);
  DMPZfactors_content(ans, facs->content);
  for (iter = facs->list; iter; iter = iter->next)
  {
    DMPZfactors_multiplicity(ans, iter->deg);
    multivariate_copy = DUPZ_to_DMPZ(iter->poly, var);
    DMPZfactors_add(ans, multivariate_copy);
    DMPZdtor(multivariate_copy);
  }
  DUPZfree(f);
  DUPZfactors_dtor(facs);
  return ans;
}

/***************************************************************************/

static int DMPZfactor_nvars(const DMPZ f)
{
  int nvars, i, *degs;

  if (f == NULL) return 0;
  nvars = 0;
  degs = (int*)MALLOC(NVARS*sizeof(int));
  DMPZdegs(degs, f);
  for (i=0; i < NVARS; i++)
    if (degs[i] != 0) nvars++;
  FREE(degs);
  return nvars;
}


static DMPZ DMPZformal_deriv(const DMPZ f, int var)
{
  DMPZ ans;

  ans = NULL;
  if (f == NULL) return ans;
  ans = DMPZformal_deriv(f->next, var);
  if (f->exps[var] == 0) return ans;
  {
    int *exps, i;
    mpz_t c;
    
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++)
      exps[i] = f->exps[i];
    exps[var]--;
    mpz_init(c);
    mpz_mul_ui(c, f->coeff, f->exps[var]);
    ans = DMPZprepend(c, exps, ans);
    mpz_clear(c);
    return ans;
  }
}


static DMPZlist DMPZsqfr_decomp(const DMPZ fin, int x)
{
  DMPZ f = DMPZcopy(fin);
  DMPZ shadow, newshadow, gcd;
  int i;
  DMPZlist ans = NULL;

#ifdef FACTOR_DEBUG
printf("DMPZsqfr_decomp: called on ");DMPZprint(fin);
#endif
  i = 1;
  {
    DMPZ df = DMPZformal_deriv(f, x);
    gcd = DMPZgcd(f, df);
    DMPZdtor(df);
  }
#ifdef FACTOR_DEBUG
printf("DMPZsqfr_decomp: first gcd is ");DMPZprint(gcd);
#endif

  shadow = DMPZdiv_exact(f, gcd);
#ifdef FACTOR_DEBUG
printf("DMPZsqfr_decomp: shadow is ");DMPZprint(shadow);
#endif
  DMPZdtor(f);
  f = gcd;
  while (DMPZdeg(shadow, x) > 0)
  {
    newshadow = DMPZgcd(f, shadow);
#ifdef FACTOR_DEBUG
printf("DMPZsqfr_decomp: newshadow is ");DMPZprint(newshadow);
#endif
    if (DMPZdeg(newshadow, x) != DMPZdeg(shadow, x))
    {
      DMPZ tmp;
      tmp = DMPZdiv_exact(shadow, newshadow);
      ans = DMPZlist_append(ans, DMPZlist_ctor(tmp, i));
    }
    i++;
    DMPZdtor(shadow);
    shadow = newshadow;
    /* f = DMPZdiv_exact(f, shadow); <-- This would leak so we do the block below instead */
    {
      DMPZ tmp = DMPZdiv_exact(f, shadow);
      DMPZdtor(f);
      f = tmp;
    }
  }
  DMPZdtor(shadow);
  DMPZdtor(f);
  return ans;
}


#if 0
static DMPZlifter DMPZbuild_tree(const DMPZ f, const DUPZlist factors)
{
  DMPZlifter ans;

  ans = (DMPZlifter)MALLOC(sizeof(struct DMPZlifter_struct));
  if (DUPZlist_length(factors) == 2)
  {
    ans = DMPZlifter_ctor(factors->poly, factors->next->poly, NULL, NULL);
  }
}
#endif

static DUPZlist DUPZlist_merge(DUPZlist a, DUPZlist b)
{
  if (a == NULL) return b;
  if (b == NULL) return a;
  if (DUPZdeg(a->poly) > DUPZdeg(b->poly))
  {
    a->next = DUPZlist_merge(a->next, b);
    return a;
  }
  b->next = DUPZlist_merge(a, b->next);
  return b;
}

static DUPZlist DUPZlist_sort(DUPZlist l)
{
  DUPZlist front, back, iter;
  int n;

  if (l == NULL || l->next == NULL) return l;
  front = l;
  iter = l;
  for (n = (DUPZlist_length(l)-1)/2; n > 0; n--) iter = iter->next;
  back = iter->next;
  iter->next = NULL;
  
  return DUPZlist_merge(DUPZlist_sort(front), DUPZlist_sort(back));
}
 


DMPZfactors DMPZfactor(const DMPZ fin)
{
//  const int range = 10;
  int i, x;
  DMPZfactors ans;
  DMPZ f;
  mpz_t content;

  /* Handle the case of univariate input specially */
  if (DMPZfactor_nvars(fin)==1)
  {
    int *degs = (int*)MALLOC(NVARS*sizeof(int));
    DMPZdegs(degs, fin);
    for(i = 0; i < NVARS; i++) if (degs[i] != 0) degs[i] = 1;
    ans = DMPZfactor_univariate(fin);
    FREE(degs);
    return ans;
  }
  ans = DMPZfactors_ctor();
  /* Tedious case: input polynomial is zero. */
  if (fin == NULL)
  {
    mpz_init(content);
    DMPZfactors_content(ans, content);
    mpz_clear(content);
    return ans;
  }
  /* Handle the trivial case of a constant polynomial */
  if (DMPZfactor_nvars(fin) == 0)
  {
    DMPZfactors_content(ans, fin->coeff);
    return ans;
  }
#ifdef FACTOR_DEBUG
  printf("DMPZfactor: not a trivial case.\n");
#endif
  {
    /* Dispose of integer content */
    mpz_t content;
    
    mpz_init(content);
    /* Need to sort f to get signs of factors right.  Order used is */
    /* deglex:  dictated by a comment in DMPZfactors.h and also     */
    /* convenient for DMPZdiv_exact.                                */
    f = DMPZsort(DMPZcopy(fin), DMPZorder_deglex);
    DMPZcontent(content, f);
    DMPZdiv2z(f, content);
    /* Life is easier with a positive leading numerical coefficient */
    if (mpz_sgn(f->coeff) < 0) { DMPZnegate(f); mpz_neg(content, content);}
    DMPZfactors_content(ans, content);
    mpz_clear(content);
  }
#ifdef FACTOR_DEBUG
  printf("DMPZfactor: integer content found.\n");
#endif
  /* Pick a "main" variable, call it x */
  for (x = 0; x < NVARS; x++)
    if (DMPZdeg(f, x) > 0) break;
#ifdef FACTOR_DEBUG
  printf("DMPZfactor: main variable is x[%d]\n", x);
#endif
  {
    /* Dispose of content wrt x */
    DMPZ contentx;
    contentx = DMPZcontent_var(f, x);
    if (mpz_sgn(contentx->coeff) < 0) { DMPZnegate(contentx); }
#ifdef FACTOR_DEBUG
    printf("DMPZfactor: content wrt x[%d] is ", x);DMPZprint(contentx);
#endif
    if (DMPZtotal_deg(contentx) > 0)
    {
      DMPZ tmp;
      DMPZfactors facs;
      DMPZlist iter;
#ifdef FACTOR_DEBUG
      printf("DMPZfactor: factorizing the content wrt x[%d] by recursion...\n", x);
#endif
      facs = DMPZfactor(contentx);
#ifdef FACTOR_DEBUG
      printf("DMPZfactor: Recursive factorization of content wrt x[%d] complete\n", x);
#endif
/* facs->content could be -1, but we sort out the sign later. */
      for (iter = facs->list; iter; iter = iter->next)
      {
        DMPZfactors_multiplicity(ans, iter->deg);
        DMPZfactors_add(ans, iter->poly);
#ifdef FACTOR_DEBUG
        printf("DMPZfactor: factor content wrt x[%d] is: ",x);
        DMPZprint(iter->poly);
#endif
      }
      DMPZfactors_dtor(facs);
      DMPZfactors_multiplicity(ans, 1);
      tmp = DMPZdiv_exact(f, contentx);
      DMPZdtor(f);
      f = tmp;
    }
    DMPZdtor(contentx);
  }

#ifdef FACTOR_DEBUG
  printf("DMPZfactor: now tackling the primitive part wrt x[%d]\n", x);  
  printf("DMPZfactor: primitive part wrt x[%d] is ",x);DMPZprint(f);
#endif
  {
    DMPZlist sqfr_components, iter;

#ifdef FACTOR_DEBUG
  printf("DMPZfactor: About to compute sqfr_decomp of prim part wrt x[%d]\n", x);  
#endif
    sqfr_components = DMPZsqfr_decomp(f, x);
    DMPZdtor(f);
#ifdef FACTOR_DEBUG
  printf("DMPZfactor:  Computed sqfr_decomp of prim part wrt x[%d]\n", x);  
#endif

    for (iter = sqfr_components; iter; iter = iter->next)
    {
      int *flag, *substitution, var;
      DUPZ f1;
      DUPZfactors facs;
      DUPZlist iter1;
    
#ifdef FACTOR_DEBUG
      printf("DMPZfactor: dealing with sqfr compt of multiplicity %d\n", iter->deg);
      printf("DMPZfactor: compt is "); DMPZprint(iter->poly);
      printf("DMPZfactor: next is %lx\n", (long)iter->next);
#endif
      if (DMPZdeg(iter->poly, x) == 0) continue;
      DMPZfactors_multiplicity(ans, iter->deg);
      f = iter->poly; /* alias */
      substitution = (int*)MALLOC(NVARS*sizeof(int));
      try_again:
      /* the following block picks a good substitution: ldcf != 0 & image sqfr */
      {
        DMPZ lcf;
        mpz_t tmp;
      
        lcf = DMPZlc(f, x);
        mpz_init(tmp);
        while(1)
        {
          for (i=0; i < NVARS; i++)
            substitution[i] = rand()%99;
          DMPZeval(tmp, lcf, substitution);
          if (mpz_sgn(tmp) != 0) break;
        }
        mpz_clear(tmp);
        DMPZdtor(lcf);
      }
      flag = (int*)MALLOC(NVARS*sizeof(int));
      for (var = 0; var < NVARS; var++)
        flag[var] = (var != x);
      {
        DMPZ tmp;
        tmp = DMPZeval_partial(f, substitution, flag);
        f1 = DMPZ_to_DUPZ(tmp, x);
        DMPZdtor(tmp);
      }
      FREE(flag);
      {
        DUPZ df = DUPZformal_deriv(f1);
        DUPZ rptdfacs = DUPZgcd(f1, df);
        int sqfr = (DUPZdeg(rptdfacs) == 0);
        DUPZfree(rptdfacs);
        DUPZfree(df);
        if (!sqfr) { DUPZfree(f1); goto try_again; }
      }
      facs = DUPZfactor(f1);
#ifdef FACTOR_DEBUG
      printf("DMPZfactor: chosen substitution is ");for (var = 0; var < NVARS; var++)printf("%d  ",substitution[var]);printf("\n");
      printf("DMPZfactor: univariate image is "); DUPZprint(f1);    
      printf("DMPZfactor: got univariate factorization\n");
      printf("DMPZfactor: degrees are ");
      facs->list = DUPZlist_sort(facs->list);
      for (iter1 = facs->list; iter1; iter1 = iter1->next) printf("%d ", DUPZdeg(iter1->poly));printf("\n");
#endif
      if (DUPZfactors_nfactors(facs) == 1)
      {
        DUPZfree(f1);
        FREE(substitution);
        DUPZfactors_dtor(facs);
        DMPZfactors_add(ans, f);
        continue;
      }
      f = DMPZcopy(f); /* up to here f was an alias, now it is a real copy */
      for (iter1 = facs->list; iter1->next; iter1 = iter1->next)
      {
        DUPZ g1, h1; /* g1 is an alias */
        DMPZlifter lifter;
        g1 = iter1->poly;
        h1 = DUPZdiv(f1, g1);
#ifdef FACTOR_DEBUG
        printf("DMPZfactor: lifting based on ");DUPZprint(g1);DUPZprint(h1);
#endif
        lifter = DMPZlifter_ctor(f, g1, h1, substitution, x);
        DMPZlift(lifter);
        if (lifter->FAILED)
        {
          DMPZfactors newfacs = DMPZfactor(f);
          DMPZlist iter;
          DMPZlifter_dtor(lifter);
          DUPZfree(h1);
#ifdef FACTOR_DEBUG
          printf("DMPZfactor: Uh oh!  Doing recursive call.\n");
#endif
          for (iter=newfacs->list; iter; iter=iter->next)
            DMPZfactors_add(ans, iter->poly);
          DMPZfactors_dtor(newfacs);
          DMPZdtor(f); f = NULL;
          break;
        }
#ifdef FACTOR_DEBUG
        printf("DMPZfactor: got the factors ");DMPZprint(lifter->g);DMPZprint(lifter->h);
#endif
        DMPZfactors_add(ans, lifter->g);
        DMPZdtor(f); f = DMPZcopy(lifter->h);
        DUPZcopy2(f1, h1); /* h1 is the new f1 multiplied by an integer */
        DUPZfree(h1);
        {
          /* h1 may have had some content that we don't want in f1... */
          mpz_t tmp;
          mpz_init(tmp);
          DUPZcontent(tmp, f1);
          DUPZdiv2z(f1, tmp);
          mpz_clear(tmp);
          /* Now f1 is content-free */
        }
        DMPZlifter_dtor(lifter);
      }
      DUPZfree(f1);
      FREE(substitution);
      DUPZfactors_dtor(facs);
      if (f != NULL) DMPZfactors_add(ans, f);
      DMPZdtor(f);
    }
    DMPZlist_dtor(sqfr_components);
    return ans;
  }
}
