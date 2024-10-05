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

/* only needed for DMPZprint -- for debugging really */
#include <stdio.h>


#include "DMPZ.h"


#include <stdlib.h>
#include "jalloc.h"
#include "jaaerror.h"

int NVARS = -1; /* silly initial value, should provoke errors if uninitialized */

/***************************************************************************/

#if 0
/* This should be normal lexicographic ordering */
static int PPcmp(int *pp1, int *pp2)
{
  int var;

  for(var=0; var < NVARS-1; var++)
    if (pp1[var] != pp2[var]) break;
  return pp1[var] - pp2[var];
}
#endif


DMPZ DMPZctor(const mpz_t coeff, int *exps)
{
  DMPZ ans;

  ans = (DMPZ)MALLOC(sizeof(struct DMPZ_struct));
  mpz_init_set(ans->coeff, coeff);
  ans->exps = exps;
  ans->next = NULL;
  return ans;
}

void DMPZdtor(DMPZ f)
{
  DMPZ THIS, next;

  for (THIS = f; THIS; THIS = next)
  {
    next = THIS->next;
    mpz_clear(THIS->coeff);
    FREE(THIS->exps);
    FREE(THIS);
  }
}


/* A useful variant of DMPZctor. */
DMPZ DMPZprepend(const mpz_t coeff, int *exps, DMPZ f)
{
  DMPZ ans;

  ans = DMPZctor(coeff, exps);
  ans->next = f;
  return ans;
}


/* Destructive reversal */
DMPZ DMPZreverse(DMPZ f)
{
  DMPZ THIS, prev, next;

  if (f == NULL || f->next == NULL) return f;
  prev = NULL;
  THIS = f;
  while (THIS)
  {
    next = THIS->next;
    THIS->next = prev;
    prev = THIS;
    THIS = next;
  }
  return prev;
}

DMPZ DMPZcopy(const DMPZ f)
{
  DMPZ ans, iter, tmp;
  int i, *exps_copy;

  ans = NULL;
  for (iter=f; iter; iter = iter->next)
  {
    exps_copy = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++) exps_copy[i] = iter->exps[i];
    tmp = DMPZctor(iter->coeff, exps_copy);
    tmp->next = ans;
    ans = tmp;
  }
  return DMPZreverse(ans);
}


DMPZ DMPZone(void)
{
  int i;
  DMPZ ans = (DMPZ)MALLOC(sizeof(struct DMPZ_struct));

  ans->next = NULL;
  ans->exps = (int*)MALLOC(NVARS*sizeof(int));
  mpz_init_set_ui(ans->coeff, 1);
  for (i=0; i < NVARS; i++) ans->exps[i] = 0;
  return ans;
}

/***************************************************************************/

int DMPZ_nterms(const DMPZ f)
{
  int i;
  DMPZ iter;

  i = 0;
  for (iter=f; iter; iter=iter->next) i++;
  return i;
}


/***************************************************************************/
/* Makes a copy of the leading coeff of f wrt the variable var.            */

DMPZ DMPZlc(const DMPZ f, int var)
{
  return DMPZcoeff(f, var, DMPZdeg(f, var));
}


/***************************************************************************/
/* Makes a copy of that part of f for which the degree of var equals deg   */
/* and divides it by var^deg.  The coefficient of var^deg in f.            */
/*
DMPZ DMPZcoeff(const DMPZ f, int var, int deg)
{
  DMPZ leading_term;
  int *exps_copy;
  int i;
  
  if (f == NULL) return NULL;
  if (f->exps[var] != deg) return DMPZcoeff(f->next, var, deg);
  exps_copy = (int*)MALLOC(NVARS*sizeof(int));
  for (i=0; i < NVARS; i++) exps_copy[i] = f->exps[i];
  exps_copy[var] = 0;
  leading_term = DMPZctor(f->coeff, exps_copy);
  leading_term->next = DMPZcoeff(f->next, var, deg);
  return leading_term;
}
*/
DMPZ DMPZcoeff(const DMPZ f, int var, int deg)
{
  DMPZ iter, ans;
  int *exps_copy;
  int i;

  ans = NULL;
  for (iter = f; iter; iter = iter->next)
  {
    if (iter->exps[var] != deg) continue;
    exps_copy = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++) exps_copy[i] = iter->exps[i];
    exps_copy[var] = 0;
    ans = DMPZprepend(iter->coeff, exps_copy, ans);    
  }
  return DMPZreverse(ans);
}


/***************************************************************************/

DMPZ* DMPZcoeffs(const DMPZ f, int x)
{
  int dx, i, j;
  DMPZ *ans, iter;

  dx = DMPZdeg(f, x);
  ans = (DMPZ*)MALLOC((1+dx)*sizeof(DMPZ));
  for (i=0; i <= dx; i++) ans[i] = NULL;
  for (iter=f; iter; iter = iter->next)
  {
    int *exps;
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (j=0; j < NVARS; j++) exps[j] = iter->exps[j];
    exps[x] = 0;
    ans[iter->exps[x]] = DMPZprepend(iter->coeff, exps, ans[iter->exps[x]]);
  }
  return ans;
}

/***************************************************************************/
/* Degree and multi-degree of a DMPZ.                                      */

void DMPZdegs(int *degs, const DMPZ f)
{
  DMPZ iter;
  int var;

  for (var=0; var < NVARS; var++) degs[var] = 0;
  for (iter=f; iter; iter = iter->next)
    for (var=0; var < NVARS; var++)
      if (degs[var] < iter->exps[var]) degs[var] = iter->exps[var];
}

int DMPZdeg(const DMPZ f, int var)
{
  DMPZ iter;
  int deg = 0;

  for (iter=f; iter; iter = iter->next)
    if (deg < iter->exps[var]) deg = iter->exps[var];
  return deg;
}

int DMPZtotal_deg(const DMPZ f)
{
  int max, curr, i;
  DMPZ iter;
  
  max = -1;
  for (iter=f; iter; iter = iter->next)
  {
    curr = 0;
    for (i=0; i < NVARS; i++)
      curr += iter->exps[i];
    if (curr > max) max = curr;
  }
  return max;
}

void DMPZmindegs(int *degs, const DMPZ f)
{
  DMPZ iter;
  int var;

  for (var=0; var < NVARS; var++) degs[var] = 0;
  if (f == NULL) return;
  for (var=0; var < NVARS; var++) degs[var] = f->exps[var];
  for (iter=f->next; iter; iter = iter->next)
    for (var=0; var < NVARS; var++)
      if (degs[var] > iter->exps[var]) degs[var] = iter->exps[var];
}


/***************************************************************************/
/* Compute integer content of a DMPZ polynomial and put result in content. */
/* The algorithm is not terribly clever.                                   */
/* If f = 0 result is zero, otherwise it is positive.                      */

void DMPZcontent(mpz_t content, const DMPZ f)
{
  DMPZ iter, min_ptr;
  unsigned int min;

  if (f == NULL) { mpz_set_ui(content, 0); return; } /* silly case, f==0 */
  /* Start content with value of the "smallest" coefficient. */
  min_ptr = f;
  min = mpz_sizeinbase(f->coeff, 2);
  for (iter=f; iter; iter = iter->next)
    if (mpz_sizeinbase(iter->coeff, 2) < min)
    {
      min = mpz_sizeinbase(iter->coeff, 2);
      min_ptr = iter;
    }
  mpz_abs(content, min_ptr->coeff);

  /* Now gcd content with all the other coefficients. */
  for (iter=f; (mpz_cmp_ui(content, 1) != 0) && iter; iter = iter->next)
    if (iter != min_ptr)
      mpz_gcd(content, content, iter->coeff);
}


/***************************************************************************/
/* Negate a DMPZ in place, i.e. "destructive".                             */
void DMPZnegate(DMPZ f)
{
  DMPZ iter;
  for (iter=f; iter; iter = iter->next)
    mpz_neg(iter->coeff, iter->coeff);
}


/***************************************************************************/
/* This function is used by DMPZdiv_exact below.                           */

static DMPZ DMPZshift_sub(DMPZ f, DMPZ g, DMPZ monom)
{
  int var, *exps;
  mpz_t coeff;
  DMPZ f_iter, g_iter, ans;
  int cmp, update_exps;
  
  ans = NULL;
  exps = NULL; /* to keep the compiler quiet */
  g_iter = g;
  f_iter = f;
  update_exps = 1;
  while (f_iter || g_iter)
  {
    if (g_iter == NULL) goto step_f;
    if (update_exps)
    {
      exps = (int*)MALLOC(NVARS*sizeof(int));
      for (var=0; var < NVARS; var++)
        exps[var] = g_iter->exps[var] + monom->exps[var];
      update_exps = 0;
    }
    if (f_iter == NULL) goto step_g;
    cmp = DMPZorder_deglex(f_iter->exps, exps);
    if (cmp > 0) goto step_f;
    if (cmp < 0) goto step_g;
    mpz_init(coeff);
    mpz_mul(coeff, g_iter->coeff, monom->coeff);
    mpz_sub(coeff, f_iter->coeff, coeff);
    if (mpz_sgn(coeff)) ans = DMPZprepend(coeff, exps, ans);
    else FREE(exps);
    mpz_clear(coeff);
    f_iter = f_iter->next;
    g_iter = g_iter->next;
    update_exps = 1;
    continue;

    step_f:
    {
      int *expf = (int*)MALLOC(NVARS*sizeof(int));
      for (var=0; var < NVARS; var++) expf[var] = f_iter->exps[var];
      ans = DMPZprepend(f_iter->coeff,expf, ans);
      f_iter = f_iter->next;
      continue;
    }
    step_g:
    mpz_init(coeff);
    mpz_mul(coeff, g_iter->coeff, monom->coeff);
    mpz_neg(coeff, coeff);
    ans = DMPZprepend(coeff, exps, ans);
    mpz_clear(coeff);
    g_iter = g_iter->next;
    update_exps = 1;
    continue;
  }

  return DMPZreverse(ans);
}

/***************************************************************************/
/* Exact division of one multivariate polynomial by another.               */
/* Result is NULL if the division is not exact (or if the first arg is 0). */

DMPZ DMPZdiv_exact(const DMPZ fin, const DMPZ gin)
{
  mpz_t q, r;
  int *exps, var;
  DMPZ f, g, quot, tmp;

  if (gin == NULL) JERROR(JERROR_DIV_BY_ZERO);
  if (fin == NULL) return NULL;
  /* If you change the order below, you must also change DMPZfactor(...) */
  f = DMPZsort(DMPZcopy(fin), DMPZorder_deglex);
  g = DMPZsort(DMPZcopy(gin), DMPZorder_deglex);
  quot = NULL;
  mpz_init(q);
  mpz_init(r);
  while(f)
  {
    mpz_fdiv_qr(q, r, f->coeff, g->coeff);
    if (mpz_sgn(r)) { DMPZdtor(quot); quot = NULL; break; }
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (var=0; var < NVARS; var++) exps[var] = f->exps[var] - g->exps[var];
    for (var=0; var < NVARS; var++)
      if (exps[var] < 0) break; /* no need to free exps, as it gets owned by quot then freed */
    quot = DMPZprepend(q, exps, quot);
    if (var < NVARS) { DMPZdtor(quot); quot = NULL; break; }
    tmp = DMPZshift_sub(f, g, quot);
    DMPZdtor(f);
    f = DMPZsort(tmp, DMPZorder_deglex);
  }
  mpz_clear(r);
  mpz_clear(q);
  DMPZdtor(f);
  DMPZdtor(g);
  return DMPZreverse(quot);
}


void DMPZdiv2z(DMPZ f, const mpz_t n)
{
  DMPZ iter;

  for (iter=f; iter; iter = iter->next)
    mpz_divexact(iter->coeff, iter->coeff, n);
}


void DMPZmul2z(DMPZ f, const mpz_t n)
{
  DMPZ iter;

  for (iter=f; iter; iter = iter->next)
    mpz_mul(iter->coeff, iter->coeff, n);
}


DMPZ DMPZmul(const DMPZ f, const DMPZ g)
{
  DMPZ ans = NULL;
  DMPZ iterf, iterg;
  mpz_t coeff;
  int i;

  if (f == NULL || g == NULL) return NULL;
  mpz_init(coeff);
  for (iterf = f; iterf; iterf = iterf->next)
    for (iterg = g; iterg; iterg = iterg->next)
    {
      int *exps = (int*)MALLOC(NVARS*sizeof(int));
      for (i=0; i < NVARS; i++)
        exps[i] = iterf->exps[i] + iterg->exps[i];
      mpz_mul(coeff, iterf->coeff, iterg->coeff);
      ans = DMPZprepend(coeff, exps, ans);
    }
  mpz_clear(coeff);
  return DMPZsort(ans, DMPZorder_lex);
}


DMPZ DMPZshift(DMPZ f, const int *shift)
{
  DMPZ iter;
  int i;

  for (iter=f; iter; iter = iter->next)
  {
    for (i=0; i < NVARS; i++)
      iter->exps[i] += shift[i];
  }
  return f;
}


/***************************************************************************/
/* A power product ordering: lexicographic.                                */
/* Higher powers of a variable are "bigger" than lower powers, if powers   */
/* equal look at next variable...                                          */

int DMPZorder_lex(int* a, int* b)
{
  int i;

  for (i=0; i < NVARS; i++)
    if (a[i] != b[i]) return (a[i] - b[i]);
  return 0;
}

int DMPZorder_deglex(int* a, int* b)
{
  int i, da, db;

  da = db = 0;
  for (i=0; i < NVARS; i++) { da += a[i]; db += b[i]; }
  if (da != db) return da-db;
  return DMPZorder_lex(a, b);
}


/***************************************************************************/
/* An implementation of "merge" sort on a linked list (= terms in the poly)*/
/* Terms are sorted into DECREASING power product order according to cmp.  */

/* // this is old WORKING VERSION */
/* static DMPZ DMPZmerge(DMPZ a, DMPZ b, DMPZexps_cmp cmp) */
/* { */
/*   int comparison; */
/*   if (a == NULL) return b; */
/*   if (b == NULL) return a; */

/*   comparison = cmp(a->exps, b->exps); */
/*   if (comparison > 0) */
/*   { */
/*     a->next = DMPZmerge(a->next, b, cmp); */
/*     return a; */
/*   } */
/*   if (comparison < 0) */
/*   { */
/*     b->next = DMPZmerge(a, b->next, cmp); */
/*     return b; */
/*   } */
/*   /\* the two exponents are equal so merge the terms *\/ */
/*   mpz_add(a->coeff, a->coeff, b->coeff); */
/*   { */
/*     DMPZ a_next, b_next; */
/*     b_next = b->next; */
/*     b->next = NULL; */
/*     DMPZdtor(b); */
/*     if (mpz_sgn(a->coeff)) { a->next = DMPZmerge(a->next, b_next, cmp); return a; } */
/*     a_next = a->next; */
/*     a->next = NULL; */
/*     DMPZdtor(a); */
/*     return DMPZmerge(a_next, b_next, cmp); */
/*   } */
/* } */


/* // hacked version 2015-09-17 -- seems to work!?! */
/* static DMPZ DMPZmerge(DMPZ a, DMPZ b, DMPZexps_cmp cmp) */
/* { */
/*   int comparison; */

/*   DMPZ ans = NULL; */
/*   DMPZ *tail = &ans; */
/*   loop: */
/*   if (a == NULL) { *tail = b; return ans; } */
/*   if (b == NULL) { *tail = a; return ans; } */

/*   comparison = cmp(a->exps, b->exps); */
/*   if (comparison > 0) */
/*   { */
/*     *tail = a; */
/*     tail = &((*tail)->next); */
/*     a = a->next; */
/*     goto loop; */
/*   } */
/*   if (comparison < 0) */
/*   { */
/*     *tail = b; */
/*     tail = &((*tail)->next); */
/*     b = b->next; */
/*     goto loop; */
/*   } */
/*   /\* the two exponents are equal so merge the terms *\/ */
/*   mpz_add(a->coeff, a->coeff, b->coeff); */
/*   { */
/*     if (mpz_sgn(a->coeff) == 0) { } */
/*     DMPZ junk; */
/*     junk = b; */
/*     b = b->next; */
/*     junk->next = NULL; */
/*     DMPZdtor(junk); */
/*     if (mpz_sgn(a->coeff)) {     *tail = a; */
/*       tail = &((*tail)->next); */
/*     a = a->next; */
/*     goto loop; */
/*     } */
/*     junk = a; */
/*     a = a->next; */
/*     junk->next = NULL; */
/*     DMPZdtor(junk); */
/*     goto loop; */
/*   } */
/* } */


/* Destructively add a and b; effectively merge the lists */
static DMPZ DMPZmerge(DMPZ a, DMPZ b, DMPZexps_cmp cmp)
{
  DMPZ ans = NULL;
  DMPZ *tail = &ans;

  while (a != NULL && b != NULL)
  {
    const int comparison = cmp(a->exps, b->exps);
    if (comparison > 0)
    { /* case: LPP(a) > LPP(b) */
      *tail = a;
      tail = &((*tail)->next);
      a = a->next;
      continue;
    }
    if (comparison < 0)
    { /* case: LPP(b) > LPP(a) */
      *tail = b;
      tail = &((*tail)->next);
      b = b->next;
      continue;
    }
    /* case: LPP(a) = LPP(b), so add the terms */
    mpz_add(a->coeff, a->coeff, b->coeff);
    { /* delete leading term of b, and advance b */
      DMPZ junk = b;
      b = b->next;
      junk->next = NULL;
      DMPZdtor(junk);
    }
    if (mpz_sgn(a->coeff) == 0) /* leading terms cancel */
    { /* delete leading term of a, and advance a */
      DMPZ junk = a;
      a = a->next;
      junk->next = NULL;
      DMPZdtor(junk);
      continue;
    }
    /* leading term of a is now the sum, move it into the answer */
    *tail = a;
    tail = &((*tail)->next);
    a = a->next;
  }
  if (a == NULL) *tail = b;
  else *tail = a;
  return ans;
}


DMPZ DMPZsort(DMPZ f, DMPZexps_cmp cmp)
{
  DMPZ front, back, iter;
  int n;

  if (f == NULL || f->next == NULL) return f;
  front = f;
  iter = f;
  for (n = (DMPZ_nterms(f)-1)/2; n > 0; n--) iter = iter->next;
  back = iter->next;
  iter->next = NULL;
  
  return DMPZmerge(DMPZsort(front, cmp), DMPZsort(back, cmp), cmp);
}
 


DMPZ DMPZadd(const DMPZ a, const DMPZ b)
{
  DMPZ a_copy = DMPZsort(DMPZcopy(a), DMPZorder_lex);
  DMPZ b_copy = DMPZsort(DMPZcopy(b), DMPZorder_lex);

  return DMPZmerge(a_copy, b_copy, DMPZorder_lex);
}


void DMPZadd3(DMPZ* dest, const DMPZ a, const DMPZ b)
{
  DMPZ tmp = DMPZadd(a, b);
  DMPZdtor(*dest);
  *dest = tmp;
}


/***************************************************************************/

void DMPZprint(const DMPZ f)
{
  DMPZ iter;
  int i;

  if (f == NULL) { printf("0\n"); return; }
  for (iter=f; iter; iter = iter->next)
  {
    if (mpz_sgn(iter->coeff) > 0) printf("+");
    mpz_out_str(stdout,10,iter->coeff);
    for (i=0; i < NVARS; i++)
      if (iter->exps[i])
      {
        printf("*x(%d)", i);
        if (iter->exps[i] > 1) printf("^%d", iter->exps[i]);
      }
  }
  printf("\n");
}

