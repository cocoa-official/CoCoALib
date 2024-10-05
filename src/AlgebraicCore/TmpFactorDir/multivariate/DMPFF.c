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

/* only needed for DMPFFprint -- for debugging really */
#include <stdio.h>


#include "DMPFF.h"


#include <stdlib.h>
#include "jalloc.h"
#include "jaaerror.h"

int DMPFF_NVARS = -1; /* silly initial value, should provoke errors if uninitialized */

/***************************************************************************/

/* This should be normal lexicographic ordering */
static int PPcmp(int *pp1, int *pp2)
{
  int var;

  for(var=0; var < DMPFF_NVARS-1; var++)
    if (pp1[var] != pp2[var]) break;
  return pp1[var] - pp2[var];
}


DMPFF DMPFFctor(int coeff, int *exps)
{
  DMPFF ans;

  ans = (DMPFF)MALLOC(sizeof(struct DMPFF_struct));
  ans->exps = exps;
  ans->next = NULL;
  return ans;
}

void DMPFFdtor(DMPFF f)
{
  DMPFF THIS, next;

  for (THIS = f; THIS; THIS = next)
  {
    next = THIS->next;
    FREE(THIS->exps);
    FREE(THIS);
  }
}


/* A useful variant of DMPFFctor. */
DMPFF DMPFFprepend(int coeff, int *exps, DMPFF f)
{
  DMPFF ans;

  ans = DMPFFctor(coeff, exps);
  ans->next = f;
  return ans;
}


/* Destructive reversal */
DMPFF DMPFFreverse(DMPFF f)
{
  DMPFF THIS, prev, next;

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

DMPFF DMPFFcopy(const DMPFF f)
{
  DMPFF ans, iter, tmp;
  int i, *exps_copy;

  ans = NULL;
  for (iter=f; iter; iter = iter->next)
  {
    exps_copy = (int*)MALLOC(DMPFF_NVARS*sizeof(int));
    for (i=0; i < DMPFF_NVARS; i++) exps_copy[i] = iter->exps[i];
    tmp = DMPFFctor(iter->coeff, exps_copy);
    tmp->next = ans;
    ans = tmp;
  }
  return DMPFFreverse(ans);
}


/***************************************************************************/

int DMPFF_nterms(const DMPFF f)
{
  int i;
  DMPFF iter;

  i = 0;
  for (iter=f; iter; iter=iter->next) i++;
  return i;
}


/***************************************************************************/
/* Makes a copy of that part of f for which the degree of var equals deg   */
/* and divide it by var^deg.  The coefficient of var^deg in f.             */
/*
DMPFF DMPFFcoeff(const DMPFF f, int var, int deg)
{
  DMPFF leading_term;
  int *exps_copy;
  int i;
  
  if (f == NULL) return NULL;
  if (f->exps[var] != deg) return DMPFFcoeff(f->next, var, deg);
  exps_copy = (int*)MALLOC(DMPFF_NVARS*sizeof(int));
  for (i=0; i < DMPFF_NVARS; i++) exps_copy[i] = f->exps[i];
  exps_copy[var] = 0;
  leading_term = DMPFFctor(f->coeff, exps_copy);
  leading_term->next = DMPFFcoeff(f->next, var, deg);
  return leading_term;
}
*/
DMPFF DMPFFcoeff(const DMPFF f, int var, int deg)
{
  DMPFF iter, ans;
  int *exps_copy;
  int i;

  ans = NULL;
  for (iter = f; iter; iter = iter->next)
  {
    if (iter->exps[var] != deg) continue;
    exps_copy = (int*)MALLOC(DMPFF_NVARS*sizeof(int));
    for (i=0; i < DMPFF_NVARS; i++) exps_copy[i] = f->exps[i];
    exps_copy[var] = 0;
    ans = DMPFFprepend(f->coeff, exps_copy, ans);    
  }
  return DMPFFreverse(ans);
}


/***************************************************************************/
/* Degree and multi-degree of a DMPFF.                                      */

void DMPFFdegs(int *degs, const DMPFF f)
{
  DMPFF iter;
  int var;

  for (var=0; var < DMPFF_NVARS; var++) degs[var] = 0;
  for (iter=f; iter; iter = iter->next)
    for (var=0; var < DMPFF_NVARS; var++)
      if (degs[var] < iter->exps[var]) degs[var] = iter->exps[var];
}


int DMPFFdeg(const DMPFF f, int var)
{
  DMPFF iter;
  int deg = 0;

  for (iter=f; iter; iter = iter->next)
    if (deg < iter->exps[var]) deg = iter->exps[var];
  return deg;
}


int DMPFFtotal_deg(const DMPFF f)
{
  int max, curr, i;
  DMPFF iter;
  
  max = -1;
  for (iter=f; iter; iter = iter->next)
  {
    curr = 0;
    for (i=0; i < DMPFF_NVARS; i++)
      curr += iter->exps[i];
    if (curr > max) max = curr;
  }
  return max;
}


void DMPFFmindegs(int *degs, const DMPFF f)
{
  DMPFF iter;
  int var;

  for (var=0; var < DMPFF_NVARS; var++) degs[var] = 0;
  if (f == NULL) return;
  for (var=0; var < DMPFF_NVARS; var++) degs[var] = f->exps[var];
  for (iter=f->next; iter; iter = iter->next)
    for (var=0; var < DMPFF_NVARS; var++)
      if (degs[var] > iter->exps[var]) degs[var] = iter->exps[var];
}



/***************************************************************************/
/* This function is used by DMPFFdiv_exact below.                           */

//static DMPFF DMPFFshift_sub(DMPFF f, DMPFF g, DMPFF monom)
//{
//  int var, *exps;
//  mpz_t coeff;
//  DMPFF f_iter, g_iter, ans;
//  int cmp, update_exps;
//  
//  ans = NULL;
//  exps = NULL; /* to keep the compiler quiet */
//  g_iter = g;
//  f_iter = f;
//  update_exps = 1;
//  while (f_iter || g_iter)
//  {
//    if (g_iter == NULL) goto step_f;
//    if (f_iter == NULL) goto step_g;
//    if (update_exps)
//    {
//      exps = (int*)malloc(DMPFF_NVARS*sizeof(int));
//      for (var=0; var < DMPFF_NVARS; var++)
//        exps[var] = g_iter->exps[var] + monom->exps[var];
//      update_exps = 0;
//    }
//    cmp = PPcmp(f_iter->exps, exps);
//    if (cmp > 0) goto step_f;
//    if (cmp < 0) goto step_g;
//    mpz_init(coeff);
//    mpz_mul(coeff, g_iter->coeff, monom->coeff);
//    mpz_sub(coeff, f->coeff, coeff);
//    ans = DMPFFprepend(coeff, exps, ans);
//    f_iter = f_iter->next;
//    g_iter = g_iter->next;
//    update_exps = 1;
//    continue;
//
//    step_f:
//    ans = DMPFFprepend(f->coeff,f->exps, ans);
//    f_iter = f_iter->next;
//    continue;
//
//    step_g:
//    mpz_init(coeff);
//    mpz_mul(coeff, g_iter->coeff, monom->coeff);
//    mpz_neg(coeff, coeff);
//    ans = DMPFFprepend(coeff, exps, ans);
//    g_iter = g_iter->next;
//    update_exps = 1;
//    continue;
//  }
//
//  return DMPFFreverse(ans);
//}

/***************************************************************************/
/* Exact division of one multivariate polynomial by another.               */
/* DESTROYS FIRST ARGUMENT!                                                */

//DMPFF DMPFFdiv_exact(DMPFF f, const DMPFF g)
//{
//  mpz_t q;
//  int *exps, var;
//  DMPFF quot, tmp;
//
//  if (g == NULL) JERROR(JERROR_DIV_BY_ZERO);
//  if (f == NULL) return NULL;
//  quot = NULL;
//  mpz_init(q);
//  while(f)
//  {
//    mpz_fdiv_q(q, f->coeff, g->coeff);
//    exps = (int*)MALLOC(DMPFF_NVARS*sizeof(int));
//    for (var=0; var < DMPFF_NVARS; var++) exps[var] = f->exps[var] - g->exps[var];
//    quot = DMPFFprepend(q, exps, quot);
//    tmp = DMPFFshift_sub(f, g, quot);
//    DMPFFdtor(f);
//    f = tmp;
//  }
//  mpz_clear(q);
//
//  return DMPFFreverse(quot);
//}

/***************************************************************************/


void DMPFFdiv2z(DMPFF f, int n)
{
  DMPFF iter;

  for (iter=f; iter; iter = iter->next)
    iter->coeff = FFdiv(iter->coeff, n);
}


void DMPFFprint(const DMPFF f)
{
  DMPFF iter;
  int i;

  if (f == NULL) { printf("0\n"); return; }
  for (iter=f; iter; iter = iter->next)
  {
    printf("+%d", iter->coeff);
    for (i=0; i < DMPFF_NVARS; i++)
      if (iter->exps[i])
      {
        printf("*x(%d)", i);
        if (iter->exps[i] > 1) printf("^%d", iter->exps[i]);
      }
  }
  printf("\n");
}
