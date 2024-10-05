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
#include "Zdet.h"

#include <math.h>
#include "jalloc.h"
#include "jaaerror.h"
#include "logi.h"
#include "FFmat.h"
#include "FFdet.h"
#include "mpz_cra_ui.h"
#include "mpz_log.h"
#include "primes.h"
#include "Zmat_det_bound.h"
#include "Zsolve.h"

/***************************************************************************/
/* Compute determinant of a (square) Zmat.  Result is placed in det.       */
/* Method used is Chinese Remaindering -- this is a bit limiting on a      */
/* 32-bit machine; NTL has a better scheme in that case.                   */

static void Zdet_cra(mpz_t det, Zmat M)
{
  int n;
  FFelem p, detp;
  FFmat Mp;
  double log_det_modulus, det_bound;
  mpz_t det_modulus;
  FF FFp;
  PrimeSource PS = PrimeSourceCtor();

  n = M->nrows;
  Mp = FFmat_ctor(n, n);

  det_bound = 1 + Zmat_det_bound(M); /* add 1 to allow for approx error */

  log_det_modulus = 0;
  mpz_init_set_ui(det_modulus, 1);
  mpz_set_ui(det, 0);
  while (log_det_modulus < det_bound)
  {
    p = NextPrime(PS);
    FFp = FFctor(p);
    FFselect(FFp);
    Zmat_to_FFmat(Mp, M);
    detp = FFdet(Mp);
    mpz_cra_ui(det, det_modulus, detp, p); /* ignore rtn value */
    mpz_mul_ui(det_modulus, det_modulus, (unsigned long)p);
    log_det_modulus += logi(p);
    FFdtor(FFp);
  }
  mpz_clear(det_modulus);
  PrimeSourceDtor(PS);
  FFmat_dtor(Mp);
}


/***************************************************************************/
/* Compute determinant of a (square) Zmat M.  Result is placed in det.     */
/* Method used is that published in ISSAC'99 (Abbott, Bronstein, Mulders). */
/* Strictly this implementation is "incomplete"; so in the rare cases when */
/* the missing bit of code is needed we fall back on Zdet_cra rather than  */
/* doing it properly -- this saves lots of implementation hassle and will  */
/* probably never be felt in practice.                                     */

void Zdet(mpz_t det, Zmat M)
{
  Zmat rhs, soln;
  mpz_t *denom, D, det_factor;
  int i, j, r, singular;
  int row, col;
  int n = M->nrows;

  if (M->nrows != M->ncols) JERROR(JERROR_MATRIX);
  if (n == 1) { mpz_set(det, M->entry[0][0]); return; }
  {
    /* If det bound is large compared to matrix size, use CRA, but     */
    /* if det bound is huge, cannot use CRA as we'd exhaust the primes */
    double hadamard = Zmat_det_bound(M);
    if (hadamard < 0.99*MAX_PRIME && 3000*n < hadamard) /* 3000 is HEURISTIC */
    {
      Zdet_cra(det, M);
      return;
    }
  }

  /* "Eliminate" random row and column by swapping with last row/col */
  /* We will swap them back at the end. */
  row = rand()%n;
  col = rand()%n;
  Zmat_swap_rows(M, row, n-1);
  Zmat_swap_cols(M, col, n-1);

  /* We are going to solve an (n-1)x(n-1) system with r=1 rhs vector */
  r = 1;
  rhs = Zmat_ctor(n-1, r);
  soln = Zmat_ctor(n-1, r);
  denom = (mpz_t*)MALLOC(r*sizeof(mpz_t));
  for(j=0; j < r; j++) mpz_init(denom[j]);
  for(i=0; i < n-1; i++) mpz_set(rhs->entry[i][0], M->entry[i][n-1]);

  {
    /* In this block we fake the size of M for the call to Zsolve */
    M->nrows=n-1;M->ncols=n-1; /* fake size as (n-1)x(n-1) */
    singular = Zsolve(soln, denom, M, rhs);
    M->nrows=n;M->ncols=n;     /* restore original size    */
  }
  if (singular == 0)
  {
    /* MISSING CODE -- WIMP OUT PATCH USED INSTEAD */
    /* This is where the "missing" piece of code should be. */
    Zmat_dtor(soln); Zmat_dtor(rhs);
    for(j=0; j < r; j++) mpz_clear(denom[j]);
    FREE(denom);
    Zmat_swap_rows(M, row, n-1);
    Zmat_swap_cols(M, col, n-1);
    Zdet_cra(det, M);
    return;
  }

  mpz_init_set(D, denom[0]);
  mpz_init(det_factor);
  mpz_mul(det_factor, M->entry[n-1][n-1], denom[0]);
  { /* use a block to limit the scope of prod */
    mpz_t prod;
    mpz_init(prod);
    for (i=0; i < n-1; i++)
    {
      mpz_mul(prod, soln->entry[i][0], M->entry[n-1][i]);
      mpz_sub(det_factor, det_factor, prod);
    }
    mpz_clear(prod);
  }
  /* free rhs, soln, and denom */
  for(j=0; j < r; j++) mpz_clear(denom[j]);
  FREE(denom);
  Zmat_dtor(rhs);
  Zmat_dtor(soln);

  mpz_set_ui(det, 0);
  if (mpz_sgn(det_factor) == 0) goto tidy_up;

  {
    PrimeSource PS = PrimeSourceCtor();
    FFmat Mp = FFmat_ctor(n-1, n-1);
    FFelem p, detp, Dp;
    FF FFp;
    mpz_t det_modulus;
    double log_det_modulus, det_bound;

    det_bound = 1 + Zmat_det_bound2(M, n-1); /* add 1 for safety */
    det_bound -= mpz_log(D);

    mpz_init_set_ui(det_modulus, 1);
    log_det_modulus = 0;

    while (log_det_modulus < det_bound)
    {
      p = NextPrime(PS);
      Dp = mpz_fdiv_ui(D, p);
      if (Dp == 0) continue; /* prime is "bad" so skip it */
      FFp = FFctor(p);
      FFselect(FFp);
      for (i=0; i < n-1; i++)
        for (j=0; j < n-1; j++)
          Mp->entry[i][j] = mpz_fdiv_ui(M->entry[i][j], (unsigned long)p);
      detp = FFdet(Mp);
      detp = FFdiv(detp, Dp);
      mpz_cra_ui(det, det_modulus, detp, p); /* ignore rtn value */
      mpz_mul_ui(det_modulus, det_modulus, (unsigned long)p);
      log_det_modulus += logi((unsigned long)p);
      FFdtor(FFp);
    }
    mpz_clear(det_modulus);
    PrimeSourceDtor(PS);
    FFmat_dtor(Mp);
  }
  mpz_mul(det, det, det_factor);

  /* juggle sign in case a row/column swap was not really a swap */
  if ((row == n-1) ^ (col == n-1)) mpz_neg(det, det);
tidy_up:
  mpz_clear(det_factor);
  mpz_clear(D);
  /* Put matrix M back as it was originally */
  Zmat_swap_rows(M, row, n-1);
  Zmat_swap_cols(M, col, n-1);
}


