//   Copyright (c)  1997-2006,2008  John Abbott,  Anna M Bigatti

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
#include <math.h>

#include "jalloc.h"
#include "Zsolve.h"
#include "FFmat.h"
#include "Zmat.h"
#include "Zmat_det_bound.h"
#include "logi.h"
#include "FFsolve.h"
#include "FFdet.h"
#include "primes.h"
#include "WGD.h"
#include "mpz_from_padic.h"
#include "jaaerror.h"

/***************************************************************************/
/* This function gets two p-adic digits at a time; the lower order one is  */
/* put into "lo" the higher order one into "hi".                           */
/* The prime for the p-adic expansion is taken from CurrentFF.prime.       */
/* The matrix M is left unchanged.                                         */

static void Zmat_to_FFmat2(FFmat lo, FFmat hi, Zmat M)
{
  int i, j;
  FFelem p = CurrentFF.prime;
  FFelem p2 = p*p;

  if (M->nrows > lo->nrows || M->ncols > lo->ncols ||
      M->nrows > hi->nrows || M->ncols > hi->ncols)
    JERROR(JERROR_MATRIX);
  for (i=0; i < M->nrows; i++)
    for (j=0; j < M->ncols; j++)
    {
      lo->entry[i][j] = mpz_fdiv_ui(M->entry[i][j], p2);
      hi->entry[i][j] = lo->entry[i][j]/p;
      lo->entry[i][j] -= p * hi->entry[i][j];
    }
}


/***************************************************************************/
/* Function to multiply two matrices with FFelem entries, computing result */
/* modulo p^2 (instead of modulo p).                                       */

static void FFmat_mul2(FFmat product, FFmat A, FFmat B)
{
  int i, j, k, m, n, r;
  FFelem sum;
  FFelem p2 = CurrentFF.prime; p2=p2*p2;

  m = A->nrows;
  n = A->ncols;
  r = B->ncols;
  for (i=0; i < m; i++)
    for (j=0; j < r; j++)
    {
      sum = 0;
      for (k=0; k < n; k++)
      {
        sum += A->entry[i][k]*B->entry[k][j];
        if (sum >= p2) sum -= p2;
      }
      product->entry[i][j] = sum;
    }
}



/***************************************************************************/
/* Solve a non-singular linear system where both matrix and RHS are integer*/
/* Method used is Hensel lifting.                                          */
/* WARNING: the matrix M must be square and invertible!                    */

int Zsolve(Zmat soln, mpz_t *denom, Zmat M, Zmat rhs)
{
  int n, r, i, j, nsteps, target_height;
  FFelem p;
  mpz_t q, tmp, denom_bound, D;
  FFmat Mp, Mp1, rhsp, rhsp2, solnp, solnp2, solnpk;
  FFmat Mp_inverse;
  Zmat Msolnp, num, den;
  FF Fp = NULL; /* to keep compiler quiet */

  if (M->nrows != M->ncols ||
      M->ncols != rhs->nrows ||
      soln->nrows != rhs->nrows ||
      soln->ncols != rhs->ncols)
    JERROR(JERROR_MATRIX);

  n = M->nrows;
  r = rhs->ncols;
  if (r!=1) {fprintf(stderr,"FIX IMPLEMENTATION OF SOLNPK FOR r>1\n");exit(1);}

  Mp = FFmat_ctor(n, n);

/* The following block searches for a suitable prime p upon which to base */
/* the p-adic lifting.  If the matrix is singular we'll discover it here. */
  {
    PrimeSource PS = PrimeSourceCtor();
    double log_modulus = 0;
    double h = Zmat_det_bound(M);
    unsigned int rhs_max;
    FFelem detp=0; /* set to zero to keep compiler quiet */
    for (p = NextPrime(PS); log_modulus < h; p = NextPrime(PS))
    {
      Fp = FFctor(p);
      FFselect(Fp);
      Zmat_to_FFmat(Mp, M);
      detp = FFdet(Mp);
      log_modulus += logi(p);
      if (detp != 0) break;
      FFdtor(Fp);
    }
    PrimeSourceDtor(PS);

    if (detp == 0)  /* matrix was singular */
    {
      FFmat_dtor(Mp);
      return 0;
    }

    rhs_max = 0;
    for (i=0; i < n; i++)
      for(j=0; j < r; j++)
        if (mpz_sizeinbase(rhs->entry[i][j], 2) > rhs_max)
          rhs_max = mpz_sizeinbase(rhs->entry[i][j], 2);
    /* Formula below comes from Hadamard's bound on Cramer's numerator */
    target_height = 1+(int)((logi(2)+logi(2)*(rhs_max+1)+logi(n)/2+h)/logi(p));
  }

  /* Compute inverse modulo p -- put result in Mp_inverse */
  Mp_inverse = FFmat_ctor(n, n);
  {
    FFmat Mp_identity = FFmat_ctor_identity(n);
    int *junk = (int*)MALLOC(n*sizeof(int));
    for (i=0; i < n; ++i) junk[i]=0;    
    Zmat_to_FFmat(Mp, M); /* note:  Mp is trashed by call below to FFsolve */
    FFsolve(Mp_inverse, junk, Mp, Mp_identity);
    FREE(junk);
    FFmat_dtor(Mp_identity);
  }

  rhsp = FFmat_ctor(n, r);
  rhsp2 = FFmat_ctor(n, r);
  solnp = FFmat_ctor(n, r);
  solnp2 = FFmat_ctor(n, r);

  Msolnp = Zmat_ctor(n, r);
  num = Zmat_ctor(n, r);
  den = Zmat_ctor(n, r);
  for (i=0; i < n; i++)
    for (j=0; j < r; j++)
      mpz_set_ui(den->entry[i][j], 1);

  solnpk = FFmat_ctor(n, target_height);

  Mp1 = FFmat_ctor(n, n);
  Zmat_to_FFmat2(Mp, Mp1, M);

  for (nsteps = 0; nsteps < target_height; nsteps++)
  {
    Zmat_to_FFmat2(rhsp, rhsp2, rhs);
    FFmat_mul(solnp, Mp_inverse, rhsp);
    FFmat_mul2(rhsp, Mp, solnp);
    for (i=0; i < n; i++)
      for (j=0; j < r; j++)
        rhsp->entry[i][j] = FFsub(rhsp2->entry[i][j], rhsp->entry[i][j]/p);
    FFmat_mul(rhsp2, Mp1, solnp);
    for (i=0; i < n; i++)
      for (j=0; j < r; j++)
        rhsp->entry[i][j] = FFsub(rhsp->entry[i][j], rhsp2->entry[i][j]);

    FFmat_mul(solnp2, Mp_inverse, rhsp);
    for (i=0; i < n; i++)
      for (j=0; j < r; j++)
        solnp->entry[i][j] += p*solnp2->entry[i][j];

    Zmat_mul_FFmat(Msolnp, M, solnp);
    for (i=0; i < n; i++)
      for (j=0; j < r; j++)
      {
        solnpk->entry[i][nsteps] = solnp->entry[i][j];
        mpz_sub(rhs->entry[i][j], rhs->entry[i][j], Msolnp->entry[i][j]);
        mpz_fdiv_q_ui(rhs->entry[i][j], rhs->entry[i][j], p*p);
      }
  }

  FFdtor(Fp);

  mpz_init(tmp);
  mpz_init(q);
  mpz_ui_pow_ui(q, p*p, nsteps);
  mpz_init_set_ui(denom_bound,1);
  mpz_mul_2exp(denom_bound, denom_bound, mpz_sizeinbase(q,2)/2);
  mpz_init(D);
  for (j=0; j < r; j++)
  {
    mpz_set_ui(denom[j], 1);
    for (i=0; i < n; i++)
    {
      mpz_from_padic(tmp, p*p, nsteps, solnpk->entry[i]);
      mpz_mul(tmp, tmp, denom[j]);
      mpz_mod(tmp,tmp,q);
      mpz_set_ui(D,1);
      mpz_mul_2exp(D, D, mpz_sizeinbase(denom_bound,2)-mpz_sizeinbase(denom[j],2));
      modular_to_rational(num->entry[i][j], den->entry[i][j], D, tmp, q);
      mpz_mul(den->entry[i][j], den->entry[i][j], denom[j]);
      mpz_set(denom[j], den->entry[i][j]);
    }
    for (i=0; i < n; i++)
    {
      mpz_divexact(D, denom[j], den->entry[i][j]);
      mpz_mul(soln->entry[i][j], num->entry[i][j], D);
    }
  }
  mpz_clear(D);
  mpz_clear(denom_bound);
  mpz_clear(tmp);
  FFmat_dtor(solnpk);
  FFmat_dtor(solnp);
  FFmat_dtor(solnp2);
  FFmat_dtor(rhsp);
  FFmat_dtor(rhsp2);
  FFmat_dtor(Mp);
  FFmat_dtor(Mp1);
  FFmat_dtor(Mp_inverse);
  Zmat_dtor(Msolnp);
  
#if 0
  printf("ONLY NEED CHECK ANSWER IF WE USE EARLY ABORT -- not currently\n");
  {
    int M_biggest;

    M_biggest = 0;
    for (i=0; i < n; i++)
      for (j=0; j < n; j++)
        if (mpz_sizeinbase(M->entry[i][j], 2) > M_biggest)
          M_biggest = mpz_sizeinbase(M->entry[i][j], 2);
    
    for (j=0; j < r; j++)
    {
      S_biggest = 0;
      for (i=0; i < n; i++)
        if (mpz_sizeinbase(soln->entry[i][j], 2) > S_biggest)
          S_biggest = mpz_sizeinbase(soln->entry[i][j], 2);
      if (M_biggest + S_biggest + logi(n)/logi(2) > log2_modulus) break;
    }
  }
#endif

  mpz_clear(q);
  Zmat_dtor(num);
  Zmat_dtor(den);

  return 1;
}
