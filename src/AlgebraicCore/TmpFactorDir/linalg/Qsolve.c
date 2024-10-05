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

#include "Qsolve.h"

#include <math.h>
#include "jalloc.h"
#include "logi.h"
#include "FF.h"
#include "FFmat.h"
#include "FFsolve.h"
#include "WGD.h"
#include "mpz_alias.h"
#include "mpz_cra_ui.h"
#include "primes.h"
#include "Zmat.h"

#if 0
/***************************************************************************/
THIS FUNCTION IS NOT USED CURRENTLY

static int FFelem_equal_mpq(FFelem r, mpq_t q)
{
  FFelem p = CurrentFF.prime;
  
  return FFmul(r, mpz_fdiv_ui(mpq_denref(q), p)) == mpz_fdiv_ui(mpq_numref(q), p);
}
#endif

/***************************************************************************/

struct Qsolver_struct
{
  Qmat M;                       /* Matrix nrows-by-ncols over Q          */
  Qmat rhs;                     /* Matrix ncols-by-r over Z              */

  FFelem p;                     /* Last prime used.                      */
  FFmat Mp;                     /* workspace, usu. M mod p               */
  Zmat res;                     /* Modular solutions.                    */
  mpz_t modulus;                /* Modulus of solns in res               */
  Qmat soln;                    /* Numerator/denominator of solns        */
  int checki, checkj;           /* Soln coord to watch (for end cond)    */
};

typedef struct Qsolver_struct *Qsolver;

/***************************************************************************/

static Qsolver Qsolver_ctor(Qmat soln, /*const*/ Qmat M, /*const*/ Qmat rhs)
{
  Qsolver ans;

  ans = (Qsolver)MALLOC(sizeof(struct Qsolver_struct));
  ans->soln = soln;
  ans->M = M;
  ans->rhs = rhs;

  ans->p = 0;
  mpz_init_set_ui(ans->modulus, 1);
  ans->Mp = FFmat_ctor(M->nrows, M->ncols);

  ans->res = Zmat_ctor(M->ncols, rhs->ncols);
  ans->checki = 0;
  ans->checkj = 0;

  return ans;
}


static void Qsolver_dtor(Qsolver THIS)
{
  mpz_clear(THIS->modulus);
  Zmat_dtor(THIS->res);
  FFmat_dtor(THIS->Mp);
  FREE(THIS);
}


/***************************************************************************/
/* Result is true iff we have already established that no solution for     */
/* vector r exists.                                                        */

static int no_soln(const Qsolver THIS, int r)
{
  return mpz_sgn(mpq_denref(THIS->soln->entry[0][r])) == 0;
}

static void tell_no_soln(Qsolver THIS, int r)
{
  mpz_set_ui(mpq_denref(THIS->soln->entry[0][r]), 0);
}

/***************************************************************************/
/* Check solution coords for "stability" without knowing determinant.      */
/* Result is non-zero if all coords are stable, otherwise zero.            */

static int Qsolve_check_coeffs1(Qsolver THIS, FFmat solnp)
{
  FFelem p = THIS->p;
  int i, j, ncols, r, max_num;
  int lcd_bound, num_bound;
  Qmat soln = THIS->soln; /* alias */
  mpz_t soln_lcd, tmp, tmp2;

  if (mpz_sizeinbase(THIS->modulus, 2) < 40) return 0;
  lcd_bound = mpz_sizeinbase(THIS->modulus, 2)/2 - 20;
  num_bound = lcd_bound;
//printf("Starting stability check:\n");
//printf("log2(modulus) = %d\n", mpz_sizeinbase(THIS->modulus, 2));
//printf("lcd and num bounds are %d\n", num_bound);
/*
  i = THIS->checki;
  j = THIS->checkj;
  modular_to_rational(mpq_numref(soln->entry[i][j]), mpq_denref(soln->entry[i][j]), D, THIS->res->entry[i][j], THIS->modulus);
  if (!FFelem_equal_mpq(solnp->entry[i][j], soln->entry[i][j]))
  {
    mpz_clear(D);
    return 0;
  }
*/
  
  /* The coord we were watching is stable, so check the others... */

  mpz_init(soln_lcd);
  mpz_init(tmp);
  mpz_init(tmp2);
  ncols = THIS->M->ncols;
  r = THIS->rhs->ncols;
  for (j=0; j < r; j++)
  {
    if (no_soln(THIS, j)) continue;
    if (solnp->entry[0][j] == p)
    {
      mpz_set_ui(mpq_denref(soln->entry[0][j]), 0);
      continue;
    }
    mpz_set_ui(soln_lcd, 1);
    max_num = -999999; /* should really be minus infinity */
    for (i=0; i < ncols; i++)
    {
      mpz_alias num = mpq_numref(soln->entry[i][j]);
      mpz_alias den = mpq_denref(soln->entry[i][j]);
      int den_bound = lcd_bound - mpz_sizeinbase(soln_lcd, 2) - 1;
      int size;

      if (mpz_sgn(THIS->res->entry[i][j]) == 0)
      {
        mpz_set_ui(num, 0);
        mpz_set_ui(den, 1);
        continue;
      }
      mpz_mul(tmp, soln_lcd, THIS->res->entry[i][j]);
      modular_to_rational2(num, den, den_bound, tmp, THIS->modulus);
//printf("Soln has size %d/%d\n", mpz_sizeinbase(num, 2), mpz_sizeinbase(soln_lcd, 2)+mpz_sizeinbase(den, 2));
      if (2 + (int)mpz_sizeinbase(num, 2) > num_bound) goto failed;
      size = 1 + mpz_sizeinbase(num, 2) - mpz_sizeinbase(den, 2) - mpz_sizeinbase(soln_lcd, 2);
      if (size > max_num) max_num = size;
      mpz_gcd(tmp, num, soln_lcd);
      mpz_divexact(num, num, tmp);
      mpz_divexact(tmp2, soln_lcd, tmp); /* GMP bug inhibits aliasing between tmp and tmp2 */
      mpz_mul(soln_lcd, soln_lcd, den);
      mpz_mul(den, den, tmp2);
//printf("New lcd has size %d, bound is %d\n", mpz_sizeinbase(soln_lcd, 2) , lcd_bound);      
      if ((int)mpz_sizeinbase(soln_lcd, 2) < lcd_bound &&
          max_num + (int)mpz_sizeinbase(soln_lcd, 2) < num_bound) continue;

      failed:
      /* Mark this coord for future watching, clean up, and return 0 */
      THIS->checki = i;
      THIS->checkj = j;
      mpz_clear(soln_lcd);
      mpz_clear(tmp);
      mpz_clear(tmp2);
      return 0;
    }
  }
  mpz_clear(soln_lcd);
  mpz_clear(tmp);
  mpz_clear(tmp2);
  return 1;
}



/***************************************************************************/
/* This function is called after the linear system has been solved in some */
/* new finite field.  It does a spot check on a particular solution coord  */
/* and if that did not change, it then checks all the other solution coords*/
/* and if they did not change this last iteration then it tries to prove the*/
/* solution obtained is correct.  If so return it, otherwise lift further. */

/* Result is 0 if more "lifting" is needed; otherwise non-zero.            */
static int Qsolve_check(Qsolver THIS, FFmat solnp)
{
  if (mpz_cmp_ui(THIS->modulus, 1) == 0) return 0; /* do nothing first time */
  if (!Qsolve_check_coeffs1(THIS, solnp)) return 0;
  return 1; /* We have proved the answer correct. */
}


/***************************************************************************/
/* Compute "height" of a Qmat.                                             */

static int Qmat_height(Qmat M)
{
  int i, j, M_height, entry_height;

  M_height = 0;
  for (i=0; i < M->nrows; i++)
  {
    for (j=0; j < M->ncols; j++)
    {
      if (mpz_sgn(mpq_numref(M->entry[i][j])) == 0) continue;
      entry_height = mpz_sizeinbase(mpq_numref(M->entry[i][j]), 2) +
                     mpz_sizeinbase(mpq_denref(M->entry[i][j]), 2);
      if (entry_height > M_height) M_height = entry_height;
    }
  }
  return M_height;
}


/***************************************************************************/
/* Compute a useful bound for a Qmat.                                      */
/* The bound is max of logarithms of #columns times infinity row norm of M */
/* times least common denominator of row.  This bound gives an idea of how */
/* much numerators can grow in a matrix-by-vector product with M.          */

static double Qsolve_bound(Qmat M)
{
  int i, j, bound, row_max, log2_Mij;
  mpz_t row_lcd;

  bound = 0;
  mpz_init(row_lcd);
  for (i=0; i < M->nrows; i++)
  {
    /* Compute least common denom for row i in row_lcd */
    mpz_set_ui(row_lcd, 1);
    for (j=0; j < M->ncols; j++)
      if (mpz_sgn(mpq_numref(M->entry[i][j])) != 0)
        mpz_lcm(row_lcd, row_lcd, mpq_denref(M->entry[i][j]));

    /* Compute in row_max the least power of two which exceeds the largest */
    /* numerator when row is written using the least common denominator.   */
    row_max = 0;
    for (j=0; j < M->ncols; j++)
      if (mpz_sgn(mpq_numref(M->entry[i][j])) != 0)
      {
        log2_Mij = 1 + mpz_sizeinbase(mpq_numref(M->entry[i][j]), 2)
                     - mpz_sizeinbase(mpq_denref(M->entry[i][j]), 2);
        if (log2_Mij > row_max) row_max = log2_Mij;
      }
    row_max += 1 + mpz_sizeinbase(row_lcd, 2);
    if (row_max > bound) bound = row_max;
  }
  mpz_clear(row_lcd);
  return logi(M->ncols) + logi(2)*bound;
}


/***************************************************************************/
/* Solve a system of linear equations over the rationals.                  */
/* The right hand side contains one or more vectors.                       */

/* Input: solns space for the solutions (i.e. a matrix ncols by r).        */
/*        nrows is the number of rows M has.                               */
/*        ncols is the number of columns M has.                            */
/*        M is the "left hand side" matrix.                                */
/*        r is the number of vectors on the right hand side.               */
/*        rhs is an nrows-by-r matrix whose columns are these vectors.     */
/* Output: the return value is the rank of M.                              */
/*        The solutions are written into soln; a zero denominator of the   */
/*        first coordinate of a solution means that no solution exists.    */

/* The method is Chinese Remaindering with "early detection".              */
/* Really most of the work (and clever bits) are in Qsolve_check above.    */
/* No attempt is made to look for block structure.                         */


int Qsolve(Qmat soln, Qmat M, Qmat rhs)
{
  int i, j, update;
  int *shape;
  FFelem p;
  FFmat Mp, rhsp, solnp;
  Qsolver system;
  FF FFp;
  int min_bits, step_size, count_down;
  PrimeSource PS = PrimeSourceCtor();
  int r = rhs->ncols;
  double bound, log_modulus;

  bound = Qsolve_bound(M);
  log_modulus = 0;
  /* Find minimum height to lift to = P*Q where P/Q is in M or rhs */
  /* min_bits = max(Qmat_height(M), Qmat_height(rhs)); */
  {int tmp=Qmat_height(M);min_bits=Qmat_height(rhs);if(tmp>min_bits)min_bits=tmp;}

  system = Qsolver_ctor(soln, M, rhs);

  Mp = system->Mp;
  rhsp = FFmat_ctor(rhs->nrows, rhs->ncols);
  solnp = FFmat_ctor(M->ncols, rhs->ncols);

  shape = (int*)MALLOC(M->ncols*sizeof(int));
  for (i=0; i < M->ncols; i++) shape[i] = 0;
  step_size = 3;
  count_down = 3;
  while(1)
  {
    p = NextPrime(PS);
    FFp = FFctor(p);
    FFselect(FFp);
    /* Skip unsuitable primes */
    if (Qmat_to_FFmat(Mp, M) == 0 ||
        Qmat_to_FFmat(rhsp, rhs) == 0)
    { FFdtor(FFp); continue; }

    update = FFsolve(solnp, shape, Mp, rhsp);
    if (update & FFsolve_new_pivot)
    {
      mpz_set_ui(system->modulus, 1);
      for (i=0; i < M->ncols; i++)
        for (j=0; j < r; j++)
        {
          mpz_set_ui(system->res->entry[i][j], 0);
          mpq_set_ui(soln->entry[i][j], 0, 1);
        }
      step_size = 3;
      count_down = 3;
    }
    if (update & FFsolve_invalid) { FFdtor(FFp); continue; }
    
    /* We found a solution mod p.  Check if we're done; if not, combine   */
    /* with earlier solutions in a Chinese Remaindering "lift", and loop. */
    system->p = p;

    for (j=0; j < r; j++)
    {
      if (no_soln(system, j)) continue;
      if (solnp->entry[0][j] == p) { tell_no_soln(system, j); continue; }
      for (i=0; i < M->ncols; i++)
	mpz_cra_ui(system->res->entry[i][j], system->modulus, solnp->entry[i][j], p);
    }
    mpz_mul_ui(system->modulus, system->modulus, p);
    log_modulus += logi(p);
    FFdtor(FFp);
    count_down--;
    if (count_down == 0) count_down = ++step_size;
    if (count_down == step_size &&
        log_modulus > 2*bound + 40 &&
        Qsolve_check(system, solnp)) break;
  }
  
  /* When we get here, we have found and verified the solution. */
  PrimeSourceDtor(PS);
  Qsolver_dtor(system);
  FFmat_dtor(rhsp);
  FFmat_dtor(solnp);

  /* Compute the rank */
  r = 0;
  for (i=0; i < M->ncols; i++) r += shape[i];
  FREE(shape);

  return r;
}


