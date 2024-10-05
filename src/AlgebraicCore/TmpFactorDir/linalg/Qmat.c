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

#ifdef FACTOR_DEBUG
#include <stdio.h>
#endif
#include "Qmat.h"
#include "jalloc.h"
#include "jaaerror.h"
#include "mpq_to_FFelem.h"

Qmat Qmat_ctor(int nrows, int ncols)
{
  Qmat THIS;
  int i, j;

  THIS = (Qmat)MALLOC(sizeof(struct Qmat_struct));
  THIS->nrows = nrows;
  THIS->ncols = ncols;
  THIS->entry = (mpq_t**)MALLOC(nrows*sizeof(mpq_t*));
  for (i=0; i < nrows; i++)
  {
    THIS->entry[i] = (mpq_t*)MALLOC(ncols*sizeof(mpq_t));
    for (j=0; j < ncols; j++)
      mpq_init(THIS->entry[i][j]);
  }
  return THIS;
}

void Qmat_dtor(Qmat THIS)
{
  int i, j;

  for (i=0; i < THIS->nrows; i++)
  {
    for (j=0; j < THIS->ncols; j++)
      mpq_clear(THIS->entry[i][j]);
    FREE(THIS->entry[i]);
  }
  FREE(THIS->entry);
  FREE(THIS);
}


void Qmat_swap_rows(Qmat M, int i, int j)
{
  mpq_t *swap;
  
  if (M->nrows < i || M->nrows < j || i < 0 || j < 0) JERROR(JERROR_MATRIX);
  if (i == j) return;
  swap = M->entry[i];
  M->entry[i] = M->entry[j];
  M->entry[j] = swap;
}


void Qmat_swap_cols(Qmat M, int i, int j)
{
  int k;
  mpq_t swap;
  
  if (M->ncols < i || M->ncols < j || i < 0 || j < 0) JERROR(JERROR_MATRIX);
  if (i == j) return;
  mpq_init(swap);
  for (k=0; k < M->nrows; k++)
  {
    mpq_set(swap, M->entry[k][i]);
    mpq_set(M->entry[k][i], M->entry[k][j]);
    mpq_set(M->entry[k][j], swap);
  }
  mpq_clear(swap);
}


int Qmat_to_FFmat(FFmat Mp, Qmat M)
{
  int i, j;

  if (Mp->nrows != M->nrows || Mp->ncols != M->ncols) JERROR(JERROR_MATRIX);
  for (i=0; i < M->nrows; i++)
    for (j=0; j < M->ncols; j++)
      if (!mpq_to_FFelem(&Mp->entry[i][j], M->entry[i][j])) return 0;
  return 1;
}

#if 0
// OLD VERSION
int Qmat_to_FFmat(FFmat Mp, Qmat M)
{
  int i, j, OK=1;
  FFelem num, den, p = CurrentFF.prime;

  if (Mp->nrows != M->nrows || Mp->ncols != M->ncols) JERROR(JERROR_MATRIX);
  for (i=0; i < M->nrows; i++)
    for (j=0; j < M->ncols; j++)
      if (!mpq_to_FFelem(&Mp->entry[i][j], M->entry[i][j])) return 0;
    {
      den = mpz_fdiv_ui(mpq_denref(M->entry[i][j]), p);
      num = mpz_fdiv_ui(mpq_numref(M->entry[i][j]), p);
      if (den == 0) OK = 0;
      else Mp->entry[i][j] = FFdiv(num, den); /* never divides by zero */
    }
  return OK;
}
#endif


void Qmat_mul(Qmat product, Qmat A, Qmat B)
{
  int i, j, k, m, n, r;
  mpq_t tmp;

  if (A->ncols != B->nrows) JERROR(JERROR_MATRIX);
  m = A->nrows;
  n = A->ncols;
  r = B->ncols;
  mpq_init(tmp);
  for (i=0; i < m; i++)
    for (j=0; j < r; j++)
    {
      mpq_set_ui(product->entry[i][j], 0, 1);
      for (k=0; k < n; k++)
      {
	mpq_mul(tmp, A->entry[i][k], B->entry[k][j]);
	mpq_add(product->entry[i][j], product->entry[i][j], tmp);
      }
    }
  mpq_clear(tmp);
}

/***************************************************************************/
/* Matrix by vector product.  This code is simple rather than fast.        */
/* ALIASING between product and vec is NOT ALLOWED (and not checked).      */
/* Assumes that the vectors have enough space available.                   */

void Qmat_mul_vec(mpq_t *product, Qmat M, mpq_t *vec)
{
  int i, j;
  mpq_t tmp;

  mpq_init(tmp);
  for (i=0; i < M->nrows; i++)
  {
    mpq_set_ui(product[i], 0, 1);
    for (j=0; j < M->ncols; j++)
    {
      mpq_mul(tmp, M->entry[i][j], vec[j]);
      mpq_add(product[i], product[i], tmp);
    }
  }
  mpq_clear(tmp);
}


#ifdef FACTOR_DEBUG
/***************************************************************************/
/* Routine to print out a Qmat as an array of numbers separated by tabs.   */

void Qmat_print(Qmat M)
{
  int ncols = M->ncols;
  int nrows = M->nrows;
  int i, j;
  for (i=0; i < nrows; i++)
  {
    for (j=0; j < ncols; j++)
    {
      mpz_out_str(stdout, 10, mpq_numref(M->entry[i][j]));
      printf("/");
      mpz_out_str(stdout, 10, mpq_denref(M->entry[i][j]));
      printf("\t");
    }
    printf("\n");
  }
}
#endif
