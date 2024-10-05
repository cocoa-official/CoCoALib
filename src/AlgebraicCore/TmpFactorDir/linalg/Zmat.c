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

#include "Zmat.h"
#include "jalloc.h"
#include "jaaerror.h"

Zmat Zmat_ctor(int nrows, int ncols)
{
  Zmat THIS;
  int i, j;

  THIS = (Zmat)MALLOC(sizeof(struct Zmat_struct));
  THIS->nrows = nrows;
  THIS->ncols = ncols;
  THIS->entry = (mpz_t**)MALLOC(nrows*sizeof(mpz_t*));
  for (i=0; i < nrows; i++)
  {
    THIS->entry[i] = (mpz_t*)MALLOC(ncols*sizeof(mpz_t));
    for (j=0; j < ncols; j++)
      mpz_init(THIS->entry[i][j]);
  }
  return THIS;
}

void Zmat_dtor(Zmat THIS)
{
  int i, j;

  for (i=0; i < THIS->nrows; i++)
  {
    for (j=0; j < THIS->ncols; j++)
      mpz_clear(THIS->entry[i][j]);
    FREE(THIS->entry[i]);
  }
  FREE(THIS->entry);
  FREE(THIS);
}


void Zmat_swap_rows(Zmat M, int i, int j)
{
  mpz_t *swap;
  
  if (M->nrows < i || M->nrows < j || i < 0 || j < 0) JERROR(JERROR_MATRIX);
  if (i == j) return;
  swap = M->entry[i];
  M->entry[i] = M->entry[j];
  M->entry[j] = swap;
}


void Zmat_swap_cols(Zmat M, int i, int j)
{
  int k;
  mpz_t swap;
  
  if (M->ncols < i || M->ncols < j || i < 0 || j < 0) JERROR(JERROR_MATRIX);
  if (i == j) return;
  mpz_init(swap);
  for (k=0; k < M->nrows; k++)
  {
    mpz_set(swap, M->entry[k][i]);
    mpz_set(M->entry[k][i], M->entry[k][j]);
    mpz_set(M->entry[k][j], swap);
  }
  mpz_clear(swap);
}


void Zmat_to_FFmat(FFmat Mp, Zmat M)
{
  int i, j;
  FFelem p;

  if (Mp->nrows != M->nrows || Mp->ncols != M->ncols) JERROR(JERROR_MATRIX);
  p = CurrentFF.prime;
  for (i=0; i < M->nrows; i++)
    for (j=0; j < M->ncols; j++)
      Mp->entry[i][j] = mpz_fdiv_ui(M->entry[i][j], p);
}



void Zmat_mul(Zmat product, Zmat A, Zmat B)
{
  int i, j, k, m, n, r;
  mpz_t tmp;

  if (A->ncols != B->nrows) JERROR(JERROR_MATRIX);
  m = A->nrows;
  n = A->ncols;
  r = B->ncols;
  mpz_init(tmp);
  for (i=0; i < m; i++)
    for (j=0; j < r; j++)
    {
      mpz_set_ui(product->entry[i][j], 0);
      for (k=0; k < n; k++)
      {
	mpz_mul(tmp, A->entry[i][k], B->entry[k][j]);
	mpz_add(product->entry[i][j], product->entry[i][j], tmp);
      }
    }
  mpz_clear(tmp);
}


void Zmat_mul_FFmat(Zmat product, Zmat A, FFmat B)
{
  int i, j, k, m, n, r;
  mpz_t tmp;

  if (A->ncols != B->nrows ||
      product->nrows != A->nrows ||
      product->ncols != B->ncols)
    JERROR(JERROR_MATRIX);

  m = A->nrows;
  n = A->ncols;
  r = B->ncols;
  mpz_init(tmp);
  for (i=0; i < m; i++)
    for (j=0; j < r; j++)
    {
      mpz_set_ui(product->entry[i][j], 0);
      for (k=0; k < n; k++)
      {
	mpz_mul_ui(tmp, A->entry[i][k], (unsigned long)B->entry[k][j]);
	mpz_add(product->entry[i][j], product->entry[i][j], tmp);
      }
    }
  mpz_clear(tmp);
}

