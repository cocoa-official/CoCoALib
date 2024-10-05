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

#include "FFmat.h"
#include "jalloc.h"
#include "jaaerror.h"

#include <stddef.h>

FFmat FFmat_ctor(int nrows, int ncols)
{
  FFmat THIS;
  int i;

  if (nrows < 1 || ncols < 1) return NULL;
  THIS = (FFmat)MALLOC(sizeof(struct FFmat_struct));
  THIS->nrows = nrows;
  THIS->ncols = ncols;
  THIS->entry = (FFelem**)MALLOC(nrows*sizeof(FFelem*));
  for (i=0; i < nrows; i++)
    THIS->entry[i] = (FFelem*)MALLOC(ncols*sizeof(FFelem));
  return THIS;
}


FFmat FFmat_ctor_identity(int n)
{
  int i, j;
  FFmat THIS = FFmat_ctor(n, n);

  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      THIS->entry[i][j] = (i == j);
  return THIS;
}


void FFmat_dtor(FFmat THIS)
{
  int i;
  for (i=0; i < THIS->nrows; i++) FREE(THIS->entry[i]);
  FREE(THIS->entry);
  FREE(THIS);
}


void FFmat_mul(FFmat product, FFmat Amat, FFmat Bmat)
{
  int i, j, k, stop;
  int m, n, r;
  FFelem **A = Amat->entry;
  FFelem **B = Bmat->entry;
  FFelem sum;
  FFelem p = CurrentFF.prime;
  FFelem limit = CurrentFF.k;
  FFelem shift = CurrentFF.shift;

  if (Amat->ncols != Bmat->nrows) JERROR(JERROR_MATRIX);
  m = Amat->nrows;
  n = Amat->ncols;
  r = Bmat->ncols;
  for (i=0; i < m; i++)
    for (j=0; j < r; j++)
    {
      sum = 0;
      for (k=0, stop=limit; stop < n; stop += limit)
      {
        for (; k < stop; k++)
          sum += A[i][k]*B[k][j];
        if (sum >= shift) sum -= shift;
      }
      for (; k < n; k++)
        sum += A[i][k]*B[k][j];
      product->entry[i][j] = sum%p;
    }
}

void FFmat_to_Zmat(Zmat A, FFmat Amodp)
{
  int i, j;
  int nrows = A->nrows;
  int ncols = A->ncols;

  if (nrows != Amodp->nrows || ncols != Amodp->ncols) JERROR(JERROR_MATRIX);
  for (i=0; i < nrows; i++)
    for (j=0; j < ncols; j++)
      mpz_set_ui(A->entry[i][j], (unsigned long)Amodp->entry[i][j]);
}


/***************************************************************************/
/* Matrix by vector product.                                               */
/* ALIASING between product and vec is NOT ALLOWED (and not checked for).  */
/* This is a stripped out copy of FFmat_mul specialized for vectors.       */

void FFmat_mul_vec(FFelem *product, FFmat Mmat, FFelem *vec)
{
  int i, k, stop;
  int m, n;
  FFelem **M = Mmat->entry;
  FFelem sum;
  FFelem p = CurrentFF.prime;
  FFelem limit = CurrentFF.k;
  FFelem shift = CurrentFF.shift;

  m = Mmat->nrows;
  n = Mmat->ncols;
  for (i=0; i < m; i++)
  {
    sum = 0;
    for (k=0, stop=limit; stop < n; stop += limit)
    {
      for (; k < stop; k++)
        sum += M[i][k]*vec[k];
      if (sum >= shift) sum -= shift;
    }
    for (; k < n; k++)
      sum += M[i][k]*vec[k];
    product[i] = sum%p;
  }
}

