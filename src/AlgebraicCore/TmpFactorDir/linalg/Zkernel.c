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

#include "Zkernel.h"
#include "jalloc.h"
#include "mpz_cmp_abs.h"

/* We use the function rand(). */
#include <stdlib.h>

#if 1
/***************************************************************************/
/* Find pivot row: we pick randomly a non-zero element from the column     */
/* avoiding the (first) element of maximal magnitude unless there is only  */
/* one non-zero element in which case we pick that.  If all entries are    */
/* zero, we return -1.  This strategy seems often to produce kernel bases  */
/* with quite short generators -- also you can call the entry function     */
/* several times and take whichever basis you like best as usually the     */
/* outputs will be different.  It is a bit slower than some methods :-(    */

static int pivot_row(Zmat M, int col, int *row_used)
{
  int i, pivot, max_row, count;
  int nrows = M->nrows;
  
  max_row = -1;
  count = 0;
  for (i=0; i < nrows; i++)
  {
    if (row_used[i] || mpz_sgn(M->entry[i][col]) == 0) continue;
    count++;
    if (max_row == -1 || mpz_cmp_abs(M->entry[max_row][col], M->entry[i][col]) < 0)
      max_row = i;
  }
  if (count < 2) return max_row; /* this is -1 if all entries are 0 */

  /* pick randomly an unused row with non-zero entry, avoiding max entry */
  count = 1 + rand()%(count-1);
  for (pivot=0; pivot < nrows; pivot++)
  {
    if (row_used[pivot] || mpz_sgn(M->entry[pivot][col]) == 0 || pivot == max_row)
      continue;
    if (--count == 0) break;
  }
  return pivot;
}
#endif

/* The following alternative implementation is simpler and probably faster */
/* but some tests show that the kernel basis vectors produced are longer.  */
/* This alternative is hidden using pre-processor directives.              */
#if 0
/***************************************************************************/
/* Find pivot row: the one having non-zero entry of smallest magnitude in  */
/* the specified column.  Havas hints at using this strategy.              */
static int pivot_row(Zmat M, int col, int *row_used)
{
  int i, pivot;
  int nrows = M->nrows;
  
  pivot = -1;
  for (i=0; i < nrows; i++)
  {
    if (row_used[i] || mpz_sgn(M->entry[i][col]) == 0) continue;
    if (pivot >= 0 && mpz_cmp_abs(M->entry[i][col], M->entry[pivot][col]) >= 0) continue;
    pivot = i;
  }
  return pivot;
}
#endif

/* The following alternative implementation is simpler and probably faster */
/* but some tests show that the kernel basis vectors produced are longer.  */
/* This alternative is hidden using pre-processor directives.              */
#if 0
/***************************************************************************/
/* Find pivot row: the one having second biggest entry  (in magnitude) in  */
/* the specified column.  Havas hints at using this strategy.              */
static int pivot_row(Zmat M, int col, int *row_used)
{
  int i, biggest, biggest2;
  int nrows = M->nrows;
  
  biggest = -1;
  biggest2 = -1;
  for (i=0; i < nrows; i++)
  {
    if (row_used[i] || mpz_sgn(M->entry[i][col]) == 0) continue;
    if (biggest < 0) { biggest = i; continue; }
    if (mpz_cmp_abs(M->entry[i][col], M->entry[biggest][col]) >= 0)
    {
      biggest2 = biggest;
      biggest = i;
      continue;
    }
    if (biggest2 < 0) { biggest2 = i; continue; }
    if (mpz_cmp_abs(M->entry[i][col], M->entry[biggest2][col]) <= 0) continue;
    biggest2 = i;
  }
  if (biggest2 >= 0) return biggest2;
  return biggest;
}
#endif

/***************************************************************************/
/* Value returned is dimension of the kernel.                              */

int Zmat_left_kernel_Zbasis(Zmat *Zbasis, Zmat M)
{
  int i, j, pivot, col, rows_remaining;
  int finished;
  int *row_used;
  mpz_t q, tmp;
  Zmat U;
  int nrows = M->nrows;
  int ncols = M->ncols;
  
  mpz_init(q);
  mpz_init(tmp);
  U = Zmat_ctor(nrows, nrows);
  for (i=0; i < nrows; i++) mpz_set_ui(U->entry[i][i], 1);
  row_used = (int*)MALLOC(nrows*sizeof(int));
  for (i=0; i < nrows; i++) row_used[i] = 0;
  rows_remaining = nrows;
  for (col=0; col < ncols; col++)
  {
    if (rows_remaining == 0) break;
    pivot = pivot_row(M, col, row_used);
    if (pivot < 0) continue;
    do
    {
      finished = 1;

      if (mpz_sgn(M->entry[pivot][col]) < 0)
      {
        for (j=col; j < ncols; j++) mpz_neg(M->entry[pivot][j], M->entry[pivot][j]);
        for (j=0;   j < nrows; j++) mpz_neg(U->entry[pivot][j], U->entry[pivot][j]);
      }
      finished = 1;
      for (i=0; i < nrows; i++)
      {
        if (row_used[i] || i == pivot) continue;
        mpz_fdiv_qr(q, M->entry[i][col], M->entry[i][col], M->entry[pivot][col]);
        if (mpz_sgn(q) == 0) continue;
        finished = 0;
        for (j=0; j < nrows; j++)
        {
          mpz_mul(tmp, q, U->entry[pivot][j]);
          mpz_sub(U->entry[i][j], U->entry[i][j], tmp);
        }
        for (j=col+1; j < ncols; j++)
        {
          mpz_mul(tmp, q, M->entry[pivot][j]);
          mpz_sub(M->entry[i][j], M->entry[i][j], tmp);
        }
      }
      pivot = pivot_row(M, col, row_used);
    } while (!finished);
    row_used[pivot] = 1;
    rows_remaining -= 1;
  }

  /* At this point dimension of kernel is rows_remaining */
  /* Move the Z-basis of the kernel into the first rows of U */
  j = 0;
  for (i=0; i < rows_remaining; i++, j++)
  {
    while (row_used[j]) j++;
    Zmat_swap_rows(U, i, j);
  }

  mpz_clear(q);
  mpz_clear(tmp);
  FREE(row_used);
  *Zbasis = U;
  return rows_remaining;
}
