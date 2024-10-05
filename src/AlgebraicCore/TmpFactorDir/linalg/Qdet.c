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
#include "Qdet.h"

#include "Zdet.h"

// ***********************************************
// *** !!! THE FN BELOW WAS NEVER FINISHED !!! ***
// ***********************************************
// (so we hide it using a CPP conditional)
#if 0

static void compute_row_lcm(mpz_t *row_lcm, Qmat M)
{
  int i, j, n;
  n = M->nrows;
  for (i=0; i < n; i++)
  {
    mpz_set_ui(row_lcm[i], 0);
    for (j=0; j < n; j++)
      mpz_lcm(row_lcm[i], row_lcm[i], mpq_denref(M->entry[i][j]));
  }
}

static void compute_col_lcm(mpz_t *col_lcm, Qmat M)
{
  int i, j, n;
  n = M->nrows;
  for (j=0; j < n; j++)
  {
    mpz_set_ui(col_lcm[j], 0);
    for (i=0; i < n; i++)
      mpz_lcm(col_lcm[j], col_lcm[j], mpq_denref(M->entry[i][j]));
  }
}


static int choose_row_or_col(Qmat M, mpz_t *row_lcm, mpz_t *col_lcm)
{
  int n, i, j, choice, rating;
  mpz_t tmp;

  mpz_init(tmp);
  n = M->nrows;
  rating = 0;
  choice = 0;
  for (i=0; i < n; i++)
  {
    int this_row;
    if (mpz_cmp_ui(row_lcm[i], 1) == 0) continue;
    mpz_set_ui(tmp, 1);
    for (j=0; j < n; j++)
      mpz_mul(tmp, tmp, mpq_denref(M->entry[i][j]));
    this_row = n*mpz_sizeinbase(row_lcm[i],2) - mpz_sizeinbase(tmp, 2);
    if (this_row > rating) continue;
    rating = this_row;
    choice = 1+i;
  }
  for (j=0; j < n; j++)
  {
    int this_col;
    if (mpz_cmp_ui(col_lcm[j], 1) == 0) continue;
    mpz_set_ui(tmp, 1);
    for (i=0; i < n; i++)
      mpz_mul(tmp, tmp, mpq_denref(M->entry[i][j]));
    this_col = n*mpz_sizeinbase(col_lcm[j],2) - mpz_sizeinbase(tmp, 2);
    if (this_col > rating) continue;
    rating = this_col;
    choice = -(1+j);
  }
  mpz_clear(tmp);
  return choice;
}



/***************************************************************************/
/* Compute determinant of a (square) Qmat.  Result is placed in det.       */
/* We compute the det by clearing denominators in a "clever" way, and then */
/* we call Zdet on the resulting integer matrix.                           */

void Qdet(mpq_t D, Qmat Min)
{
  Qmat M;
  int n, i, j, k;
  mpz_t tmp, Dnum, Dden;

  n = Min->nrows;
  if (n == 1) { mpq_set(det, Min->entry[0][0]); return; }
  M = Qmat(n, n);
  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      mpq_set(M->entry[i][j], Min->entry[i][j]);
  compute_row_lcm(row_lcm, M);
  compute_col_lcm(col_lcm, M);
  mpz_init(tmp);
  mpz_init(Dnum);
  mpz_init_set_ui(Dden, 1);
  while ((k = choose_row_or_col(M, row_lcm, col_lcm)))
  {
    if (k > 0) /* the choice was row k-1 */
    {
      i = k-1;
      mpz_mul(Dden, Dden, row_lcm[i]);
      for (j=0; j < n; j++)
      {
        mpz_divexact(tmp, row_lcm[i], mpq_denref(M->entry[i][j]));
        mpz_mul(mpq_numref(M->entry[i][j]), mpq_numref(M->entry[i][j]), tmp);
      }
      compute_col_lcm(col_lcm, M);
    }
    else /* the choice was col -k-1 */
    {
      j = -k-1;
      mpz_mul(Dden, Dden, col_lcm[j]);
      for (i=0; i < n; i++)
      {
        mpz_divexact(tmp, col_lcm[j], mpq_denref(M->entry[i][j]));
        mpz_mul(mpq_numref(M->entry[i][j]), mpq_numref(M->entry[i][j]), tmp);
      }
      compute_row_lcm(row_lcm, M);
    }
  }
  mpz_clear(tmp);
  {
    Zmat M_integer;
    M_integer = Zmat_ctor(n, n);
    for (i=0; i < n; i++)
      for (j=0; j < n; j++)
        mpz_set(M_integer[i][j], mpq_numref(M->entry[i][j]));
    Zdet(Dnum, M_integer);
    Zmat_dtor(M_integer);
    mpz_set(mpq_numref(D), Dnum);
    mpz_set(mpq_denref(D), Dden);
    mpq_canonicalize(D);
  }
  Qmat_dtor(M);
}

#endif
