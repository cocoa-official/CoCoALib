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

#include "Zmat_det_bound.h"
#include "mpz_log.h"
#include "add_logs.h"

/*******************************************************************/
/* This computes the log (approx) of Hadamard's determinant bound  */
/* of the leading n-by-n minor of the matrix M.  Assumes that both */
/* matrix dimensions are at least n -- otherwise probably crashes. */
/* We compute both for rows and columns then take the smaller.     */

static double hadamard(Zmat M, int n)
{
  int i, j;
  double this_row, this_col;
  double row_bd, col_bd;

  row_bd = 0;
  for (i=0; i < n; i++)
  {
    this_row = -99; /* log of "zero" */
    for (j=0; j < n; j++)
      if (mpz_sgn(M->entry[i][j]))
        this_row = add_logs(this_row, 2*mpz_log(M->entry[i][j]));
    row_bd += this_row/2;
  }

  col_bd = 0;
  for (i=0; i < n; i++)
  {
    this_col = -99; /* log of "zero" */
    for (j=0; j < n; j++)
      if (mpz_sgn(M->entry[j][i]))
        this_col = add_logs(this_col, 2*mpz_log(M->entry[j][i]));
    col_bd += this_col/2;
  }

  if (col_bd < row_bd) return col_bd;
  return row_bd;
}


double Zmat_det_bound(Zmat M)
{
  if (M->nrows != M->ncols) return 0;
  return hadamard(M, M->nrows);
}


double Zmat_det_bound2(Zmat M, int n)
{
  if (n > M->nrows || n > M->ncols) return 0;
  return hadamard(M, n);
}
