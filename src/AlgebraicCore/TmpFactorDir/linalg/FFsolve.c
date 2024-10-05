//   Copyright (c)  1997-2007  John Abbott,  Anna M Bigatti

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

#include "FFsolve.h"
#include "jalloc.h"
#include "jaaerror.h"

/* Two bits used in the return value */
const int FFsolve_new_pivot = 2;
const int FFsolve_invalid = 1;



/* Find a solution to a system of linear equations (in place).          */
/* Several right hand side vectors may be specified.                    */
/* Method is just direct Gaussian elimination -- fine for finite fields */
/* The value returned is the rank of M.                                 */
/* The result is placed in soln; if there is no solution for some rhs   */
/* vector then its solution has first coordinate equal to p, the modulus.*/

int FFsolve(FFmat soln, int *shape, FFmat matrix, FFmat RHS)
{
  int nrows, ncols, r, i, j, k, trows, *pivot_row, result;
  unsigned int count;
  FFelem *swap, invMii, Mji;
  FFelem **M = matrix->entry;
  FFelem **rhs = RHS->entry;
  FFelem p = CurrentFF.prime;
  FFelem limit = CurrentFF.k;
  FFelem shift = CurrentFF.shift;

  if (matrix->nrows != RHS->nrows ||
      soln->nrows != matrix->ncols ||
      soln->ncols != RHS->ncols)
  {
    JERROR(JERROR_MATRIX);
    return 0;
  }
  nrows = matrix->nrows;
  ncols = matrix->ncols;
  r = RHS->ncols;
  trows = 0;
  result = 0;
  count = 0;
  pivot_row = (int*)MALLOC(ncols*sizeof(int));
  for (i=0; i < ncols; i++)
  {
    for (j=trows; j < nrows; j++)
      if ((M[j][i] %= p)) break; /* assignment inside the if!! */
    if (j == nrows)
    {
      if (shape[i]) { FREE(pivot_row); return FFsolve_invalid; }
      continue;
    }
    if (!shape[i])
    {
      shape[i] = 1;
      result |= FFsolve_new_pivot;
      for (k=i+1; k < ncols; ++k) shape[k] = 0;
    }
    pivot_row[trows] = i;
    swap = M[j]; M[j] = M[trows]; M[trows] = swap;
    swap = rhs[j]; rhs[j] = rhs[trows]; rhs[trows] = swap;
    invMii = FFdiv(1, M[trows][i]);
    for (j=i; j < ncols; j++) M[trows][j] = FFmul(M[trows][j]%p, invMii);
    for (j=0; j < r; j++) rhs[trows][j] = FFmul(rhs[trows][j]%p, invMii);
    for (j=0; j < nrows; j++)
    {
      if (j == trows) continue;
      Mji = M[j][i]%p;
      if (Mji == 0) continue;
      Mji = p-Mji;
      for (k=0; k < r; k++) rhs[j][k] += Mji*rhs[trows][k];
      for (k=i+1; k < ncols; k++) M[j][k] += Mji*M[trows][k];
      M[j][i] = 0;
    }
    trows++;
    if (++count < limit) continue;
    count = 0;
    for (j=0; j < nrows; j++)
    {
      for (k=0; k < r; k++)
        if (rhs[j][k] >= shift) rhs[j][k] -= shift;
      for (k=i+1; k < ncols; k++)
        if (M[j][k] >= shift) M[j][k] -= shift;
    }
  }

  for (j=0; j < r; j++)
  {
    /* If any coord beyond trows is non-zero, no solution exists. */
    for (i=trows; i < nrows; i++)
      if ((rhs[i][j] %= p)) break;
    if (i != nrows) { soln->entry[0][j] = p; continue; } /* no soln; mark it & go on */
    k = trows-1;
    for (i=ncols-1; i>=0; i--)
    {
      if (k >= 0 && i == pivot_row[k])
      { soln->entry[i][j] = rhs[k][j]%p; k--; }
      else soln->entry[i][j] = 0;
    }
  }
  FREE(pivot_row);

  return result;
}
