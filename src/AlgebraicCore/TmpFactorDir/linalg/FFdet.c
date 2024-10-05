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

#include "FFdet.h"
#include "jaaerror.h"

/************************************************************************/
/* Compute determinant of a square matrix over a (small prime) finite   */
/* Method is just direct Gaussian elimination in place -- fine for      */
/* finite fields.  For speed we reduce mod p only when we have to.      */
/* Result is determinant of the (LHS) matrix.                           */
/* The input matrix is overwritten by this routine; on return its       */
/* entries may not be reduced to the range 0..(p-1).                    */


FFelem FFdet(FFmat matrix)
{
  int i, j, k, n;
  FFelem det, *swap, **M, Mii, Mji;
  unsigned int count;
  FFelem p = CurrentFF.prime;
  FFelem limit = CurrentFF.k;
  FFelem shift = CurrentFF.shift;


  if (matrix->nrows != matrix->ncols) JERROR(JERROR_MATRIX);
  n = matrix->nrows;
  M = matrix->entry;
  if (n < 1) return 0;
  det = 1;
  count = 0;
  for (i=0; i < n; i++)
  {
    for (j=i; j < n; j++)
      if (M[j][i]%p != 0) break;
    if (j == n) return 0;
    swap = M[j]; M[j] = M[i]; M[i] = swap;
    for (k=i; k < n; k++) M[i][k] %= p;
    Mii = M[i][i];
    det = (det*Mii)%p;
    if (i != j) det = p - det;
    Mii = p - FFdiv(1, Mii);
    for (j=i+1; j < n; j++)
    {
      Mji = M[j][i]%p;
      if (Mji == 0) continue;
      Mji = (Mji*Mii)%p;
      /* Now add Mji times row i to row k.  Naively we could do this:  */
      /* for (k=i+1; k < n; k++)                                       */
      /*   M[j][k] += Mji*M[i][k];                                     */
      /* Instead we give the compiler a helping hand with the following*/
      {
        FFelem *Mjk = &M[j][i+1],
               *Mik = &M[i][i+1],
               *end = &M[i][n];
        while (Mik != end)
          *Mjk++ += Mji * *Mik++;
      }
    }
    if (++count < limit) continue;
    count = 0;
    for (j=i+1; j < n; j++)
      for (k=i+1; k < n; k++)
        if (M[j][k] >= shift) M[j][k] -= shift;
  }
  return det;
}
