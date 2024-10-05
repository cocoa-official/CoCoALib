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

#include "FFfind_block.h"

/* Function to search for blocks in a matrix.                           */
/* Input: matrix M with nrows rows and ncols columns.                   */
/*        space for the output in the arrays rows and cols.             */
/* Output: in arrays rows[0..nrows-1], cols[0..ncols-1]                 */
/* Those rows/columns belonging to the block found have their entries   */
/* set to 1 all others have entry 0.                                    */
/* Matrix M is left unchanged.                                          */

void FFfind_block(int *rows, int *cols, int nrows, int ncols, uint32 **M)
{
  int i, j, changed;

  for (i=0; i < nrows; i++) rows[i] = 0;
  for (j=0; j < ncols; j++) cols[j] = 0;

  rows[0] = 2;
  while(1)
  {
    changed = 0;
    for (i=0; i < nrows; i++)
    {
      if (rows[i] != 2) continue;
      rows[i] = 1;
      for (j=0; j < ncols; j++)
        if (M[i][j] && cols[j] == 0) cols[j] = changed = 2;
    }
    if (!changed) return;

    changed = 0;
    for (j=0; j < ncols; j++)
    {
      if (cols[j] != 2) continue;
      cols[j] = 1;
      for (i=0; i < nrows; i++)
        if (M[i][j] && rows[i] == 0) rows[i] = changed = 2;
    }
    if (!changed) return;
  }
}
