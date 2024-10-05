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

#include "FFkernel.h"

#include "jalloc.h"

int FFkernel(FFelem **B, int nrows, int ncols, FFelem *M)
{
  FFelem *Mi, *Mj, *c, *ans, *x;
  int i, j, k, s, *d;
  FFelem q, Mik;
  int KerDim;
  unsigned int count;
  FFelem shift = CurrentFF.shift;
  FFelem limit = CurrentFF.k;
  FFelem p = CurrentFF.prime;

  c = (FFelem*)MALLOC(nrows*sizeof(FFelem));
  d = (int*)MALLOC(ncols*sizeof(int));
  for (i=0; i < nrows; i++) c[i] = 0;

  KerDim = 0;
  count = 0;
  for (k=0; k < ncols; k++)
  {
    for (j=0; j < nrows; j++) M[j*ncols+k] %= p;
    for (j=0; j < nrows; j++)
      if (c[j] == 0 && M[j*ncols+k] != 0) break;
    if (j == nrows) { KerDim++; d[k] = -1; continue; }
    Mj = M+(ncols*j);
    q = p - FFdiv(1, Mj[k]);
    for (i=k+1; i < ncols; i++) Mj[i] = FFmul(Mj[i]%p, q);
    for (i=0; i < nrows; i++)
    {
      if (i == j) continue;
      Mi = M+(i*ncols);
      Mik = Mi[k];
      if (Mik == 0) continue;
      Mi[k] = 0;
      /* for (s=k+1; s < ncols; s++) Mi[s] += Mik*Mj[s]; */
      {
        FFelem *Mis = &Mi[k+1],
               *Mjs = &Mj[k+1],
               *last = &Mi[ncols-1];
        while (Mis <= last) *(Mis++) += Mik * *(Mjs++);
      }
    }
    c[j] = k+1; /* this is only a boolean */
    d[k] = j;
    if (++count < limit) continue;
    count = 0;
    for (i=0; i < nrows; i++)
    {
      Mi = M+(i*ncols);
      for (j=k+1; j < ncols; j++)
        if (Mi[j] >= shift) Mi[j] -= shift;
    }
  }
  
  ans = (FFelem*)MALLOC(KerDim*ncols*sizeof(FFelem));
  s = 0;
  x = ans;
  for (k=0; k < ncols; k++)
  {
    if (d[k] != -1) continue;
    for (i=0; i < ncols; i++)
    {
      if (d[i] == -1) x[i] = (i == k);
      else x[i] = M[d[i]*ncols + k];
    }
    x += ncols;
    s++;
  }
  FREE(c);
  FREE(d);

  *B = ans;
  return KerDim;
}
