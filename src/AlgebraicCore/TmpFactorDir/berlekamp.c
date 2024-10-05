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

/* Essentially copied from H Cohen "A Course in Computational Algebraic
   Number Theory" (2nd, corrected printing) pages 131 and 57 */

#include <stdlib.h>
#include "jalloc.h"
#include "berlekamp.h"
#include "DUPFFmod.h"
#include "FFkernel.h"

static void MakeQminusI(FFelem **Q, const DUPFF f)
{
  DUPFF xpi, xp, x;
  int df;
  FFelem p;
  FFelem *ans;
  int i, j;

  df = DUPFFdeg(f);
  if (df > 32767) { *Q = NULL; return; } /* ludicrously big polynomial */
  p = CurrentFF.prime;
  /* create and clear the output matrix */
  ans = (FFelem*)MALLOC(df*df*sizeof(FFelem));
  for (i=df*df-1; i >= 0; i--) ans[i] = 0;

  xp = DUPFFnew(2*df);
  xpi = DUPFFnew(2*df);
  x = DUPFFnew(1);
  x->deg = 1;
  x->coeffs[1] = 1;
  x->coeffs[0] = 0;
  DUPFFexpt3mod(xp, x, p, f);
  xpi->deg = 0;
  xpi->coeffs[0] = 1;
  for (i=1; i < df; i++)
  {
    DUPFFmul3mod(xpi, xpi, xp, f);
    for (j=DUPFFdeg(xpi); j >= 0; j--)
      ans[j*df+i] = xpi->coeffs[j];
    ans[i*df+i] = FFsub(ans[i*df+i], 1);  /* subtract identity matrix */
  }
  DUPFFfree(x);
  DUPFFfree(xp);
  DUPFFfree(xpi);
  *Q = ans;

  return;
}


#if 0
DUPFFlist Berlekamp(const DUPFF f)
{
  FFelem *Q, *B;
  int df, i, j, k, p, r, s;
  DUPFF T, tmp, tmp2;
  DUPFFlist ans, iter, newfactors;

  df = DUPFFdeg(f);
  MakeQminusI(&Q, f);
  r = FFkernel(&B, df, df, Q);
  FREE(Q);

  p = CurrentFF.prime;
  T = DUPFFnew(df);
  tmp = DUPFFnew(2*df);
  tmp2 = DUPFFnew(2*df);
  k = 1;
  ans = DUPFFlist_ctor(DUPFFcopy(f), 1);
  for (j=1; k < r; j++)  /* deliberately skip j=0 */
  {
    for (i=0; i < df; i++)
      T->coeffs[i] = B[j*df+i];
    i = df-1;
    while (T->coeffs[i] == 0) i--;
    T->deg = i;
    for (s=0; s<p && k < r; s++)
    {
      newfactors = NULL;
      for (iter=ans; iter; iter = iter->next)
      {
	DUPFFcopy2(tmp, T);
	tmp->coeffs[0] = s;
	DUPFFcopy2(tmp2, iter->poly);
	DUPFFgcd2(tmp, tmp2);
	if (DUPFFdeg(tmp) == 0 || DUPFFdeg(tmp) == DUPFFdeg(iter->poly))
	{
	  /*	  putchar('|');*/
	  continue;
	}
	/*	if (i > 0) putchar('*'); else putchar('-');*/
        newfactors = DUPFFlist_append(newfactors, DUPFFlist_ctor(DUPFFcopy(tmp), 1));
        DUPFFdiv4(iter->poly, tmp2, iter->poly, tmp);
        if (++k == r) break;
      }
      ans = DUPFFlist_append(ans, newfactors);
    }
  }

  FREE(B);
  DUPFFfree(T);
  DUPFFfree(tmp);
  DUPFFfree(tmp2);

  return ans;
}
#endif

#if 1
DUPFFlist Berlekamp(const DUPFF f)
{
  FFelem *Q, *B;
  int df, i, k, p, r, s;
  DUPFF T, tmp, tmp2;
  DUPFFlist ans, iter, newfactors;

  df = DUPFFdeg(f);
  MakeQminusI(&Q, f);
  r = FFkernel(&B, df, df, Q);
  FREE(Q);

  p = CurrentFF.prime;
  T = DUPFFnew(df);
  tmp = DUPFFnew(2*df);
  tmp2 = DUPFFnew(2*df);
  k = 1;
  ans = DUPFFlist_ctor(DUPFFcopy(f), 1);
  while (k < r)
  {
    for (i=0; i < df; i++)
      T->coeffs[i] = 0;
    T->coeffs[0] = rand()%p;
    for (s=1; s < r; s++)
    {
      FFelem a = rand()%p;
      for (i=0; i < df; i++)
	T->coeffs[i] = FFadd(T->coeffs[i], FFmul(a, B[s*df+i]));
    }
    for (i = df-1; i > 0 && T->coeffs[i] == 0; i--) {}
    if (i < 1) continue;
    T->deg = i;
    for (s=0; s < p && s < 4 && k < r; s++) /* 4 is heuristic */
    {
      newfactors = NULL;
      for (iter=ans; iter; iter = iter->next)
      {
        if (p == 2) /* special case, extra simple */
        {
          DUPFFcopy2(tmp, T);
          tmp2->coeffs[0] = s;
        }
        else /* general case of p != 2 */
        {
          DUPFFcopy2(tmp2, T);
          tmp2->coeffs[0] = FFadd(tmp2->coeffs[0], s);;
          DUPFFrem2(tmp2, iter->poly);
          DUPFFexpt3mod(tmp, tmp2, (p-1)/2, iter->poly);
          if (DUPFFdeg(tmp) < 1) continue;
          tmp->coeffs[0] = FFadd(1, tmp->coeffs[0]);
        }
        DUPFFcopy2(tmp2, iter->poly);
        DUPFFgcd2(tmp, tmp2);
        if (DUPFFdeg(tmp) == 0 || DUPFFdeg(tmp) == DUPFFdeg(iter->poly))
          continue;
        newfactors = DUPFFlist_append(newfactors, DUPFFlist_ctor(DUPFFcopy(tmp), 1));
        DUPFFdiv4(iter->poly, tmp2, iter->poly, tmp);
        if (++k == r) break;
      }
      ans = DUPFFlist_append(ans, newfactors);
    }
  }

  FREE(B);
  DUPFFfree(T);
  DUPFFfree(tmp);
  DUPFFfree(tmp2);

  return ans;
}

#endif


