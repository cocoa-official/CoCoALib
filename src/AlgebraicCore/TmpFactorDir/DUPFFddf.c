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


/* This is essentially copied from "A Course in Computational Algebraic
   Number Theory" by H Cohen, Springer GTM 138, page 125 (second corrected
   printing) */

/* need stddef for NULL */
#include <stddef.h>
#include "jalloc.h"
#include "DUPFFddf.h"
#include "DUPFFmod.h"

void DUPFFddf(DUPFF *cmpts, int **deg, const DUPFF f) /* f MUST BE SQUAREFREE */
{
  int df, d, p, n, i, ncmpts, dcmpt;
  DUPFF g, x, h, tmp, tmp2;
  DUPFF cmpt, junk; /* just aliases */

  FFelem *copy;
  int *degree;
  DUPFF component;

  p = CurrentFF.prime;
  df = DUPFFdeg(f);

  /* make space for the answer; n is 1 bigger than max poss number of cmpts */
  n = 1; while (n*(n+1) <= 2*df) n++;
  component = (DUPFF)MALLOC(n*sizeof(struct DUPFFstruct));
  degree = (int*)MALLOC(n*sizeof(int));
  /* initialize the arrays just for extra safety - not strictly necessary */
  for (i=0; i<n; i++)
  {
    (component+i)->coeffs = NULL;
    (component+i)->maxdeg = -1;
    (component+i)->deg = -1;
    degree[i] = 0;
  }

  ncmpts = 0;
  tmp = DUPFFnew(df); cmpt = tmp;
  tmp2 = DUPFFnew(df); junk = tmp2;
  h = DUPFFnew(2*df);
  x = DUPFFnew(1);
  x->deg = 1;
  x->coeffs[0] = 0;
  x->coeffs[1] = 1;
  
  g = DUPFFcopy(f);
  DUPFFcopy2(h, x);
  for(d=1; d <= DUPFFdeg(g)/2; d++)
  {
    DUPFFcopy2(tmp, h);
    DUPFFexpt3mod(h, tmp, p, g);   /* h = h^p mod g */
    DUPFFsub3(tmp, h, x);
    DUPFFcopy2(tmp2, g);
    DUPFFgcd2(tmp, tmp2);          /* cmpt = gcd(x^(p^d)-x, g) */
    dcmpt = DUPFFdeg(cmpt);
    if (dcmpt == 0) continue;
    DUPFFdiv4(g, junk, g, tmp);
    DUPFFrem2(h, g);               /* h = h mod g */
    degree[ncmpts] = d;
    copy = (FFelem*)MALLOC((1+dcmpt)*sizeof(FFelem));
    for (i=dcmpt; i>=0; i--) copy[i] = cmpt->coeffs[i];

    (component+ncmpts)->coeffs = copy;
    (component+ncmpts)->deg = dcmpt;
    (component+ncmpts)->maxdeg = dcmpt;
    ncmpts++;
  }


  dcmpt = DUPFFdeg(g);
  if (dcmpt > 0)
  {
    degree[ncmpts] = dcmpt;
    copy = (FFelem*)MALLOC((1+dcmpt)*sizeof(FFelem));
    for (i=dcmpt; i>=0; i--) copy[i] = g->coeffs[i];

    (component+ncmpts)->coeffs = copy;
    (component+ncmpts)->deg = dcmpt;
    (component+ncmpts)->maxdeg = dcmpt;
    ncmpts++;
  }

  DUPFFfree(tmp);
  DUPFFfree(tmp2);
  DUPFFfree(x);
  DUPFFfree(h);
  DUPFFfree(g);

  degree[ncmpts] = 0;
  *cmpts = component;
  *deg = degree;
  return;
}
