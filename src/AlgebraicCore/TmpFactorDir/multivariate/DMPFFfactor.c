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

#include "DMPFFfactor.h"
#include "DUPFFfactor.h"
#include "jalloc.h"


/***************************************************************************/
/* The following two functions deal specially with the case that the       */
/* polynomial is really univariate in the "variable" var.  Note that var   */
/* is in fact a power product.                                             */

static DMPFF DUPFF_to_DMPFF(const DUPFF f, const int *var)
{
  DMPFF ans;
  int i, j, *exps;

  ans = NULL;
  for (i=0; i <= DUPFFdeg(f); i++)
  {
    if (f->coeffs[i] == 0) continue;
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (j=0; j < NVARS; j++) exps[j] = i*var[j];
    ans = DMPFFprepend(f->coeffs[i], exps, ans);
  }
  return ans;
}

/* fin MUST be a polynomial in var (which is itself really a power product) */

static DMPFFfactors DMPFFfactor_univariate(const DMPFF fin, const int *var)
{
  DMPFFfactors ans;
  DMPFF term, multivariate_copy;
  DUPFF f;
  DUPFFfactors facs;
  DUPFFlist iter;
  int index, df, e;

  for (index = 0; var[index] == 0; index++);
  df = DMPFFdeg(fin, index) / var[index];
  f = DUPFFnew(df);
  /* DMPFFs are not canonical, so same exponent may appear several times */
  for (term = fin; term; term = term->next)
  {
    e = term->exps[index]/var[index];
    f->coeffs[e] += term->coeff;
    if (f->coeffs[e] >= p) f->coeffs[e] -= p;
  }
  while (df >= 0 && f->coeffs[df] == 0) df--;
  f->deg = df;
  ans = DMPFFfactors_ctor();
  if (df < 1) { DMPFFfactors_content(ans, f->coeffs[0]); return ans; }
  facs = DUPFFfactor(f);
  DMPFFfactors_content(ans, facs->content);
  for (iter = facs->list; iter; iter = iter->next)
  {
    DMPFFfactors_multiplicity(ans, iter->deg);
    multivariate_copy = DUPFF_to_DMPFF(iter->poly, var);
    DMPFFfactors_add(ans, multivariate_copy);
    DMPFFdtor(multivariate_copy);
  }
  DUPFFfree(f);
  DUPFFfactors_dtor(facs);
  return ans;
}

/***************************************************************************/

static int DMPFFfactor_nvars(const DMPFF f)
{
  int nvars, i, *degs;

  nvars = 0;
  degs = (int*)MALLOC(NVARS*sizeof(int));
  DMPFFdegs(degs, f);
  for (i=0; i < NVARS; i++)
    if (degs[i] != 0) nvars++;
  FREE(degs);
  return nvars;
}

DMPFFfactors DMPFFfactor(const DMPFF fin)
{
  int i, *degs;
  DMPFFfactors ans;

  /* Handle the case of univariate input specially */
  if (DMPFFfactor_nvars(fin)==1)
  {
    degs = (int*)MALLOC(NVARS*sizeof(int));
    DMPFFdegs(degs, fin);
    for(i = 0; i < NVARS; i++) if (degs[i] != 0) degs[i] = 1;
    ans = DMPFFfactor_univariate(fin, degs);
    FREE(degs);
    return ans;
  }
  ans = DMPFFfactors_ctor();
  return ans;
}
