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

#include "FF.h"
#include "DMPZfactor_modp.h"
#include "DUPFFfactor.h"
#include "jalloc.h"


/***************************************************************************/
/* The following two functions deal specially with the case that the       */
/* polynomial is really univariate in the "variable" var.  Note that var   */
/* is in fact a power product.                                             */

static DMPZ DUPFF_to_DMPZ(const DUPFF f, const int *var)
{
  DMPZ ans;
  int i, j, *exps;
  mpz_t coeff;

  ans = NULL;
  mpz_init(coeff);
  for (i=0; i <= DUPFFdeg(f); i++)
  {
    if (f->coeffs[i] == 0) continue;
    mpz_set_ui(coeff, f->coeffs[i]);
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (j=0; j < NVARS; j++) exps[j] = i*var[j];
    ans = DMPZprepend(coeff, exps, ans);
  }
  mpz_clear(coeff);
  return ans;
}

/* fin MUST be a polynomial in var (which is itself really a power product) */

static DMPZfactors DMPZfactor_modp_univariate(const DMPZ fin, const int *var)
{
  DMPZfactors ans;
  DMPZ term, multivariate_copy;
  DUPFF f;
  DUPFFlist facs, iter;
  FFelem lcf;
  mpz_t LCF;
  int index, df, e, i;
  const FFelem p = CurrentFF.prime;

  for (index = 0; var[index] == 0; index++) {}
  df = DMPZdeg(fin, index) / var[index];
  f = DUPFFnew(df);
  /* DMPZs are not canonical, so same exponent may appear several times */
  for (i=0; i <= df; i++) f->coeffs[i] = 0;
  for (term = fin; term; term = term->next)
  {
    e = term->exps[index]/var[index];
    f->coeffs[e] += mpz_fdiv_ui(term->coeff, p);
    if (f->coeffs[e] >= p) f->coeffs[e] -= p;
  }
  while (df >= 0 && f->coeffs[df] == 0) df--;
  f->deg = df;
  lcf = DUPFFlc(f);
  ans = DMPZfactors_ctor();
  /* Content = leading coeff because we will force all factors to be monic */
  mpz_init_set_ui(LCF, lcf);
  DMPZfactors_content(ans, LCF);
  mpz_clear(LCF);
  if (df < 1) return ans; /* Silly trivial case */
  facs = DUPFFfactor(f);
  for (iter = facs; iter; iter = iter->next)
  {
    DMPZfactors_multiplicity(ans, iter->deg);
    DUPFFmake_monic(iter->poly);
    multivariate_copy = DUPFF_to_DMPZ(iter->poly, var);
    DMPZfactors_add(ans, multivariate_copy);
    DMPZdtor(multivariate_copy);
  }
  DUPFFfree(f);
  DUPFFlist_dtor(facs);
  return ans;
}

/***************************************************************************/

static int DMPZfactor_nvars(const DMPZ f)
{
  int nvars, i, *degs;

  nvars = 0;
  degs = (int*)MALLOC(NVARS*sizeof(int));
  DMPZdegs(degs, f);
  for (i=0; i < NVARS; i++)
    if (degs[i] != 0) nvars++;
  FREE(degs);
  return nvars;
}

DMPZfactors DMPZfactor_modp(const DMPZ fin, FFelem p)
{
  int i, *degs;
  DMPZfactors ans;
  FF Fp;

  Fp = FFctor(p);
  FFselect(Fp);
  /* Handle the case of univariate input specially */
  if (DMPZfactor_nvars(fin)==1)
  {
    degs = (int*)MALLOC(NVARS*sizeof(int));
    DMPZdegs(degs, fin);
    for(i = 0; i < NVARS; i++) if (degs[i] != 0) degs[i] = 1;
    ans = DMPZfactor_modp_univariate(fin, degs);
    FREE(degs);
    FFdtor(Fp);
    return ans;
  }
  ans = DMPZfactors_ctor();
  FFdtor(Fp);
  return ans;
}
