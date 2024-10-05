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

#include <stddef.h>
#include "DUPZpadic_expand.h"
#include "DUPIsubtract_product.h"
#include "jalloc.h"
#include "DUPI_DUPFF.h"

/***************************************************************************/

/* this divides f by p, and returns the remainder modulo p */
DUPFF DUPZpadic_shift(DUPZ f, int p)
{
  DUPFF ans;
  int i, df;
  mpz_t junk;

  mpz_init(junk);
  df = DUPZdeg(f);
  ans = DUPFFnew(df);
  for (i=0; i <= df; i++)
    ans->coeffs[i] = mpz_fdiv_qr_ui(f->coeffs[i], junk, f->coeffs[i], p);
  /* fix degree of ans */
  for (i=df; (i >= 0) && (ans->coeffs[i] == 0); i--);
  ans->deg = i;
  /* fix degree of f */
  for (i=df; (i >= 0) && (mpz_sgn(f->coeffs[i]) == 0); i--);
  f->deg = i;
  mpz_clear(junk);
  return ans;
}

/***************************************************************************/
/* This mess is used for computing the p-adic expansion of f while         */
/* simultaneously dividing by its leading coefficient.                     */


DUPZpadic_expansion DUPZpadic_expand_init(const DUPZ f, int p, int kmax)
{
  DUPZpadic_expansion ans;

  ans = (DUPZpadic_expansion)MALLOC(sizeof(struct DUPZpadic_struct));
  ans->p = p;
  ans->kmax = kmax;
  ans->k = 0;
  mpz_init_set(ans->lcf_copy, f->coeffs[DUPZdeg(f)]);
  mpz_init(ans->tmp);
  ans->f_copy = DUPZcopy(f);
  ans->lcf_cmpts = (int*)MALLOC(kmax*sizeof(int));
  ans->f_cmpts = (DUPFF*)MALLOC(kmax*sizeof(DUPFF));
  ans->monic_cmpts = (DUPFF*)MALLOC(kmax*sizeof(DUPFF));
  ans->deg0 = DUPFFnew(0); ans->deg0->deg = 0; ans->deg0->coeffs[0] = 1;
  ans->E = DUPInew(DUPZdeg(f));
  ans->fk = DUPInew(DUPZdeg(f));

  return ans;
}



int DUPZpadic_expand_step(DUPZpadic_expansion z)
{
  int i, k;

  k = z->k;
  if (k >= z->kmax) return k; /* silly safety check */
  mpz_fdiv_qr_ui(z->lcf_copy, z->tmp, z->lcf_copy, z->p);
  z->lcf_cmpts[k] = mpz_get_si(z->tmp);
  z->f_cmpts[k] = DUPZpadic_shift(z->f_copy, z->p);
  DUPFF_to_DUPI2(z->fk, z->f_cmpts[k]);
  DUPIadd3(z->E, z->E, z->fk);
  for (i=0; i < k; i++)
  {
    if (z->lcf_cmpts[k-i] == 0) continue;
    z->deg0->coeffs[0] = z->lcf_cmpts[k-i];
    DUPIsubtract_product(z->E, z->monic_cmpts[i], z->deg0);
  }
  z->monic_cmpts[k] = DUPI_to_DUPFF(z->E);
  DUPFFdiv2ff(z->monic_cmpts[k], z->lcf_cmpts[0]);
  z->deg0->coeffs[0] = z->lcf_cmpts[0];
  DUPIsubtract_product(z->E, z->monic_cmpts[k], z->deg0);
  DUPIdiv2z(z->E, z->p);
  z->k++;

  return z->k;
}


void DUPZpadic_expand_end(DUPZpadic_expansion z)
{
  int i;

  if (z == NULL) return;
  DUPZfree(z->f_copy);
  mpz_clear(z->lcf_copy);
  mpz_clear(z->tmp);
  FREE(z->lcf_cmpts);
  DUPFFfree(z->deg0);
  DUPIfree(z->E);
  DUPIfree(z->fk);
  for (i=0; i < z->k; i++)
  {
    DUPFFfree(z->f_cmpts[i]);
    DUPFFfree(z->monic_cmpts[i]);
  }
  FREE(z->f_cmpts);
  FREE(z->monic_cmpts);
  FREE(z);
}
