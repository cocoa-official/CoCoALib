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

#include "DMPZmap_to_univariate.h"
#include <stddef.h>



DUPZ DMPZmap_to_univariate(const DMPZ f, int var, int *a)
{
  DMPZ iter;
  DUPZ ans;
  int i, deg;
  mpz_t tmp, tmp2;

//printf("Mapping down to univariate in x[%d]:\n", var);
  deg = DMPZdeg(f, var);
  if (deg < 0) return NULL;
  ans = DUPZnew(deg);
  mpz_init(tmp);
  mpz_init(tmp2);
  for (iter = f; iter; iter=iter->next)
  {
    mpz_set(tmp, iter->coeff);
//printf("Doing leading term of ");DMPZprint(iter);
    for (i=0; i < NVARS; i++)
    {
      if (i == var) continue;
      if (iter->exps[i] == 0) continue;
      if (a[i] == 0) { mpz_set_ui(tmp, 0); break; }
      mpz_ui_pow_ui(tmp2, a[i], iter->exps[i]);
      mpz_mul(tmp, tmp, tmp2);
    }
//printf("coeff is ");mpz_out_str(stdout,10,tmp);printf("\n");    
    mpz_add(ans->coeffs[iter->exps[var]], ans->coeffs[iter->exps[var]], tmp);
  }
  for (i=deg; i >= 0; i--) if (mpz_sgn(ans->coeffs[i])) break;
  ans->deg = i;
  mpz_clear(tmp2);
  mpz_clear(tmp);
  return ans;
}
