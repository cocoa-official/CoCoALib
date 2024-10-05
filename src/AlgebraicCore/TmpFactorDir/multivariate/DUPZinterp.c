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

#include "DUPZinterp.h"

#include "jaaerror.h"
#include "DUPZevaluate.h"

void DUPZinterp(DUPZ f, int deg, int *xi, mpz_t *yi)
{
  mpz_t fx, mx, tmp, rem;
  DUPZ modulus, linear;
  int i;
  
//#ifdef FACTOR_DEBUG
//printf("Interpolating univariate of degree %d\n", deg);
//printf("Interpolation points/values are:\n");
//for(i=0;i<=deg;i++){printf("%4d <--> ",xi[i]);mpz_out_str(stdout,10,yi[i]);printf("\n");}
//#endif
  if (deg > f->maxdeg) JERROR(JERROR_DEG_TOO_LOW);
  linear = DUPZnew(1); mpz_set_ui(linear->coeffs[1], 1); linear->deg = 1;
  modulus = DUPZnew(1+deg); mpz_set_ui(modulus->coeffs[0], 1); modulus->deg = 0;
  mpz_init(fx); mpz_init(mx); mpz_init(tmp); mpz_init(rem);
  f->deg = -1;

  for (i=0; i <= deg; i++)
  {
    DUPZevaluate(fx, f, xi[i]);
    DUPZevaluate(mx, modulus, xi[i]);
    mpz_sub(tmp, yi[i], fx);
    mpz_fdiv_qr(tmp, rem, tmp, mx);
#ifdef FACTOR_DEBUG
    if (mpz_sgn(rem) != 0) printf("non integral\n");
#endif
    mpz_neg(tmp, tmp);
    DUPZshift_sub(f, modulus, 0, tmp);
    mpz_set_si(linear->coeffs[0], -xi[i]);
    DUPZmul3(modulus, modulus, linear);
  }

  DUPZfree(linear);  DUPZfree(modulus);
  mpz_clear(tmp);  mpz_clear(fx);  mpz_clear(mx);  mpz_clear(rem);
}
