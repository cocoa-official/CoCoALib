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

#include "DMPZeval.h"
#include <stddef.h>
#include "jalloc.h"

void DMPZeval_power_product(mpz_t value, int *pp, const int *point)
{
  int i;
  mpz_t power;
  mpz_init(power);
  
  mpz_set_ui(value, 1);
  for (i=0; i < NVARS; i++)
  {
    if (pp[i] == 0) continue;
    mpz_ui_pow_ui(power, abs(point[i]), pp[i]);
    if (point[i] < 0 && (pp[i]&1)) mpz_neg(power, power);
    mpz_mul(value, value, power);
  }
  mpz_clear(power);
}


static void DMPZeval_term(mpz_t value, const DMPZ term, const int *point)
{
  int i, j;
  
  mpz_set(value, term->coeff);
  for (i=0; i < NVARS; i++)
    for (j=0; j < term->exps[i]; j++)
      if (point[i] >= 0) mpz_mul_ui(value, value, point[i]);
      else { mpz_mul_ui(value, value, -point[i]); mpz_neg(value, value); }
}


void DMPZeval(mpz_t value, const DMPZ poly, const int *point)
{
  DMPZ iter;
  mpz_t term_value;

  mpz_init(term_value);
  mpz_set_ui(value, 0);
  for (iter = poly; iter; iter = iter->next)
  {
    DMPZeval_term(term_value, iter, point);
    mpz_add(value, value, term_value);
  }
  mpz_clear(term_value);
}


/***************************************************************************/

static void DMPZeval_term_partial(mpz_t value, const DMPZ term, const int *point, const int* flag)
{
  int i, j;
  
  mpz_set(value, term->coeff);
  for (i=0; i < NVARS; i++)
  {
    if (flag[i] == 0) continue;
    for (j=0; j < term->exps[i]; j++)
    {
      if (point[i] >= 0) mpz_mul_ui(value, value, point[i]);
      else { mpz_mul_ui(value, value, -point[i]); mpz_neg(value, value); }
    }
  }
}


DMPZ DMPZeval_partial(const DMPZ poly, const int *point, const int *flag)
{
  DMPZ iter, image;
  mpz_t coeff;
  int i, *exps;

  mpz_init(coeff);
  image = NULL;
  for (iter = poly; iter; iter = iter->next)
  {
    DMPZeval_term_partial(coeff, iter, point, flag);
    if (mpz_sgn(coeff) == 0) continue;
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++)
      if (flag[i]) exps[i] = 0;
      else exps[i] = iter->exps[i];
    image = DMPZprepend(coeff, exps, image);
  }
  mpz_clear(coeff);
  return DMPZsort(image, DMPZorder_lex);
}
