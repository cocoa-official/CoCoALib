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

#include "DMPZevalQ.h"

#include "jalloc.h"
#include <stddef.h>
#include <stdlib.h>


/***************************************************************************/
/* Code for evaluating a distributed multivariate polynomial (DMPZ) at a   */
/* rational point.  The method used is Horner evaluation of a recursive    */
/* representation: this representation can be deduced fairly easily once   */
/* the polynomial is sorted into lexicographic power product order.        */
/* Rational arithmetic is avoided, though this complicates the code :-(    */

/* DMPZevalQ takes a copy of the polynomial, sorts it, creates some        */
/* workspace, then calls DMPZevalQ_iter to traverse and evaluate the poly. */
/* Plenty of scope for refinement, though it seems reasonably fast as is.  */
/* Perhaps put variable D in "stack" instead alloc/free every call?        */

static void DMPZevalQ_iter(mpz_t *stack, DMPZ poly, mpq_t *point, int var, int *degs)
{
  int e, i;
  DMPZ THIS, prev, start;
  mpz_t D;
  
  if (var == NVARS) { mpz_set(stack[var], poly->coeff); return; }
  mpz_init(D);
  mpz_pow_ui(D, mpq_denref(point[var]), degs[var] - poly->exps[var]);
  mpz_set_ui(stack[var], 0);
  THIS = poly;
  while (THIS != NULL)
  {
    start = THIS;
    e = THIS->exps[var];
    prev = NULL; /* just to keep compiler quiet */
    for (; THIS && THIS->exps[var] == e; THIS = THIS->next)
      prev = THIS;
    prev->next = NULL;
    DMPZevalQ_iter(stack, start, point, var+1, degs);
    prev->next = THIS;
    mpz_mul(stack[var+2], stack[var+1], D);
    mpz_add(stack[var], stack[var], stack[var+2]);
    if (THIS != NULL) e -= THIS->exps[var];
    for (i=e; i > 0; i--)
    {
      mpz_mul(stack[var], stack[var], mpq_numref(point[var]));
      mpz_mul(D, D, mpq_denref(point[var])); /* wasteful when THIS == NULL */
    }
  }
  mpz_clear(D);
}


/***************************************************************************/
/* Evaluate a multivariate polynomial with integer coefficients at a       */
/* rational point.  Resulting value is placed in first argument.           */
/* The array "point" must contain at least NVARS elements.                 */

void DMPZevalQ(mpq_t value, DMPZ poly, mpq_t *point)
{
  DMPZ copy;
  mpz_t *stack;
  int i, j, *degs;
  mpz_t D;

  if (poly == NULL) { mpq_set_ui(value, 0, 1); return; }
  copy = DMPZcopy(poly); /* can be clever here if coord = 0 or +1 or -1 */
  copy = DMPZsort(copy, DMPZorder_lex);
  degs = (int*)MALLOC(NVARS*sizeof(int));
  DMPZdegs(degs, copy);
  mpz_init_set_ui(D, 1);
  for (i=0; i < NVARS; i++)
  {
    if (mpz_sgn(mpq_numref(point[i])) == 0) continue;
    for (j=0; j < degs[i]; j++)
      mpz_mul(D, D, mpq_denref(point[i]));
  }
  stack = (mpz_t*)MALLOC((2+NVARS)*sizeof(mpz_t));
  for (i=0; i < 2+NVARS; i++) mpz_init(stack[i]);
  DMPZevalQ_iter(stack, copy, point, 0, degs);
  DMPZdtor(copy);
  mpz_set(mpq_numref(value), stack[0]);
  mpz_set(mpq_denref(value), D);
  mpz_clear(D);
  mpq_canonicalize(value);
  for (i=0; i < 2+NVARS; i++) mpz_clear(stack[i]);
  FREE(stack);
  FREE(degs);
}
