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

#ifndef DMPZ_H
#define DMPZ_H

#include "gmp.h"

/* Distributed multivariate polynomials over the integers. */

extern int NVARS; /* Number of variables (see exps member in struct) */

struct DMPZ_struct;
typedef struct DMPZ_struct *DMPZ;

struct DMPZ_struct
{
  mpz_t coeff;
  int *exps;      /* space for NVARS machine integer exponents */
  DMPZ next;
};

DMPZ DMPZctor(const mpz_t coeff, int *exps);
void DMPZdtor(DMPZ f);
DMPZ DMPZprepend(const mpz_t coeff, int *exps, DMPZ f);
DMPZ DMPZreverse(DMPZ f);
DMPZ DMPZcopy(const DMPZ f);
DMPZ DMPZone(void);
void DMPZcontent(mpz_t content, const DMPZ f);
int DMPZ_nterms(const DMPZ f);
DMPZ DMPZlc(const DMPZ f, int var);
DMPZ DMPZcoeff(const DMPZ f, int var, int deg);
DMPZ* DMPZcoeffs(const DMPZ f, int var);
void DMPZdegs(int *degs, const DMPZ f);
int DMPZdeg(const DMPZ f, int var);
int DMPZtotal_deg(const DMPZ f);
void DMPZmindegs(int *degs, const DMPZ f);
void DMPZnegate(DMPZ f);
DMPZ DMPZdiv_exact(const DMPZ f, const DMPZ g);
void DMPZdiv2z(DMPZ f, const mpz_t n);
void DMPZmul2z(DMPZ f, const mpz_t n);
DMPZ DMPZmul(const DMPZ f, const DMPZ g);
DMPZ DMPZshift(DMPZ f, const int *shift);
DMPZ DMPZadd(const DMPZ a, const DMPZ b);
void DMPZadd3(DMPZ* dest, const DMPZ a, const DMPZ b);

/***************************************************************************/
/* Functions for reordering the terms in a multivariate polynomial (DMPZ). */
/* Recommended usage is:  poly = DMPZ_sort(poly, ordering);                */
/* since DMPZ_sort alters the structure of its first argument.             */
/***************************************************************************/

/***************************************************************************/
/* Define the type of a DMPZ power product order, and define some specific */
/* orders.  Positive return value means first arg is "bigger" than second; */
/* negative return value means first is "smaller" than second; zero return */
/* value means the power products are equal.  The ordering must be total.  */

typedef int (*DMPZexps_cmp)(int*, int*);
int DMPZorder_lex(int* a, int* b);    /* lex order    */
int DMPZorder_deglex(int* a, int* b); /* deglex order */

/***************************************************************************/
/* Sort terms in a polynomial into DECREASING power product order as       */
/* specified by the ordering "cmp".  The structure of the first arg is     */
/* changed by this call!                                                   */
/* Recommended use is:  f = DMPZsort(f, DMPZorder_lex);                    */

DMPZ DMPZsort(DMPZ f, DMPZexps_cmp cmp); /* CHANGES STRUCTURE OF f */


void DMPZprint(const DMPZ f);

#endif
