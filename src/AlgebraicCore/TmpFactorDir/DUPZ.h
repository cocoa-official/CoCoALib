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


/* Standard "read once" trick */
#ifndef DUPZ_H
#define DUPZ_H
/***************************************************************************/
/* File really starts here                                                 */

/* DUPZs are dense univariate polynomials with big integer coefficients.   */


/* we have to include stdio.h before gmp.h otherwise i/o functions vanish */
#include <stdio.h>
#include "gmp.h"

struct DUPZstruct
{
  int maxdeg;
  int deg;
  mpz_t *coeffs;
};

typedef struct DUPZstruct *DUPZ;

DUPZ DUPZnew(int maxdeg);
DUPZ int_to_DUPZ(int i);
void DUPZinit(DUPZ x, int maxdeg);
void DUPZfree(DUPZ x);
void DUPZassign(DUPZ lhs, const DUPZ rhs);
void DUPZswap(DUPZ x, DUPZ y);
DUPZ DUPZcopy(const DUPZ f);
void DUPZcopy2(DUPZ dest, const DUPZ src);
int DUPZequal(const DUPZ x, const DUPZ y);

int DUPZdeg(const DUPZ f);
#define DUPZlc(f) (f->coeffs[f->deg])
void DUPZshift(DUPZ f, int xpower);
void DUPZreverse(DUPZ f);

void DUPZneg1(DUPZ f);
void DUPZadd3(DUPZ sum, const DUPZ x, const DUPZ y);
void DUPZsub3(DUPZ diff, const DUPZ x, const DUPZ y);
DUPZ DUPZadd(const DUPZ x, const DUPZ y);
DUPZ DUPZsub(const DUPZ x, const DUPZ y);
DUPZ DUPZmul(const DUPZ x, const DUPZ y);
void DUPZmul3(DUPZ ans, const DUPZ x, const DUPZ y);
void DUPZsquare(DUPZ f);
void DUPZmul2z(DUPZ f, const mpz_t c);
void DUPZmod2z(DUPZ f, const mpz_t c);

void DUPZshift_add(DUPZ f, const DUPZ g, const int deg, const mpz_t coeff);
void DUPZshift_sub(DUPZ f, const DUPZ g, const int deg, const mpz_t coeff);

/* DUPZdiv2z does "floor" division, and ignores non-zero remainders */
void DUPZdiv2z(DUPZ f, mpz_t n);
/* All the division functions below return "prematurely" if an inexact division is found */
DUPZ DUPZdiv(const DUPZ num, const DUPZ den);
DUPZ DUPZrem(const DUPZ num, const DUPZ den);
void DUPZrem2(DUPZ x, const DUPZ m);
void DUPZdiv4(DUPZ quot, DUPZ rem, const DUPZ num, const DUPZ den);

void DUPZmod2(DUPZ f, DUPZ g, mpz_t Q);
void DUPZmdiv2(DUPZ f, DUPZ g, mpz_t Q);
void DUPZmmul3z(DUPZ product, DUPZ f, mpz_t c, mpz_t Q);
void DUPZlinear_shift(DUPZ f, mpz_t a);

/* To avoid trouble with Macintoshes which cannot printf... */
#ifdef FACTOR_DEBUG
#include "DUPZprint.h"
#endif

#endif

