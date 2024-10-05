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

#ifndef DUPZcra_h
#define DUPZcra_h

#include "gmp.h"
#include "DUPZ.h"
#include "DUPFF.h"

/* This Chinese Remaindering function computes a symmetric modular value. */
/* i.e. the value is in the range -(m1*m2)/2 to +(m1*m2)/2.               */

/* f1 is a DUPZ with coefficients in the range -m1/2 to m1/2,             */
/* m1 is a positive (odd?) number (the modulus for f1),                   */
/* f2 is a DUPFF over the current finite field with m2 elements,          */
/* m2 is the characteristic of the current finite field (odd small prime) */
/* It is necessary that m1 and m2 be coprime.                             */
/* The argument f1 is modified to have coefficients in the range -m1*m2/2 */
/* to +m1*m2/2, so that it equals f1 modulo m1, and equals f2 modulo m2.  */
/* Result is 0 if the coeffs of f1 were unchanged, otherwise none-zero.   */

int DUPZcra(DUPZ f1, const mpz_t m1, const DUPFF f2, FFelem m2);

/* BUGS: the argument m2 appears to be quite superfluous.                 */

#endif
