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

#ifndef WGD_h
#define WGD_h

#include "gmp.h"

/* Reconstruction of a rational from its modular image.                   */
/* Input:  res is the modular image of the rational.                      */
/*         mod is the corresponding modulus.                              */
/*         Q is the denominator bound (must have 2*Q < mod).              */
/*         Implicit numerator bound is mod/(2*Q).                         */
/* Output: value returned is non-zero if the result is wrong (meaning     */
/*         that there is no rational with denominator at most Q for that  */
/*         (residue, modulus) pair.                                       */
/*         num is the numerator reconstructed                             */
/*         den is the denominator reconstructed (always positive)         */

int modular_to_rational(mpz_t num, mpz_t den, const mpz_t Q, const mpz_t res, const mpz_t mod);
int modular_to_rational2(mpz_t num, mpz_t den, int log2_den, const mpz_t res, const mpz_t mod);

#endif
