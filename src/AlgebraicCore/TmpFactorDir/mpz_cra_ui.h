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

#ifndef mpz_cra_ui_h
#define mpz_cra_ui_h

#include "jaa.h"
#include "gmp.h"

/* These Chinese Remaindering functions compute a symmetric modular value */
/* i.e. the value is in the range -(m1*m2)/2 to +(m1*m2)/2.               */

/**************************************************************************/
/* r1 is a GMP integer in the range -m1/2 to m1/2,                        */
/* m1 is a non-negative (odd?) GMP integer (the modulus for r1),          */
/* r2 is an element of the current finite field of m2 elements, hence     */
/* m2 must be a small odd prime.  m1 and m2 must be coprime.              */
/* The value of r1 is altered so that it is equivalent to r1 modulo m1,   */
/* and equivalent to r2 modulo m2.  The new value lies between -m1*m2/2   */
/* and +m1*m2/2.                                                          */
/* The value returned is 0 if r1 was unchanged,and non-zero otherwise.    */
/* The raw function uses tmp as a workspace, and m1modm2 should be equal  */
/* to m1 modulo m2.  This variant is useful if many Chinese Remaindering  */
/* operations have to be performed with the same values of m1 and m2.     */
/**************************************************************************/

int mpz_cra_ui(mpz_t r1, const mpz_t m1, FFelem r2, FFelem m2);
int mpz_cra_ui_raw(mpz_t r1, const mpz_t m1, FFelem r2, FFelem m2,
		   mpz_t tmp, FFelem m1modm2);

/* BUGS: the argument m2 appears to be quite superfluous.                 */

#endif
