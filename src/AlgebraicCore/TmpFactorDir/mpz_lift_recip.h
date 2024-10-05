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

#ifndef mpz_lift_recip_h
#define mpz_lift_recip_h

#include "gmp.h"

/* Lift a reciprocal value using Hensel lifting                          */
/* Assume that inputs satisfy for some integer M:                        */
/*   recip*value = 1 mod M  and                                          */
/*   M | modulus  and                                                    */
/*   modulus | M^2                                                       */
/* Upon return recip will be altered so that recip*value = 1 mod modulus */

void mpz_lift_recip(mpz_t recip, mpz_t value, mpz_t modulus);

#endif
