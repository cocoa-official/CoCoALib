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

/* Email sent to GMP people hoping that this will be implemented more */
/* sensibly/efficiently in some later release.  28 Oct 1997.          */

#ifndef MPZ_CMP_ABS_H
#define MPZ_CMP_ABS_H

#include "gmp.h"

int mpz_cmp_abs(mpz_t A, mpz_t B);
int mpz_cmp_abs_ui(mpz_t A, long B);

#endif
