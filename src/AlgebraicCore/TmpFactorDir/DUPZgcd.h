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

#ifndef DUPZGCD_H
#define DUPZGCD_H

#include "DUPZ.h"

/* WARNING: these functions will divide by zero if the input is */
/* so large that not enough good primes can be found.           */

DUPZ DUPZgcd(const DUPZ fin, const DUPZ gin);
void DUPZgcd3(DUPZ ans, const DUPZ fin, const DUPZ gin);

#endif
