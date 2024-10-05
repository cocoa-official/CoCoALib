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

#ifndef ZDET_H
#define ZDET_H

#include "gmp.h"
#include "Zmat.h"

/* Compute determinant of (square) Zmat M, and put the value in det.        */
/* Failure can occur with matrices having large Hadamard bound (e.g. greater*/
/* 10^20000 on 32-bit machines, limit is much higher on 64-bit machines).   */
void Zdet(mpz_t det, Zmat M);


#endif
