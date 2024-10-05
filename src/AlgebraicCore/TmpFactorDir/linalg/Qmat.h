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

#ifndef Qmat_h
#define Qmat_h

#include "Qmat_type.h"
#include "FFmat_type.h"

Qmat Qmat_ctor(int n, int m);
void Qmat_dtor(Qmat M);

void Qmat_swap_rows(Qmat M, int i, int j);
void Qmat_swap_cols(Qmat M, int i, int j);

int Qmat_to_FFmat(FFmat Mp, Qmat M); /* returns 1 if OK, o/w 0 */

void Qmat_mul(Qmat product, Qmat A, Qmat B);
void Qmat_mul_vec(mpq_t *product, Qmat M, mpq_t *vec);

#endif
