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

#ifndef Zmat_h
#define Zmat_h

#include "Zmat_type.h"
#include "FFmat_type.h"

Zmat Zmat_ctor(int n, int m);
void Zmat_dtor(Zmat M);

void Zmat_swap_rows(Zmat M, int i, int j);
void Zmat_swap_cols(Zmat M, int i, int j);

void Zmat_to_FFmat(FFmat Mp, Zmat M);

void Zmat_mul(Zmat product, Zmat A, Zmat B);
void Zmat_mul_FFmat(Zmat product, Zmat A, FFmat B);

#endif
