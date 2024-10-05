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

#ifndef Zmat_det_bound_h
#define Zmat_det_bound_h

#include "Zmat.h"

/************************************************************************/
/* Two functions for computing (Hadamard's) determinant bound of a      */
/* matrix or of one of its principal minors.                            */

/************************************************************************/
/* Compute logarithm of a bound for the determinant of M.               */
/* If M is not square, value returned is 0.                             */
double Zmat_det_bound(Zmat M);

/************************************************************************/
/* This function computes a determinant bound for the leading nxn minor */
/* If n exceeds number of rows or cols in M then 0 is returned.         */
double Zmat_det_bound2(Zmat M, int n);

#endif
