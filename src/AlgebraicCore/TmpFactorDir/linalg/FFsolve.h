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

#ifndef FFsolve_h
#define FFsolve_h

#include "FFmat.h"

/* Find a solution to a system of linear equations (in place).          */
/* The finite field is specified using FFselect (see file ../FF.h).     */
/* Input: matrix (not necessarily square)                               */
/*        RHS "right hand side" vectors held as columns of a matrix     */
/* Output: the value returned is the rank of matrix.                    */
/*         the solutions are placed in the columns of soln.             */
/*         If no solution exists for some vector then its first coord   */
/*         is set to p (the characteristic).                            */
/* With luck, if the solution space has positive dimension then the     */
/* solutions found should be correspond to the same rational solution   */
/* provided rank(M) mod p is equal to rank(M) without modulus.          */

/* Two bits used in the return value */
extern const int FFsolve_new_pivot;
extern const int FFsolve_invalid;

int FFsolve(FFmat soln, int* shape, FFmat matrix, FFmat RHS);

#endif
