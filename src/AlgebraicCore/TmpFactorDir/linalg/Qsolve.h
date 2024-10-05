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

#ifndef Qsolve_h
#define Qsolve_h

#include "Qmat.h"

/* Solve a system of linear equations over the rationals.                  */
/* The right hand side contains one or more vectors (as cols of a matrix). */

/* Input: solns space for the solutions (i.e. M->nrows by rhs->ncols).     */
/*        M is the "left hand side" matrix.                                */
/*        rhs is a matrix whose columns are these vectors.                 */
/* Output: the return value is the rank of M.                              */
/*        The solutions are written into soln; a zero denominator of the   */
/*        first coordinate of a solution means that no solution exists.    */

int Qsolve(Qmat soln, Qmat M, Qmat rhs);


#endif
