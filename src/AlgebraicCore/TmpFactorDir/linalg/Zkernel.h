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

#ifndef Zkernel_h
#define Zkernel_h

#include "Zmat.h"

/* Function to compute a Z basis of the kernel of a matrix with integer  */
/* entries.  A randomized algorithm is used so several calls with the    */
/* same arguments may produce different bases; the algorithm was inspired*/
/* by work by G. Havas.  A Z basis is clearly also a Q basis.            */

/* Input:  nrows, ncols -- the size of the matrix M.                     */
/* Output: value returned is dimension of the kernel.                    */
/* Basis vectors are put in Zbasis[0], Zbasis[1], etc.                   */
/* Note that when freeing the space for Zbasis, it was created with space*/
/* for ncols vectors (even if the answer has smaller dimension).         */
/* The matrix M is generally altered by calling this function.           */

int Zmat_left_kernel_Zbasis(Zmat *Zbasis, Zmat M);

#endif
