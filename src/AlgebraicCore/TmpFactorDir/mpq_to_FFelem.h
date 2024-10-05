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

#ifndef mpq_to_FFelem_h
#define mpq_to_FFelem_h

#include "gmp.h"
#include "FF.h"

/***************************************************************************/
/* Function to convert a rational num (mpq_t) to a FFelem using the current*/
/* finite field.  Return value is 0 if the conversion failed (i.e. the     */
/* denominator reduced to zero modulo p), and otherwise non-zero.          */
/* If conversion fails then the value of *r is left unaltered.             */

int mpq_to_FFelem(FFelem *r, mpq_t q);

#endif
