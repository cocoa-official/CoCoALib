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

#ifndef DUPI_DUPFF_H
#define DUPI_DUPFF_H

#include "DUPI.h"
#include "DUPFF.h"

/* Functions to convert between DUPIs and DUPFFs.                   */
/* These assume that FFselect has already been called appropriately */

DUPI DUPFF_to_DUPI(const DUPFF f);
void DUPFF_to_DUPI2(DUPI f, const DUPFF g);
DUPFF DUPI_to_DUPFF(const DUPI f);
void DUPI_to_DUPFF2(DUPFF f, const DUPI g);

#endif
