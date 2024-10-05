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

#ifndef DUPZ_DUPFF_H
#define DUPZ_DUPFF_H

#include "DUPZ.h"
#include "DUPFF.h"

DUPZ DUPFF_to_DUPZ(const DUPFF f);
void DUPFF_to_DUPZ2(DUPZ dest, const DUPFF f);
DUPFF DUPZ_to_DUPFF(const DUPZ f);
void DUPZ_to_DUPFF2(DUPFF fbar, const DUPZ f);

#endif
