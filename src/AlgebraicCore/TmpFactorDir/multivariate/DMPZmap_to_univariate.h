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

#ifndef DMPZmap_to_univariate_h
#define DMPZmap_to_univariate_h

/* Evaluate multivariate polynomial f at the given point except for */
/* the variable var which is retained in the univariate image.      */

#include "DUPZ.h"
#include "DMPZ.h"

DUPZ DMPZmap_to_univariate(const DMPZ f, int var, int *point);

#endif
