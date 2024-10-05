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

#ifndef DUPZSQFRD_H
#define DUPZSQFRD_H

#include "DUPZ.h"
#include "DUPZlist.h"

int DUPZsqfr(const DUPZ f);

/* Compute square-free decomposition (after discarding any content).    */
/* W A R N I N G: if leading coeff of f is negative then the product    */
/* W A R N I N G: of the sqfr components will be -f (and not f).        */
DUPZlist DUPZsqfrd(const DUPZ f);

#endif
