#ifndef CoCoA_long32or64_H
#define CoCoA_long32or64_H

//   Copyright (c)  2024  John Abbott,  Anna M. Bigatti

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


#include <climits>

// Macro CoCoA_32BIT_LONG is set iff ulong is less than 64 bits
// (here we assume either 32 or 64 bit)
#if ((((ULONG_MAX/1024)/1024)/1024)/1024)/1024 < 16383
#define CoCoA_32BIT_LONG
#endif

#endif
