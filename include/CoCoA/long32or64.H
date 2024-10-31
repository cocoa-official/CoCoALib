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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/long32or64.H,v 1.3 2024/04/12 16:55:52 abbott Exp $
// $Log: long32or64.H,v $
// Revision 1.3  2024/04/12 16:55:52  abbott
// Summary: Changed name of CPP symbol which signifies that longs are 32 bits (redmine 1661/1804)
//
// Revision 1.2  2024/03/23 13:05:44  abbott
// Summary: Corrected (previously it always gave error.. oops)
//
// Revision 1.1  2024/03/23 11:45:57  abbott
// Summary: First attempt at a macro for detecting 32-bitters (redmine 1661)
//
//

#endif
