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

#ifndef ERROR_H
#define ERROR_H

void JERROR(int errcode);

extern int LAST_ERROR;   /* 0 means no error */
extern int PRINT_ERRORS; /* zero means don't print, non-zero means print */

#define JERROR_DIV_BY_ZERO          1
#define JERROR_FF_BAD_P             2
#define JERROR_NO_MEM               3
#define JERROR_DEG_TOO_LOW          4
#define JERROR_UNSPECIFIED          5
#define JERROR_ZERO_TO_POWER_ZERO   6
#define JERROR_DIV4_ARGS            7
#define JERROR_DUPZFACTORS          8
#define JERROR_DUPISUBTRACT_PRODUCT 9
#define JERROR_DUPZ_REVERSE        10
#define JERROR_HENSEL              11
#define JERROR_EXPT_MOD            12
#define JERROR_ALIASING            13
#define JERROR_MATRIX              14
#define JERROR_PRIMES_EXHAUSTED    15

#endif
