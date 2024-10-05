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


/* Standard "read once" trick */
#ifndef DUPI_H
#define DUPI_H
/***************************************************************************/
/* File really starts here                                                 */

/* DUPIs are dense univariate polynomials with machine precision integer   */
/* coefficients.  There are ABSOLUTELY NO SAFEGUARDS against integer       */
/* overflow.  Division functions assume all integer divisions are exact.   */

#include "int32.h"

struct DUPIstruct
{
  int maxdeg;
  int deg;
  int32 *coeffs;
};

typedef struct DUPIstruct *DUPI;

DUPI DUPInew(const int maxdeg);
void DUPIfree(DUPI x);
void DUPIassign(DUPI lhs, const DUPI rhs);
void DUPIswap(DUPI x, DUPI y);
DUPI DUPIcopy(const DUPI f);
void DUPIcopy2(DUPI dest, const DUPI src);
int DUPIequal(const DUPI x, const DUPI y);

int DUPIdeg(const DUPI f);

DUPI DUPIadd(const DUPI x, const DUPI y);
void DUPIadd3(DUPI sum, const DUPI x, const DUPI y);
DUPI DUPIsub(const DUPI x, const DUPI y);
void DUPIsub3(DUPI sum, const DUPI x, const DUPI y);
DUPI DUPImul(const DUPI x, const DUPI y);
void DUPImul3(DUPI ans, const DUPI x, const DUPI y);
void DUPIsquare(DUPI f);
DUPI DUPIexpt(const DUPI base, const int power);
void DUPIexpt3(DUPI ans, const DUPI base, const int power);


void DUPIshift_add(DUPI f, const DUPI g, const uint32 deg, const int32 coeff);

/* WARNING: none of these divide functions checks for inexact division */
/* For polynomial division, it is best that the denominator/modulus be monic */
void DUPIdiv2z(DUPI f, int n);
DUPI DUPIdiv(const DUPI num, const DUPI den);
DUPI DUPIrem(const DUPI num, const DUPI den);
void DUPIrem2(DUPI x, const DUPI m);
void DUPIdiv4(DUPI quot, DUPI rem, const DUPI num, const DUPI den);

/* To avoid difficulties with Macintoshes which cannot printf... */
#ifdef FACTOR_DEBUG
#include "DUPIprint.h"
#endif

/***************************************************************************/
/* Standard "read once" trick (tail end) */
#endif

