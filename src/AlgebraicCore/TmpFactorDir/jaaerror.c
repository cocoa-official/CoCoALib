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

#include <stdio.h>
#include "jaaerror.h"

int LAST_ERROR = 0;    /* 0 means no error */
int PRINT_ERRORS = 1;  /* zero means don't print, non-zero means print */


static const char *err_msg[] =
{
  "no error has occurred",
  "division by zero",
  "argument to FFctor not an odd prime or greater than 65535",
  "memory space exhausted",
  "maximum deg in destination polynomial too low",
  "unspecified error",
  "zero to power zero",
  "DUPFFdiv4, illegal arg sharing",
  "DUPZfactors_add: not enough space",
  "DUPIsubtract_product: illegal arg sharing",
  "DUPZreverse: polynomial has zero constant coefficient",
  "Hensel error: factors not coprime",
  "DUPFFexpt3mod: power not positive",
  "Illegal argument aliasing",
  "Matrix problem (wrong size)",
  "CRA failed: prime list exhausted"
};

static int MAXERRCODE = 15;


void JERROR(int errcode)
{
  LAST_ERROR = errcode;
  if (!PRINT_ERRORS) return;
  if (errcode < 1 || errcode > MAXERRCODE)
    fprintf(stderr, "ERROR: unknown error number, %d\n", errcode);
  else
    fprintf(stderr, "ERROR: %s\n", err_msg[errcode]);
  
}
