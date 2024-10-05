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

#include "DUPZfactor_refine_fds.h"

/***************************************************************************/
/* This function updates a "factor degree set", an array of booleans.      */
/* If slot n is true then all modular factorizations so far permit the     */
/* polynomial to have a factor of degree n.                                */
/* If fact we use two bits in each slot during updating, and then tidy up  */
/* the array before returning.                                             */
/* The result is non-zero iff we have proved irreducibility.               */

int DUPZfactor_refine_fds(int *fds, DUPFFlist factors)
{
  int i, max, d;
  const int bit=2, both=3;
  DUPFFlist iter;
  
  fds[0] |= bit;
  max = 0;
  for (iter=factors; iter; iter = iter->next)
  {
    d = DUPFFdeg(iter->poly);
    for (i=max; i >= 0; i--) if (fds[i] & bit) fds[i+d] |= bit;
    max += d;
  }
  d = 0;
  for (i=0; i <= max; i++)
    d += (fds[i] = (fds[i] == both));
/*
  printf("Poss degs: ");
  for (i=0; i <= max; i++) if (fds[i]) printf("%d ",i);printf("\n");
*/
  return (d == 2);
}

