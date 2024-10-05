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

#ifndef DUPFFbivariate_lift3_h
#define DUPFFbivariate_lift3_h

// ************************************
// ***** EVERYTHING COMMENTED OUT *****
// ************************************
#if 0

#include "DUPFF.h"
#include "DUPFFlist.h"
//#include "DUPFFfactor_info1.h"

/* The factorization to be lifted is arranged in a "balanced" binary tree */
/* using the struct below.                                                */

struct DUPFFbiv_lift_struct;
typedef struct DUPFFbiv_lift_struct *DUPFFbiv_lifter;

struct DUPFFbiv_lift_struct
{
  DUPFF lcg, lch;  /* leading coeffs of g and h (resp.) or NULL if unknown*/
  DUPFF *g, *h;          /* p-adic expansions of the factors            */
  int r;                 /* Expansion of g,h up to and incl p^(r-1)     */
  DUPFF E;               /* The "error" in component p^r                */
  DUPFF grecip, hrecip;  /* auxiliary values used during lifting        */
  DUPFFbiv_lifter g_lifter, h_lifter;
};


/***************************************************************************/

//DUPFFbiv_lifter DUPFFbiv_lifter_ctor(DUPFF g, DUPFF h, DUPFFbiv_lifter g_lifter, DUPFFbiv_lifter h_lifter, int rmax);
DUPFF** DUPFFbivariate_lift(DUPFF *f, int dy, DUPFFlist factors, unsigned int alpha);
DUPFFbiv_lifter DUPFFbiv_lift_init(DUPFF *f, int dy, DUPFFlist factors, DUPFFlist lcs);
void DUPFFbiv_lift_step(DUPFFbiv_lifter, DUPFF f_cmpt_r);
int DUPFFbiv_lift_output(DUPFF **factors, const DUPFFbiv_lifter THIS, int i);
void DUPFFbiv_lifter_dtor(DUPFFbiv_lifter);

//void DUPFFbiv_lift_revise(DUPFFfactor_info info);
//void DUPFFbiv_lift_revise1(DUPFFfactor_info info, int height);


#endif  // end of commented out section

#endif
