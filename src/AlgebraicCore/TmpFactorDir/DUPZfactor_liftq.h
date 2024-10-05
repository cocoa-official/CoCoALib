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

#ifndef DUPZfactor_liftq_h
#define DUPZfactor_liftq_h

#include "DUPZ.h"
#include "DUPFF.h"
#include "DUPFFlist.h"
#include "DUPZfactor_info1.h"

/* The factorization to be lifted is arranged in a binary tree            */
/* using the struct below.                                                */

struct DUPZfactor_lift_struct;
typedef struct DUPZfactor_lift_struct *DUPZfactor_lifter;

struct DUPZfactor_lift_struct
{
  DUPZ g, h;       /* p-adic approximations of the factors              */
  DUPZ corr_factor;/* correction factor for g                           */
  DUPZfactor_lifter g_lifter, h_lifter;
};


/***************************************************************************/

DUPZfactor_lifter DUPZfactor_lifter_ctor(DUPFF g, DUPFF h, DUPZfactor_lifter g_lifter, DUPZfactor_lifter h_lifter);
void DUPZfactor_lift_init(DUPZfactor_info THIS);
void DUPZfactor_lift_step(DUPZfactor_lifter, DUPZ f, mpz_t Q);
int DUPZfactor_lift_output(DUPZ **factors, const DUPZfactor_lifter THIS, int i);
void DUPZfactor_lift_dtor(DUPZfactor_lifter);

void DUPZfactor_lift_revise(DUPZfactor_info info);
/*****void DUPZfactor_lift_revise1(DUPZfactor_info info, int height);*****/


#endif
