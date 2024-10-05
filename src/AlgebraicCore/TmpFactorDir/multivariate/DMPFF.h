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

#ifndef DMPFF_H
#define DMPFF_H

#include "FF.h"

/* Distributed multivariate polynomials over a small prime finite field. */

extern int DMPFF_NVARS; /* Number of variables (see exps member in struct) */

struct DMPFF_struct;
typedef struct DMPFF_struct *DMPFF;

struct DMPFF_struct
{
  int coeff;
  int *exps;      /* space for DMPFF_NVARS machine integer exponents */
  DMPFF next;
};

DMPFF DMPFFctor(int coeff, int *exps);
void DMPFFdtor(DMPFF f);
DMPFF DMPFFprepend(int coeff, int *exps, DMPFF f);
DMPFF DMPFFreverse(DMPFF f);
DMPFF DMPFFcopy(const DMPFF f);
int DMPFF_nterms(const DMPFF f);
DMPFF DMPFFcoeff(const DMPFF f, int var, int deg);
void DMPFFdegs(int *degs, const DMPFF f);
int DMPFFdeg(const DMPFF f, int var);
int DMPFFtotal_deg(const DMPFF f);
void DMPFFmindegs(int *degs, const DMPFF f);
DMPFF DMPFFdiv_exact(DMPFF f, const DMPFF g);
void DMPFFdiv2ff(DMPFF f, int n);

void DMPFFprint(const DMPFF f);

#endif
