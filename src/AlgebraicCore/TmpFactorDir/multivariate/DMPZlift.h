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

#ifndef DMPZlift_h
#define DMPZlift_h

#include "DUPZ.h"
#include "DMPZ.h"

struct DMPZlifter_struct
{
  DMPZ f, f2, flcf;
  int x, y;
  int *degs;
  int *substitution;
  
  DUPZ g1;
  DMPZ g, lcg, g_skeleton, g2;
  int *g_degs;

  DUPZ h1;
  DMPZ h, lch, h_skeleton, h2;
  int *h_degs;

  int *lifted;
  int FAILED;
};
typedef struct DMPZlifter_struct* DMPZlifter;

DMPZlifter DMPZlifter_ctor(DMPZ f, DUPZ g1, DUPZ h1, int *substitution, int x);
void DMPZlifter_dtor(DMPZlifter THIS);

void DMPZlift(DMPZlifter THIS);

#endif
