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

#ifndef DUPZFACTOR_BOUND_H
#define DUPZFACTOR_BOUND_H

#include "DUPZ.h"

struct DUPZfactor_bound_struct
{
  int df;
  double leading_coeff;
  double trailing_coeff;
  int dmax;
  double Bombieri2;
  double Mahler;
  double root_bound;
  double reverse_root_bound;
  double all_factors; /* a bound for all factors except perhaps one */
  double lift_bound;
};

typedef struct DUPZfactor_bound_struct *bound_info;

bound_info DUPZfactor_bound_ctor(const DUPZ f, int dmax);
void DUPZfactor_bound_dtor(bound_info bd);
void DUPZrefine_bound(bound_info bd, const DUPZ f, int dmax);
bound_info DUPZcopy_refine_bound(bound_info bd, const DUPZ f, int dmax);


#endif
