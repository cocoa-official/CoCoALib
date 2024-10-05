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

#ifndef DUPFFLIST_H
#define DUPFFLIST_H
 
#include "DUPFF.h"
 
struct DUPFFlist_struct;

typedef struct DUPFFlist_struct *DUPFFlist;

struct DUPFFlist_struct
{
  DUPFF poly;
  int deg;
  DUPFFlist next;
};
 
int DUPFFlist_length(DUPFFlist l);

DUPFFlist DUPFFlist_append(DUPFFlist front, DUPFFlist back);
DUPFFlist DUPFFlist_sort(DUPFFlist l);

DUPFFlist DUPFFlist_ctor(DUPFF f, int n);
void DUPFFlist_elem_dtor(DUPFFlist elem);
void DUPFFlist_dtor(DUPFFlist l);


#endif
