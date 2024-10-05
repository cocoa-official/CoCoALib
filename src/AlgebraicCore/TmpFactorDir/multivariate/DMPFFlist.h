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

#ifndef DMPFFLIST_H
#define DMPFFLIST_H
 
#include "DMPFF.h"
 
struct DMPFFlist_struct;

typedef struct DMPFFlist_struct *DMPFFlist;

struct DMPFFlist_struct
{
  DMPFF poly;
  int deg;
  DMPFFlist next;
};
 
int DMPFFlist_length(DMPFFlist l);

DMPFFlist DMPFFlist_append(DMPFFlist front, DMPFFlist back);

DMPFFlist DMPFFlist_ctor(DMPFF f, int n);
void DMPFFlist_elem_dtor(DMPFFlist elem);
void DMPFFlist_dtor(DMPFFlist l);


#endif
