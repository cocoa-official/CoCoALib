//  Copyright (c)  2013  John Abbott,  Anna M Bigatti

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


#ifndef PPstream_h
#define PPstream_h

struct PPstream_struct
{
  int not_yet_started;
  int nvars;
  int *last_pp;
  int (*cmp)(const void*, const void*);
  int frontier_nelems, frontier_max_sz;
  int **frontier;
  int avoid_last;
  int avoid_nelems, avoid_max_sz;
  int **avoid;
};

typedef struct PPstream_struct *PPstream;

/***************************************************************************/
PPstream PPstream_ctor(int nvars, int (*cmp)(const void*, const void*));
void PPstream_dtor(PPstream self);
int *PPstream_next(PPstream self);
void PPstream_avoid_last(PPstream self);

#endif
