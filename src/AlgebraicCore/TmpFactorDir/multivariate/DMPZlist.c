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

#include <stdlib.h>
#include "DMPZlist.h"
#include "jalloc.h"

/* we give the compiler a hand here with explicit tail recursion */
static int DMPZlist_length1(DMPZlist l, int n)
{
  if (l == NULL) return n;
  return DMPZlist_length1(l->next, n+1);
}

int DMPZlist_length(DMPZlist l)
{
  return DMPZlist_length1(l, 0);
}

/***************************************************************************/

DMPZlist DMPZlist_append(DMPZlist front, DMPZlist back)
{
  DMPZlist tmp;

  if (front == NULL) return back;
  for(tmp = front; tmp->next; tmp = tmp->next) {}
  tmp->next = back;
  return front;
}

/***************************************************************************/

DMPZlist DMPZlist_ctor(DMPZ f, int n)
{
  DMPZlist ans;
  
  ans = (DMPZlist)MALLOC(sizeof(struct DMPZlist_struct));
  ans->poly = f;
  ans->deg = n;
  ans->next = NULL;
  return ans;
}


void DMPZlist_elem_dtor(DMPZlist elem)
{
  DMPZdtor(elem->poly);
  FREE(elem);
}


void DMPZlist_dtor(DMPZlist l)
{
  if (l == NULL) return;
  DMPZlist_dtor(l->next);
  DMPZlist_elem_dtor(l);
}
