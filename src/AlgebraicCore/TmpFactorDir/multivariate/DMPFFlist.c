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
#include "DMPFFlist.h"
#include "jalloc.h"

/* we give the compiler a hand here with explicit tail recursion */
int DMPFFlist_length1(DMPFFlist l, int n)
{
  if (l == NULL) return n;
  return DMPFFlist_length1(l->next, n+1);
}

int DMPFFlist_length(DMPFFlist l)
{
  return DMPFFlist_length1(l, 0);
}

/***************************************************************************/

DMPFFlist DMPFFlist_append(DMPFFlist front, DMPFFlist back)
{
  DMPFFlist tmp;

  if (front == NULL) return back;
  for(tmp = front; tmp->next; tmp = tmp->next);
  tmp->next = back;
  return front;
}

/***************************************************************************/


DMPFFlist DMPFFlist_ctor(DMPFF f, int n)
{
  DMPFFlist ans;
  
  ans = (DMPFFlist)MALLOC(sizeof(struct DMPFFlist_struct));
  ans->poly = f;
  ans->deg = n;
  ans->next = NULL;
  return ans;
}


void DMPFFlist_elem_dtor(DMPFFlist elem)
{
  DMPFFdtor(elem->poly);
  FREE(elem);
}


void DMPFFlist_dtor(DMPFFlist l)
{
  if (l == NULL) return;
  DMPFFlist_dtor(l->next);
  DMPFFlist_elem_dtor(l);
}
