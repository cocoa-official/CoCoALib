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

#include <stddef.h>
#include <stdlib.h>
#include "DUPZlist.h"
#include "jalloc.h"

/* we give the compiler a hand here with explicit tail recursion */
static int DUPZlist_length1(DUPZlist l, int n)
{
  if (l == NULL) return n;
  return DUPZlist_length1(l->next, n+1);
}

int DUPZlist_length(DUPZlist l)
{
  return DUPZlist_length1(l, 0);
}

/***************************************************************************/

DUPZlist DUPZlist_append(DUPZlist front, DUPZlist back)
{
  DUPZlist tmp;

  if (front == NULL) return back;
  for(tmp = front; tmp->next; tmp = tmp->next) {}
  tmp->next = back;
  return front;
}

/***************************************************************************/

DUPZlist DUPZlist_ctor(DUPZ f, int n)
{
  DUPZlist ans;
  
  ans = (DUPZlist)MALLOC(sizeof(struct DUPZlist_struct));
  ans->poly = f;
  ans->deg = n;
  ans->next = NULL;
  return ans;
}


void DUPZlist_elem_dtor(DUPZlist elem)
{
  DUPZfree(elem->poly);
  FREE(elem);
}


void DUPZlist_dtor(DUPZlist l)
{
  if (l == NULL) return;
  DUPZlist_dtor(l->next);
  DUPZlist_elem_dtor(l);
}
