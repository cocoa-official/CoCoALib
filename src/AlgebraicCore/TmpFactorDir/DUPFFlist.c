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
#include "DUPFFlist.h"
#include "jalloc.h"

/* we give the compiler a hand here with explicit tail recursion */
static int DUPFFlist_length1(DUPFFlist l, int n)
{
  if (l == NULL) return n;
  return DUPFFlist_length1(l->next, n+1);
}

int DUPFFlist_length(DUPFFlist l)
{
  return DUPFFlist_length1(l, 0);
}

/***************************************************************************/

DUPFFlist DUPFFlist_append(DUPFFlist front, DUPFFlist back)
{
  DUPFFlist tmp;

  if (front == NULL) return back;
  for(tmp = front; tmp->next; tmp = tmp->next) {}
  tmp->next = back;
  return front;
}

/***************************************************************************/

static DUPFFlist DUPFFlist_merge(DUPFFlist a, DUPFFlist b)
{
  if (a == NULL) return b;
  if (b == NULL) return a;
  if (DUPFFdeg(a->poly) < DUPFFdeg(b->poly))
  {
    a->next = DUPFFlist_merge(a->next, b);
    return a;
  }
  b->next = DUPFFlist_merge(a, b->next);
  return b;
}

DUPFFlist DUPFFlist_sort(DUPFFlist l)
{
  DUPFFlist front, back, iter;
  int n;

  if (l == NULL || l->next == NULL) return l;
  front = l;
  iter = l;
  for (n = (DUPFFlist_length(l)-1)/2; n > 0; n--) iter = iter->next;
  back = iter->next;
  iter->next = NULL;
  
  return DUPFFlist_merge(DUPFFlist_sort(front), DUPFFlist_sort(back));
}
 
/***************************************************************************/

DUPFFlist DUPFFlist_ctor(DUPFF f, int n)
{
  DUPFFlist ans;
  
  ans = (DUPFFlist)MALLOC(sizeof(struct DUPFFlist_struct));
  ans->poly = f;
  ans->deg = n;
  ans->next = NULL;
  return ans;
}


void DUPFFlist_elem_dtor(DUPFFlist elem)
{
  DUPFFfree(elem->poly);
  FREE(elem);
}


void DUPFFlist_dtor(DUPFFlist l)
{
  if (l == NULL) return;
  DUPFFlist_dtor(l->next);
  DUPFFlist_elem_dtor(l);
}
