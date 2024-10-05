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

#include "DMPFFsort.h"
#include <stddef.h>
#include <stdlib.h>

/***************************************************************************/
/* A power product ordering: lexicographic.                                */
/* Higher powers of a variable are "bigger" than lower powers, if powers   */
/* equal look at next variable...                                          */

int DMPFForder_lex(DMPFF a, DMPFF b)
{
  int i;

  for (i=0; i < DMPFF_NVARS; i++)
    if (a->exps[i] != b->exps[i]) return (a->exps[i] - b->exps[i]);
  return 0;
}


/***************************************************************************/
/* Another ordering: first exponent of main_var, then lex.                 */

static int main_var = 0;

void DMPFForder_main_var(int var)
{
// ***BUG*** !!!CHECK VALUE OF var!!!  if (var < 0 || var >= DMPFF_NVARS) printf("DMPFForder_main_var: bad var index\n");
  main_var = var;
}

int DMPFForder_main_var_lex(DMPFF a, DMPFF b)
{
  int i;

  i = main_var;
  if (a->exps[i] != b->exps[i]) return (a->exps[i] - b->exps[i]);
  for (i=0; i < DMPFF_NVARS; i++)
    if (a->exps[i] != b->exps[i]) return (a->exps[i] - b->exps[i]);
  return 0;
}



/***************************************************************************/
/* An implementation of "merge" sort on a linked list (= terms in the poly)*/
/* Terms are sorted into DECREASING power product order according to cmp.  */

static DMPFF DMPFFmerge(DMPFF a, DMPFF b, DMPFFterm_cmp cmp)
{
  if (a == NULL) return b;
  if (b == NULL) return a;
  if (cmp(a, b) > 0)
  {
    a->next = DMPFFmerge(a->next, b, cmp);
    return a;
  }
  b->next = DMPFFmerge(a, b->next, cmp);
  return b;
}


DMPFF DMPFFsort(DMPFF f, DMPFFterm_cmp cmp)
{
  DMPFF front, back, iter;
  int n;

  if (f == NULL || f->next == NULL) return f;
  front = f;
  iter = f;
  for (n = (DMPFF_nterms(f)-1)/2; n > 0; n--) iter = iter->next;
  back = iter->next;
  iter->next = NULL;
  
  return DMPFFmerge(DMPFFsort(front, cmp), DMPFFsort(back, cmp), cmp);
}
 
