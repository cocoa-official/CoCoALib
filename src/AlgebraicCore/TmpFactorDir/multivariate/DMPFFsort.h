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

#ifndef DMPFFsort_h
#define DMPFFsort_h

/***************************************************************************/
/* Functions for reordering the terms in a multivariate polynomial (DMPFF).*/
/* Recommended usage is:  poly = DMPFFsort(poly, ordering);                */
/* since DMPFF_sort alters the structure of its first argument.            */
/***************************************************************************/

#include "DMPFF.h"

/***************************************************************************/
/* Define the type of a DMPFF power product order, and define some specific*/
/* orders.  Positive return value means first arg is "bigger" than second; */
/* negative return value means first is "smaller" than second; zero return */
/* value means the power products are equivalent in the ordering.          */

typedef int (*DMPFFterm_cmp)(DMPFF, DMPFF);
int DMPFForder_lex(DMPFF a, DMPFF b);          /* lex order.            */
void DMPFForder_main_var(int var);             /* select main variable  */
int DMPFForder_main_var_lex(DMPFF a, DMPFF b); /* wrt main var then lex */

/***************************************************************************/
/* Sort terms in a polynomial into DECREASING power product order as       */
/* specified by the ordering "cmp".  The structure of the first arg is     */
/* changed by this call!                                                   */

DMPFF DMPFFsort(DMPFF f, DMPFFterm_cmp cmp); /* CHANGES STRUCTURE OF f */

#endif
