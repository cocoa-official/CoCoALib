#ifndef RUM_H
#define RUM_H

//   Copyright (c)  2006  Anna Bigatti,  John Abbott

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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/* #ifdef JDEBUG */
/* #include "jalloc.h" */
/* #endif */

/* #ifdef COCOA_WITH_RUM */
/* #include "rum.h" */
/* #endif */

#define RUM_MAX_ALLOCATION 10000
#define RUM_N_STACKS 801

/* when stack is empty reload with RUM_RELOAD pieces */
#define RUM_RELOAD 1000
/* when stack is empty charge RUM_INIT_LOAD pieces */
#define RUM_INIT_LOAD 1000

#define RUM_MAX_SIZE 10000
#define RUM_STD_SIZE 2000

/* every RUM_STAT_FREQ rum-operations, statistics are automatically displayed */
/* if RUM_STAT_FREQ = -1 then statistics are displayed only manually  */
#define RUM_STAT_FREQ -1

void rum_init_all(void);
/* initialize all the rum-stacks 0 slots */

void rum_init(int rum, int size);
/* initialize the rum-stack rum for size slots */

void * rum_malloc(int size);
void rum_free(int size, void * p);

void *rum_realloc(void *oldptr, size_t old_size, size_t new_size, const char * m);
void *rum_calloc(int len, size_t size, const char * m);

#endif  /*  RUM_H  */
