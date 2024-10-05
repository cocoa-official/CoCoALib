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

#ifndef JALLOC_H
#define JALLOC_H

#include <stdlib.h>

/*  If DEBUG_MALLOC is not set then we use ordinary malloc.    */
/*  Otherwise we use our own version with some error checking. */

/***************************************************************************/
#ifndef DEBUG_MALLOC

#define MALLOC    malloc
#define FREE      free
#define REALLOC   realloc
#define MEMCHECK

/***************************************************************************/
#else

#define MALLOC(sz)               jmalloc_fn(sz, "malloc", __FILE__, __LINE__)
#define FREE(ptr)                jfree(ptr)
#define REALLOC(ptr, sz)         jrealloc_fn(ptr, sz, "realloc", __FILE__, __LINE__)

extern unsigned int TOTAL_MALLOCKED;
extern unsigned int TOTAL_FREED;
extern unsigned int PEAK_MEM;

#define jmalloc(sz)             jmalloc_fn(sz, "malloc", __FILE__, __LINE__)
#define jmalloc2(sz, str)       jmalloc_fn(sz, str, __FILE__, __LINE__)
#define jcalloc(nels, sz)       jcalloc_fn(nels, sz, "calloc", __FILE__, __LINE__)
#define jcalloc2(nels, sz, str) jcalloc_fn(nels, sz, str, __FILE__, __LINE__)
#define jrealloc(ptr, sz)       jrealloc_fn(ptr, sz, "realloc", __FILE__, __LINE__)
#define jrealloc2(ptr, sz, str) jrealloc_fn(ptr, sz, str, __FILE__, __LINE__)

void *jmalloc_fn(size_t sz, char *str, char *file, int line);
void *jcalloc_fn(int nels, size_t sz, char *str, char *file, int line);
void *jrealloc_fn(void *ptr, int new_sz, char *str, char *file, int line);
void jfree(void *ptr);
void CHECKMARGINS(void *ptr);
#define MEMCHECK printf("TOTAL_MALLOCKED=%d, TOTAL_FREED=%d, PEAK_MEM=%d\n", TOTAL_MALLOCKED, TOTAL_FREED, PEAK_MEM);
#endif

#endif
