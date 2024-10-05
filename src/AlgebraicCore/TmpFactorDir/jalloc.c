//   Copyright (c)  1997-2006,2024  John Abbott,  Anna M Bigatti

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
#include <stdio.h>
#include "jalloc.h"

unsigned int TOTAL_MALLOCKED = 0;
unsigned int TOTAL_FREED = 0;
unsigned int CURRENT_MEM = 0;
unsigned int PEAK_MEM = 0;

/* MARGIN is the number of "ints" allocated immediately before and after  */
/* each block.  These margins are checked for integrity when the space is */
/* freed, or when CHECKMARGINS is explicitly called.  The margins are     */
/* invisible to the caller.  Bigger margins give better safety but waste  */
/* more space.                                                            */
#define MARGIN 8

void *jmalloc_fn(size_t sz, char *str, char *file, int line)
{
  int intsize, i;
  int *ans;

  intsize = 1+sz/sizeof(int);
  ans = (int*)malloc((intsize+2*MARGIN)*sizeof(int));
  ans[0] = sz;
  for (i=1; i<MARGIN; i++)  ans[i] = 1234567890;
  for (i=MARGIN; i < MARGIN+intsize; i++)  ans[i] = -999999999;
  for (i=MARGIN+intsize; i < intsize+2*MARGIN; i++)  ans[i] = 1234567890;
  TOTAL_MALLOCKED += sz;
  CURRENT_MEM += sz;
  if (CURRENT_MEM > PEAK_MEM) PEAK_MEM = CURRENT_MEM;
  printf("%s at line %d of file %s\n", str, line, file);
  /* casts to int below are to keep the compiler quiet */
  printf("ALLOC %p (%lld bytes) (block is %p-%p) \tM=%dk F=%dk D=%dk MAX=%dk\n", (void*)(ans+MARGIN), (long long)sz, (void*)ans, (void*)(ans+(intsize+2*MARGIN)), TOTAL_MALLOCKED/1024, TOTAL_FREED/1024, CURRENT_MEM/1024, PEAK_MEM/1024);
  return &ans[MARGIN];
}

void *jcalloc_fn(int nels, size_t sz, char *str, char *file, int line)
{
  size_t i;
  char *space = (char*)jmalloc_fn(nels*sz, str, file, line);

  for (i=0 ; i<nels*sz ; i++) space[i] = 0;
  
  return space;
}

static void trouble(void* addr, void* block_lo, void* block_hi, int sz)
{
  printf("CHECKMARGINS: trouble at %p for block %p-%p (size=%d)\n",
         addr, block_lo, block_hi, sz);
}

void CHECKMARGINS(void *ptr)
{
  int *block;
  int sz, intsize, i;

  block = (int*)((char*)ptr-MARGIN*sizeof(int));
  sz = block[0];
  if (sz < 0)
  {
    printf("CHECKMARGINS: negative size at %p\n", ptr);
    return;
  }
  if (sz > 1000000)
  {
    printf("CHECKMARGINS: suspiciously large block (sz=%d) at %p\n", sz, ptr);
    return;
  }
  for (i=1; i<MARGIN; i++)
    if (block[i] != 1234567890)
      trouble(&block[i], ptr, ((char*)ptr+sz), sz);
  intsize = 1+sz/sizeof(int);
  for (i=MARGIN+intsize; i < intsize+2*MARGIN; i++)
    if (block[i] != 1234567890)
      trouble(&block[i], ptr, ((char*)ptr+sz), sz);
}

void jfree(void *ptr)
{
  int *block;
  int sz, intsize, i;

  if (ptr == NULL)  { printf("CALLED FREE ON ZERO POINTER\n"); return; }
  block = (int*)((char*)ptr-MARGIN*sizeof(int));
  sz = block[0];
  if (sz < 0)  { printf("FREE: block size negative, %d\n", sz); return; }
  intsize = 1+sz/sizeof(int);
  CHECKMARGINS(ptr);
  for (i=0; i < intsize+2*MARGIN; i++)  block[i] = 1122334455;
  printf("FREEING %p (actually %p-%p)\n", ptr, (void*)block, (void*)(block+(intsize+2*MARGIN)));
  TOTAL_FREED += sz;
  CURRENT_MEM -= sz;
  free(block);
}


void *jrealloc_fn(void *ptr, int new_sz, char *str, char *file, int line)
{
  int old_intsize, new_intsize, i;
  int old_sz;
  char *oldblock = (char*)ptr;
  char *newblock;

  old_sz = *((int*)ptr-MARGIN);
  new_intsize = 1+new_sz/sizeof(int);
  old_intsize = 1+old_sz/sizeof(int);
  if (new_intsize <= old_intsize) return ptr;
  newblock = (char *)jmalloc_fn(new_sz, str, file, line);
  for (i=0 ; i < old_sz ; i++ ) newblock[i] = oldblock[i];
  jfree(ptr);
  
  return newblock;
}
