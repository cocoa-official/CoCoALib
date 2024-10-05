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

/* debugging allocation functions for GMP */

void *debug_malloc(size_t sz)
{
  void *ans;

  ans = malloc(sz);
  printf("ALLOC %p (%lld bytes)\n", ans, (long long)sz);
  return ans;
}

void debug_free(void *ptr, size_t sz)
{
  printf("FREEING %p (%lld bytes)\n", ptr, (long long)sz);
  free(ptr);
}

void *debug_realloc(void *ptr, size_t old_sz, size_t new_sz)
{
  void *ans;
  ans = realloc(ptr, new_sz);
  if (ans == ptr)  return ans;
  printf("FREEING %p (%lld bytes) realloc\n", ptr, (long long)old_sz);
  printf("ALLOC %p (%lld bytes) realloc\n", ans, (long long)new_sz);
  return ans;
}
