//  Copyright 2013 (c)  John Abbott,  Anna M Bigatti

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
#include "jalloc.h"
#include "PPstream.h"

static int divides(int *pp1, int *pp2, int nvars)
{
  int i;

  for (i=0; i < nvars; i++) if (pp1[i] > pp2[i]) return 0;
  return 1;
}



PPstream PPstream_ctor(int nvars, int (*cmp)(const void*, const void*))
{
  PPstream ans;
  int i;

  ans = (PPstream)MALLOC(sizeof(struct PPstream_struct));
  ans->nvars = nvars;
  ans->not_yet_started = 1;
  ans->avoid_last = 0;
  ans->cmp = cmp;
  ans->frontier_nelems = 1;
  ans->frontier_max_sz = 8;
  ans->frontier = (int**)MALLOC(ans->frontier_max_sz*sizeof(int*));
  ans->frontier[0] = (int*)MALLOC(nvars*sizeof(int));
  for (i=0; i < nvars; i++) ans->frontier[0][i] = 0;

  ans->avoid_nelems = 0;
  ans->avoid_max_sz = 8;
  ans->avoid = (int**)MALLOC(ans->avoid_max_sz*sizeof(int*));
  return ans;
}

void PPstream_dtor(PPstream self)
{
  int i;

  for (i=0; i < self->frontier_nelems; i++)
    FREE(self->frontier[i]);
  FREE(self->frontier);
  for (i=0; i < self->avoid_nelems; i++)
    FREE(self->avoid[i]);
  FREE(self->avoid);
  FREE(self);
}


static void step(PPstream self)
{
  int i, j;
  int nvars = self->nvars;
  int **frontier = self->frontier;
  int **avoid = self->avoid;
  int buf_sz, out, skip, buf_ptr, frontier_ptr;
  int **buffer;

  buf_sz = 0;
  buffer = (int**)MALLOC((1+nvars)*sizeof(int*));
  buffer[buf_sz] = (int*)MALLOC(nvars*sizeof(int));
  for (i=0; i < nvars; i++)
  {
    for (j=0; j < nvars; j++)
      buffer[buf_sz][j] = frontier[0][j];
    buffer[buf_sz][i]++;
    skip = 0;
    for (j=0; !skip && j < self->avoid_nelems; j++)
      if (divides(avoid[j], buffer[buf_sz], nvars)) skip = 1;
    if (skip) continue;
    /* Do not test divisibility by first elem of frontier!! */
    for (j=1; !skip && j < self->frontier_nelems; j++)
      if (divides(frontier[j], buffer[buf_sz], nvars)) skip = 1;
    if (skip) continue;
    buf_sz++;
    buffer[buf_sz] = (int*)MALLOC(nvars*sizeof(int));
  }
  FREE(buffer[buf_sz]);
  qsort((char*)buffer, buf_sz, sizeof(int*), self->cmp);
  if (self->frontier_nelems + buf_sz >= self->frontier_max_sz)
  {
    self->frontier_max_sz = 2*(self->frontier_nelems + buf_sz);
    self->frontier = (int**)REALLOC(self->frontier, self->frontier_max_sz*sizeof(int*));
    frontier = self->frontier;
  }
  out = self->frontier_nelems + buf_sz - 1;
  frontier_ptr = self->frontier_nelems-1;
  buf_ptr = buf_sz-1;
  while (buf_ptr >= 0)
  {
    if (frontier_ptr < 1 || self->cmp(&buffer[buf_ptr], &frontier[frontier_ptr]) > 0)
    {
      frontier[out] = buffer[buf_ptr];
      buffer[buf_ptr] = NULL;
      out--;
      buf_ptr--;
    }
    else
    {
      frontier[out] = frontier[frontier_ptr];
      frontier[frontier_ptr] = NULL; /* not really necessary but safer */
      frontier_ptr--;
      out--;
    }
  }
  FREE(buffer);
  self->frontier_nelems += buf_sz-1;
  FREE(frontier[0]);
  for (i=0; i < self->frontier_nelems; i++)
    frontier[i] = frontier[i+1];
}


int *PPstream_next(PPstream self)
{
  int *ans;

  if (self->frontier_nelems == 0) return NULL;
  if (self->not_yet_started) goto output;
  if (self->avoid_last)
  {
    int i = self->avoid_nelems;
    self->avoid_last = 0;
    if (i+1 == self->avoid_max_sz)
    {
      self->avoid_max_sz *= 2;
      self->avoid = (int**)REALLOC(self->avoid, self->avoid_max_sz*sizeof(int*));
    }

    self->avoid[i] = self->frontier[0];
    self->avoid_nelems++;
    for (i=0; i < self->frontier_nelems-1; i++)
      self->frontier[i] = self->frontier[i+1];
    self->frontier[self->frontier_nelems-1] = NULL; /* for safety */
    if (--self->frontier_nelems == 0) return NULL;
    goto output;
  }
  step(self);
  if (self->frontier_nelems == 0) return NULL;
  output:
  self->not_yet_started = 0;
  {
    int i;
    ans = (int*)MALLOC(self->nvars*sizeof(int));
    for (i=0; i < self->nvars; i++)
      ans[i] = self->frontier[0][i];
  }

  return ans;
}


void PPstream_avoid_last(PPstream self)
{
  if (self->not_yet_started) return;
  self->avoid_last = 1;
}


