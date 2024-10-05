//  Copyright (c)  2015  John Abbott,  Anna M Bigatti

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


/***************************************************************************/
/* Work still to do: projective->affine reduction when appropriate         */
/*                   handle rational points directly (BM_check_final)      */
/*                                                   (reduce_mod_p)        */
/*                   remove edge_pp condition???                           */
/*                   act sensibly when BM_check_final needs more accuracy  */
/*                   last useful pp in projective modular case             */
/*                   dynamic use of degree compatibility as far as poss    */
/*                                                                         */
/***************************************************************************/

/* This include is just for development... */
#include "jalloc.h"

#include <stddef.h>
#include <gmp.h>
#include "mpz_alias.h"
#include "mpz_cra_ui.h"
#include "mpz_log.h"
#include "WGD.h"
#include "BM.h"
#include "PPstream.h"
#include "primes.h"
#include "FF.h"
#include "DUPFF.h"
#include "logi.h"

/***************************************************************************/
/* Two global flags indicating whether we're in affine/projective case     */
/* At any one time exactly one of these should be true: see BM_affine      */
/* and BM_projective below.                                                */

static int affine_case = 0;
static int projective_case = 0;


/***************************************************************************/
/* Standard degree of a power product                                      */

static int deg(int nvars, int *pp)
{
  int ans, i;

  ans = 0;
  for (i=0; i < nvars; i++) ans += pp[i];
  return ans;
}


/***************************************************************************/
/* Evaluate a power product at a given point (modulo p).                   */
/* This is a very naive version which does not need log/exp tables but     */
/* will be rather slow.  It would be more sensible to compute the new      */
/* evaluation vector as a multiple of an earlier vector by a single "indet"*/

static FFelem eval_pp_naive(int nvars, int *pp, FFelem *point)
{
  int i;
  FFelem value = 1;
  
  for (i=0; i < nvars; i++)
  {
    if (pp[i] == 0) continue; /* here we regard 0^0 as being 1 */
    if (point[i] == 0) return 0;
    value = FFmul(value, FFpower(point[i], pp[i]));
  }
  return value;
}


/***************************************************************************/
/* Evaluate a power product at a given point (modulo p).                   */
/* Since we have cheap finite field log/exp available, this is easy.       */
/* A faster version of this function is in multivariate/DMPZeval.c         */

static FFelem eval_pp(int nvars, int *pp, FFelem *point)
{
  const int p = CurrentFF.prime;
  int i;
  FFelem sum = 0;

  if (CurrentFF.LogTable == NULL) return eval_pp_naive(nvars, pp, point);
  for (i=0; i < nvars; i++)
  {
    if (pp[i] == 0) continue; /* here we regard 0^0 as being 1 */
    if (point[i] == 0) return 0;
    sum += FFlog(point[i])*pp[i];
  }
  sum %= (p-1);
  return FFexp(sum);
}


/***************************************************************************/
/* Compute a new row of the BM matrix: evaluate many power products at the */
/* same point, and put these values in the row.                            */

static void eval_row(FFelem *row, int nvars, int *pp, int npoints, FFelem **points)
{
  int i;
  
  for (i=0; i < npoints; i++)
    row[i] = eval_pp(nvars, pp, points[i]);
}


/***************************************************************************/
/* Test whether one power product divides another.                         */
/* Result is 1 if first arg divides the second, otherwise 0.               */

static int divides(int nvars, int *pp1, int *pp2)
{
  int i;

  for (i=0; i < nvars; i++) if (pp1[i] > pp2[i]) return 0;
  return 1;
}


/***************************************************************************/
/* Constructor and destructor for the structure upon which BM_mod_p acts.  */
/* Initially we do not know how big to make the arrays, so we guess, and   */
/* offer the function BM_grow if our guess was too small.                  */
/* We only ever create one of these structures, and the fields pp, issep,  */
/* and GBsize keep their values between successive calls to BM_mod_p.      */
/* Furthermore, the array "det" is used only in the projective case (well, */
/* the affine case uses just det[0]); similarly projective_sep is use only */
/* in the projective case (set by the function fill_projective_sep, q.v.). */
/* npp_max tells us how big the arrays GB, issep, M, and pp are, while npp */
/* tells us how many of the elements are actually used: indices 0..(npp-1).*/

static BM BM_ctor(int npoints, int nvars, int (*cmp)(const void*,const void*))
{
  int i;
  int npp_max = 2*npoints;
  BM self;

  self = (BM)MALLOC(sizeof(struct BuchbergerMoeller_struct));
  self->npoints = npoints;
  self->nvars = nvars;
  self->npp_max = npp_max;
  self->npp = 0;
  self->det = (int*)MALLOC((1+npoints)*sizeof(int));
  self->GB = (int*)MALLOC(npp_max*sizeof(int));
  self->sep = (int*)MALLOC(npoints*sizeof(int));
  self->projective_sep = NULL;
  self->issep = (int*)MALLOC(npp_max*sizeof(int));
  self->pp = (int**)MALLOC(npp_max*sizeof(int*));
  self->M = (FFelem**)MALLOC(npp_max*sizeof(FFelem*));
  self->cmp = cmp;
  for (i=0; i < npp_max; i++)
  {
    self->issep[i] = npoints;
    self->pp[i] = NULL; /* as PPstream creates space for each PP; but note */
                        /* that bm_mod_p2/append_new_pps also call MALLOC. */
    self->M[i] = (FFelem*)MALLOC(2*npoints*sizeof(FFelem));
  }
  self->GBsize = nvars;
  return self;
}


void BM_dtor(BM self)
{
  int i;
  int npp_max = self->npp_max;

  /* The following line may call FREE on NULL as pp[i] could be NULL. */
  for (i=0; i < npp_max; i++) { FREE(self->pp[i]); FREE(self->M[i]); }
  FREE(self->M);
  FREE(self->pp);
  FREE(self->issep);
  if (self->projective_sep) FREE(self->projective_sep);
  FREE(self->sep);
  FREE(self->GB);
  FREE(self->det);
  FREE(self);
}

/* For simplicity we merely double the maximum number of power products  */
/* which the structure can handle.  Repeated calls to BM_grow repeatedly */
/* double the space available.                                           */

static void BM_grow(BM self)
{
  int i;
  int npp_max = self->npp_max;
  
  self->npp_max = 2*npp_max;
  self->GB = (int*)REALLOC(self->GB, 2*npp_max*sizeof(int));
  self->issep = (int*)REALLOC(self->issep, 2*npp_max*sizeof(int));
  self->pp = (int**)REALLOC(self->pp, 2*npp_max*sizeof(int*));
  self->M = (FFelem**)REALLOC(self->M, 2*npp_max*sizeof(FFelem*));
  for (i=npp_max; i < 2*npp_max; i++)
  {
    self->issep[i] = self->npoints;
    self->pp[i] = NULL; /* See comment in BM_ctor. */
    self->M[i] = (FFelem*)MALLOC(2*self->npoints*sizeof(FFelem));
  }
}

/***************************************************************************/
/* Check for linear dependence of the ith row in BM->M on the earlier rows */
/* This function typically consumes almost all the execution time (in the  */
/* call to DUPFFshift_add_raw).                                            */
/* Note that the matrix BM->M has 2*npoints columns: the first npoints cols*/
/* contain the row-reduced values of the power product, the last npoints   */
/* cols contain the linear combination of separator PPs which was used to  */
/* reduce that row.  As side effect, this procedure accumulates a sort of  */
/* "determinant" in the third arg -- this will be the determinant of the   */
/* original separator rows if a full set of separators is found.           */
/* The value returned is -1 if the row reduced to zero, otherwise it is the*/
/* index of the first column with a non-zero entry (from 0 up to npoints-1)*/
/* and the row has been scaled so that this entry is 1.                    */

static int lin_indep(int i, BM self, int *det)
{
  int pivot, j;
  FFelem *row, c;
  FFelem *v = self->M[i];
  int npoints = self->npoints;
  int *sep = self->sep;

  for (j=0; j < npoints; j++) v[j+npoints] = 0;
  pivot = -1;
  for (j=0; j < npoints; j++)
  {
    if (v[j] == 0) continue;
    if (pivot < 0 && sep[j] < 0)
    {
      sep[j] = i;
      v[j+npoints] = 1;
      pivot = j;
      continue;
    }
    if (sep[j] < 0) continue;
    row = self->M[sep[j]];
    c = CurrentFF.prime - FFdiv(v[j], row[j]);
    DUPFFshift_add_raw(&v[j], &row[j], &row[2*npoints-1], c);
  }
  if (pivot < 0) return npoints;
  /* Got a new "separator", rescale it, and reduce previous separators by it */
  *det = FFmul(*det, v[pivot]);
  c = v[pivot];
  for (j=0; j < 2*npoints; j++) v[j] = FFdiv(v[j], c);
  for (j=0; j < npoints; j++)
  {
    if (sep[j] < 0) continue;
    if (j == pivot) continue; /* don't reduce self! */
    row = self->M[sep[j]];
    if (row[pivot] == 0) continue;
    c = CurrentFF.prime - FFdiv(row[pivot], v[pivot]);
    DUPFFshift_add_raw(&row[0], &v[0], &v[2*npoints-1], c);
  }

  return pivot;
}


/***************************************************************************/
/* This function determines if a PP ordering is degree compatible up to a  */
/* given degree: i.e. if all power products up to and including degree     */
/* deg+1 are written down in increasing order then their degrees are       */
/* non-decreasing.                                                         */

static int is_deg_compatible1(int (*cmp)(const void*,const void*), int nvars, int deg)
{
  int i, xmin, xmax, ans;
  int *PP1, *PP2;

  if (nvars == 1) return 1;
  PP1 = (int*)MALLOC(nvars*sizeof(int));
  for (i=0; i < nvars; i++) PP1[i] = 0;
  PP2 = (int*)MALLOC(nvars*sizeof(int));
  for (i=0; i < nvars; i++) PP2[i] = 0;
  /* Firstly find the smallest and largest variables... xmin & xmax */
  xmin = 0;
  xmax = 0;
  for (i=1; i < nvars; i++)
  {
    PP1[i] = 1;
    PP2[xmin] = 1;
    if (cmp(&PP1, &PP2) < 0) xmin = i;
    PP2[xmin] = 0;
    PP2[xmax] = 1;
    if (cmp(&PP1, &PP2) > 0) xmax = i;
    PP2[xmax] = 0;
    PP1[i] = 0;
  }
  /* Now conduct the test (quite simple really!) */
  PP1[xmax] = deg;
  PP2[xmin] = deg+1;
  ans = (cmp(&PP1, &PP2) < 0);
  FREE(PP1);
  FREE(PP2);
  return ans;
}


/***************************************************************************/
/* For the degree compatible case, generate those PPs of the next degree   */
/* which must be considered.  We generate PPs in a way which avoids        */
/* duplicates, but we must still check for divisibility by leading PPs of  */
/* GBasis elements, and must also sort the resulting list of PPs according */
/* to the chosen ordering.                                                 */

static int append_new_pps(BM self, int lo, int hi)
{
  int i, j, k;
  int newlo = hi;
  int newhi = hi;
  int **pp = self->pp;
  int nvars = self->nvars;
  int npoints = self->npoints;

  for (j=0; j < nvars; j++)
  {
    for (i=lo; i < hi; i++)
    {
      if (self->issep[i] == npoints) continue; /* if GB elem, skip */
      if (newhi >= self->npp_max) { BM_grow(self); pp = self->pp; }
      for (k=0; k < j; k++) if (pp[i][k]) break;
      if (k != j) continue;
      if (pp[newhi] == NULL) pp[newhi] = (int*)MALLOC(nvars*sizeof(int));
      for (k=0; k < nvars; k++) pp[newhi][k] = pp[i][k];
      pp[newhi][j]++;
      for (k=0; k < hi; k++)
        if (self->issep[k] == npoints && divides(nvars, pp[k], pp[newhi])) break;
      if (k == hi) newhi++;
    }
  }
  if (newhi == newlo) return newhi;
  qsort((char*)&pp[newlo], newhi - newlo, sizeof(int*), self->cmp);
  return newhi;
}


/***************************************************************************/
/* The normal Buchberger-Moeller algorithm (modulo p).                     */
/* We call this function only if we do not have a previous (plausible)     */
/* reduction path to follow -- this is decided by BM_mod_p.                */
/* This function will work with any power product ordering, though may be  */
/* a lot slower than BM_mod_p2 if there are many variables and the ordering*/
/* is degree compatible.                                                   */
/* In principle, this function could use the partial reduction path already*/
/* followed by BM_mod_p, but the extra coding complication was deemed to be*/
/* prohibitive, and the likely performance loss negligible.                */

static int BM_mod_p1(BM self, FFelem **points)
{
  int i, new_sep;
  int nvars = self->nvars;
  int npoints = self->npoints;
  int *pp;
  PPstream PPs;
  
  /* Clear all separators */
  for (i=0; i < npoints; i++) self->sep[i] = -1;
  self->det[0] = 1;
  self->maxdeg = 0;
  PPs = PPstream_ctor(nvars, self->cmp);
  self->GBsize = 0;
  i = 0;
  while ((pp = PPstream_next(PPs)))
  {
    if (i == self->npp_max) { BM_grow(self); }
    self->pp[i] = pp;
    eval_row(self->M[i], nvars, self->pp[i], npoints, points);
    new_sep = lin_indep(i, self, &self->det[0]);
    self->issep[i] = new_sep;
    if (new_sep == npoints)
    {
      PPstream_avoid_last(PPs);
      self->GB[self->GBsize] = i;
      self->GBsize++;
    }
    i++;
  }
  PPstream_dtor(PPs);
  self->npp = i;
  /* Check we have separated all points. */
  for (i=0; i < npoints; i++)
    if (self->sep[i] == -1) break;
  if (i == npoints) return 1;
  self->issep[0] = npoints; /* erase bad path */
  return -1;
}

/***************************************************************************/
/* The degree-compatible Buchberger-Moeller algorithm (modulo p).          */
/* We call this function only if we do not have a previous (plausible)     */
/* reduction path to follow -- this is decided by BM_mod_p.                */
/* In using append_new_pps this function ASSUMES A DEGREE COMPATIBLE       */
/* ORDERING.  If there are many variables and the ordering is degree       */
/* compatible, this function is likely to be much faster than BM_mod_p1.   */
/* In principle, this function could use the partial reduction path        */
/* already followed by BM_mod_p, but the extra coding complication was     */
/* deemed to be prohibitive, and the likely performance loss negligible.   */

static int BM_mod_p2(BM self, FFelem **points)
{
  int hi, lo, i, j, new_sep;
  int **pp = self->pp;
  int nvars = self->nvars;
  int npoints = self->npoints;
  
  /* Clear all separators */
  for (i=0; i < npoints; i++) self->sep[i] = -1;
  self->det[0] = 1;
  self->maxdeg = 0;
  /* Start with the power product 1 */
  if (pp[0] == NULL) pp[0] = (int*)MALLOC(nvars*sizeof(int));
  for (i=0; i < nvars; i++) pp[0][i] = 0;
  lo = 0;
  hi = 1;
  self->GBsize = 0;
  while (hi > lo)
  {
    for (i=lo; i < hi; i++)
    {
      eval_row(self->M[i], nvars, self->pp[i], npoints, points);
      new_sep = lin_indep(i, self, &self->det[0]);
      self->issep[i] = new_sep;
      if (new_sep == npoints)
      {
        self->GB[self->GBsize] = i;
        self->GBsize++;
      }
    }
    j = append_new_pps(self, lo, hi);
    lo = hi; hi = j;
  }
  self->npp = hi;
  /* Check we have separated all points. */
  for (i=0; i < npoints; i++)
    if (self->sep[i] == -1) break;
  if (i == npoints) return 1;
  self->issep[0] = npoints; /* erase bad path */
  return -1;
}


/***************************************************************************/
/* In the projective case we choose that separator of minimal degree for   */
/* each point.  This function goes through all the (partial) separators to */
/* determine the true ones of minimal degree for each point; their indices */
/* are stored in the array projective_sep, or if that array has already    */
/* been filled it checks that these indices agree.                         */
/* Three return values are possible:                                       */
/*      0  means that this result is compatible with the previous ones     */
/*      1  means that this result is correct and previous ones are not     */
/*     -1  means that this result is incorrect and should be discarded     */

static int fill_projective_sep(BM self)
{
  int i, j, k, nonzero_count, flag;
  int npoints = self->npoints;
  int first_time = self->projective_sep == NULL;

  if (first_time) self->projective_sep = (int*)MALLOC(npoints*sizeof(int));
  /* Replace separator table by table of "minimal" separators */
  flag = 0;
  for (i=0; i < npoints; i++)
  {
    for (k=0; k < self->npp; k++)
    {
      if (self->issep[k] != i) continue;
      nonzero_count = 0;
      for (j=0; j < npoints; j++)
        nonzero_count += (self->M[k][j] != 0);
      if (nonzero_count == 1) break;
    }
    if (k == self->npp) exit(13); /* error, no separator found for point #i */
    if (first_time) { self->projective_sep[i] = k; continue; }
    if (self->projective_sep[i] == k) continue;
    if (self->projective_sep[i] > k) return -1; /* self prime is bad */
    /* if (self->projective_sep[i] < k) ... */
    flag = 1;
    self->projective_sep[i] = k;
    first_time = 1;
  }
  return flag;
}


/***************************************************************************/
/* The projective Buchberger-Moeller algorithm (modulo p).                 */
/* We call this function only if we do not have a previous (plausible)     */
/* reduction path to follow -- this is decided by BM_mod_p.                */
/* It may be noticeably faster to reduce to an affine problem than calling */
/* this routine directly.                                                  */
/* See also the function fill_projective_sep (for the projective case).    */

static int edge_pp(int *first, int *second, int *pp, int nvars)
{
  int i;
  
  for (i = 0; i < nvars; i++) if (pp[i]) break;
  if (i == nvars) return 0; /* self should never happen */
  *first = i;
  for (i++; i < nvars; i++) if (pp[i]) break;
  if (i == nvars) { *second = *first; return 1; }
  *second = i;
  for (i++; i < nvars; i++) if (pp[i]) break;
  if (i == nvars) return 1;
  return 0;
}

static int BM_mod_projective_finished(BM self, int lo, int hi)
{
  int i, j, k, l, r, s, d, count, **pp, *corner, deg, finished;
  int npoints = self->npoints;
  int nvars = self->nvars;

  if (nvars == 1) return 1; /* silly case */
  corner = (int*)MALLOC(npoints*sizeof(int));
  pp = (int**)MALLOC(npoints*sizeof(int*));
  for (i=0; i < npoints; i++) corner[i] = -1;
  i = 0;
  for (j=lo; j < hi; j++)
    if (self->issep[j] != npoints) pp[i++] = self->pp[j];
  deg = 0;
  for (i=0; i < nvars; i++) deg += pp[0][i];
  finished = 0; /* set to true only at the very end */
  /* Now check to see if any vertex is connected to the opposite face... */
  for (i=0; i < nvars; i++)                 /* For each vertex of the simplex... */
  {
    for (j=0; j < npoints; j++)
      if (pp[j][i] == deg) break;
    if (j == npoints) continue;             /* skip it if it is in the GBasis    */
    /* Here we could just check quickly for an empty co-degree for more speed.*/
    corner[j] = j;
    for (d=deg-1; d >= 0; d--)              /* Work away from vertex by co-degree*/
    {
      count = 0;
      for (k=0; k < npoints; k++)           /* For each point in the simplex...  */
      {
        if (pp[k][i] != d) continue;        /* of the correct co-degree...       */
                                            /* check that some single step toward*/
                                            /* vertex i goes to a point connected*/
                                            /* to vertex i.                      */

        for (l=0; l < nvars; l++)           /* For each possible "direction" (or */
                                            /* neighbour) towards vertex i...    */
        {
          if (l == i || pp[k][l] == 0) continue;
          for (r=0; r < npoints; r++)       /* Look for neighbour outside the GB */
          {
            if (pp[r][i] != d+1 || pp[r][l] != pp[k][l]-1) continue;
            for (s=0; s < nvars; s++)       
            {
              if (s == i || s == l) continue;
              if (pp[r][s] != pp[k][s]) break;
            }
            if (s == nvars) break;          /* OK, we have found the neighbour   */
          }
          if (r == npoints) continue;       /* [the neighbour was in the GB, so  */
                                            /* try next direction]               */
          if (corner[r] != j) continue;     /* that neighbour wasn't connected, try another */
          corner[k] = j;                    /* OK, we are connected.  Break.     */
          count++;
          break;
        }
      }
      if (count == 0) break;                /* We have proof of disconnection.   */
    }
    if (d < 0) goto tidy_up;                /* We reached the other side.        */
  }
  finished = 1;
  tidy_up:
  FREE(corner);
  FREE(pp);
  return finished;
}


static int BM_mod_projective(BM self, FFelem **points)
{
  int hi, lo, i, newhi, new_sep, nseps, deg;
  int **pp = self->pp;
  int nvars = self->nvars;
  int npoints = self->npoints;
  int j, k, maxedges, nedges, **edge;
  
  nedges = 0;
  maxedges = (nvars*(nvars-1))/2;
  edge = (int**)MALLOC(nvars*sizeof(int*));
  for (i=1; i < nvars; i++)
  {
    edge[i] = (int*)MALLOC(i*sizeof(int));
    for (j=0; j < i; j++) edge[i][j] = 0;
  }
      
  deg = 0;
  /* Start with the power product 1 */
  if (pp[0] == NULL) pp[0] = (int*)MALLOC(nvars*sizeof(int));
  for (i=0; i < nvars; i++) pp[0][i] = 0;
  lo = 0;
  hi = 1;
  self->GBsize = 0;
  while (hi > lo)
  {
    /* Clear all (pseudo-)separators */
    for (i=0; i < npoints; i++) self->sep[i] = -1;
    nseps = 0;
    self->det[deg] = 1;
    for (i=lo; i < hi; i++)
    {
      eval_row(self->M[i], nvars, self->pp[i], npoints, points);
      new_sep = lin_indep(i, self, &self->det[deg]);
      self->issep[i] = new_sep;
      if (new_sep != npoints) { nseps++; continue; }
      self->GB[self->GBsize] = i;
      self->GBsize++;
      if (nedges == maxedges) continue;
      if (!edge_pp(&j, &k, self->pp[i], nvars)) continue;
      if (j < k) { if (!edge[k][j]) nedges++; edge[k][j] = 1; continue; }
      for (k=0; k < j; k++) { if (!edge[j][k]) nedges++; edge[j][k] = 1; }
      for (k=j+1; k < nvars; k++) { if (!edge[k][j]) nedges++; edge[k][j] = 1; }
    }
    if (nedges == maxedges && nseps == npoints &&
        BM_mod_projective_finished(self, lo, hi))
      break;
    newhi = append_new_pps(self, lo, hi);
    lo = hi; hi = newhi;
    deg++;
  }
  for (i=1; i < nvars; i++) FREE(edge[i]);
  FREE(edge);
  self->npp = hi;
  self->maxdeg = deg;
  fill_projective_sep(self);
  return 1;
}


/***************************************************************************/
/* Apply the BM algorithm modulo p, using previous reduction path if poss. */
/* This is the main function we use to compute many GBases modulo various  */
/* p.  Normally it will follow the reduction path of all previous calls    */
/* (and return a value of 0) unless it detects a discrepancy.  If this     */
/* prime is "bad" a value of -1 is returned, otherwise BM_mod_p1 is called */
/* to determine an improved reduction path, and a value of 1 is returned.  */
/* Return value is 0 if result is compatible with previous path,           */
/*                -1 if this result is wrong                               */
/*                 1 if this result proves all previous ones wrong         */

static int BM_mod_p(BM self, FFelem **points)
{
  int i, j, new_sep, d;
  int npoints = self->npoints;
  int nvars = self->nvars;

  d = 0;
  self->det[d] = 1; /* the determinant will accumulate in here */
  /* Clear all separators */
  for (i=0; i < npoints; i++) self->sep[i] = -1;
  for (i=0; i < self->npp; i++)
  {
    if (projective_case && i > 0 && deg(nvars, self->pp[i]) > deg(nvars, self->pp[i-1]))
    {
      d = deg(nvars, self->pp[i]);
      self->det[d] = 1;
      for (j=0; j < npoints; j++) self->sep[j] = -1;
    }
    eval_row(self->M[i], nvars, self->pp[i], npoints, points);
    new_sep = lin_indep(i, self, &self->det[d]);
    if (new_sep == self->issep[i]) continue;
    if (new_sep > self->issep[i]) return -1; /* self prime is bad        */
    if (self->projective_sep) FREE(self->projective_sep);
    self->projective_sep = NULL;
    break;
  }
  if (projective_case) self->maxdeg = d; else self->maxdeg = 0;
  if (self->npp > 0 && i == self->npp)
  {
    if (affine_case) return 0; /* compatible with previous */
    return fill_projective_sep(self);             /* can be -1, 0 or 1 */
  }

  /* Current prime has shown that all previous ones were bad. */
  /* Although wasteful, the simplest is to compute everything from the start */

  if (projective_case) return BM_mod_projective(self, points);
  if (is_deg_compatible1(self->cmp, nvars, npoints))
    return BM_mod_p2(self, points);
  return BM_mod_p1(self, points);
}


/***************************************************************************/
/* Two functions for checking that a set of points mod p are distinct.     */
/* One function handles the affine case, the other the projective.         */
/* A return value of 1 means that the points are all distinct, 0 means that*/
/* at least two are equal.                                                 */

static int affine_check_mod_p(int nvars, int npoints, FFelem **points)
{
  int i, j, k;

  for (i=1; i < npoints; i++)
    for (j=0; j < i; j++)
    {
      for (k=0; k < nvars; k++)
        if (points[i][k] != points[j][k]) break;
      if (k == nvars) return 0; /* two points are equal */
    }
  return 1; /* OK, all points are distinct */
}


static int projective_check_mod_p(int nvars, int npoints, FFelem **points)
{
  int i, j, k, c;
  
  for (i=0; i < npoints; i++)
  {
    for (j=0; j < nvars; j++)
      if (points[i][j] != 0) break;
    if (j == nvars) return 0; /* disallow (0, 0, ..., 0) in projective case */
    for (j=0; j < i; j++)
    {
      for (k=0; k < nvars; k++)
        if (points[i][k] != 0 && points[j][k] != 0) break;
      if (k == nvars) continue;
      c = FFdiv(points[i][k], points[j][k]);
      for (k=0; k < nvars; k++)
        if (FFmul(c, points[j][k]) != points[i][k]) break;
      if (k == nvars) return 0; /* two projective points are equal */
    }
  }
  return 1; /* OK, all projective points are distinct */
}


/***************************************************************************/
/* Reduce the input points modulo the current prime.                       */
/* Return value is 0 if conversion was unsuccessful, otherwise non-zero.   */
/* The values of the points modulo p are placed in points_modp.            */

static int reduce_mod_p(FFelem **points_modp, int nvars, int npoints, mpz_t **points)
{
  int i, j;
  int p = CurrentFF.prime;
  
  for (i=0; i < npoints; i++)
    for (j=0; j < nvars; j++)
      points_modp[i][j] = mpz_fdiv_ui(points[i][j], p);

  if (affine_case) return affine_check_mod_p(nvars, npoints, points_modp);
  return projective_check_mod_p(nvars, npoints, points_modp);
}

/***************************************************************************/
/* Constructor, destructor etc for the struct we use for building the      */
/* non-modular result.  We use a struct to avoid using global variables or */
/* passing lots of parameters.  Supposedly memory management should be     */
/* easier this way...                                                      */


static BMGB BMGB_ctor(int npoints, int GBsize)
{
  int i, j;
  BMGB self;

  self = (BMGB)MALLOC(sizeof(struct BMGB_struct));
  
  self->npoints = npoints;
  self->GBsize = GBsize;
  self->compute_GB = 1;
  self->compute_separators = 1;
  self->last_fail = 0;
  self->target_height = 2; /* could use a more intelligent value here */
  mpz_init_set_ui(self->modulus, 1);
  self->det = (mpz_t*)MALLOC((1+npoints)*sizeof(mpz_t));
  self->GB_mod = (mpz_t**)MALLOC(GBsize*sizeof(mpz_t*));
  self->GB = (mpq_t**)MALLOC(GBsize*sizeof(mpq_t*));
  for (i=0; i <= npoints; i++) mpz_init(self->det[i]);
  for (i=0; i < GBsize; i++)
  {
    self->GB_mod[i] = (mpz_t*)MALLOC(npoints*sizeof(mpz_t));
    self->GB[i] = (mpq_t*)MALLOC(npoints*sizeof(mpq_t));
    for (j=0; j < npoints; j++)
    {
      mpz_init(self->GB_mod[i][j]);
      mpq_init(self->GB[i][j]);
    }
  }
  self->sep_mod = (mpz_t**)MALLOC(npoints*sizeof(mpz_t*));
  self->sep = (mpq_t**)MALLOC(npoints*sizeof(mpq_t*));
  for (i=0; i < npoints; i++)
  {
    self->sep_mod[i] = (mpz_t*)MALLOC(npoints*sizeof(mpz_t));
    self->sep[i] = (mpq_t*)MALLOC(npoints*sizeof(mpq_t));
    for (j=0; j < npoints; j++)
    {
      mpz_init(self->sep_mod[i][j]);
      mpq_init(self->sep[i][j]);
    }
  }
  return self;
}

void BMGB_dtor(BMGB self)
{
  int i, j;

  if (self == NULL) return;
  mpz_clear(self->modulus);
  for (i=0; i < self->GBsize; i++)
  {
    for (j=0; j < self->npoints; j++)
    {
      mpz_clear(self->GB_mod[i][j]);
      mpq_clear(self->GB[i][j]);
    }
    FREE(self->GB_mod[i]);
    FREE(self->GB[i]);
  }
  for (i=0; i <= self->npoints; i++) mpz_clear(self->det[i]);

  for (i=0; i < self->npoints; i++)
  {
    for (j=0; j < self->npoints; j++)
    {
      mpz_clear(self->sep_mod[i][j]);
      mpq_clear(self->sep[i][j]);
    }
    FREE(self->sep_mod[i]);
    FREE(self->sep[i]);
  }
  FREE(self->sep_mod);
  FREE(self->sep);
  FREE(self->GB_mod);
  FREE(self->GB);
  FREE(self->det);
  FREE(self);
}

/***************************************************************************/
/* Compute log of value of a power product at a given point.               */
/* log_point contains the logs of the absolute values of the coordinates of*/
/* the point.  This function assumes that the point has integer coords.    */
/* A zero coord is indicated by a negative logarithm, and if the value of  */
/* the power product is zero the result is negative.                       */

static double log_eval(int nvars, int *pp, double *log_point)
{
  double logval;
  int i;

  logval = 0;
  for (i=0; i < nvars; i++)
  {
    if (pp[i] > 0 && log_point[i] < 0) return -1.0e20;
    logval += pp[i] * log_point[i];
  }
  return logval;
}


/***************************************************************************/
/* This function performs a final verification of the separators/GBasis.   */
/* It does this by checking that each GBasis element evaluates to zero at  */
/* every point.  This is already known to be true modulo modulus*p, and in */
/* most cases we need only verify that this implies the values are indeed  */
/* zero without modulus; analogously for the separators.  We use only a    */
/* crude check as this usually suffices.                                   */
/* The value returned is non-zero if the verification succeeded, and       */
/* zero otherwise (in which case char0->target_height is set).             */

static void hack_sep_table(BM modp, int degree)
{
  int i, j, top, bottom;
  int nvars = modp->nvars;
  int npoints = modp->npoints;

  for (top = modp->npp - 1; deg(nvars, modp->pp[top]) > degree; top--) {}
  if (top == 0) return; /* only happens if there is just a single point */
  for (bottom = top-1; deg(nvars, modp->pp[bottom]) == degree; bottom--) {}
  for (i=0; i < npoints; i++)
  {
    for (j=top; j > bottom; j--) if (modp->issep[j] == i) break;
    if (j == bottom) j = -99999999;
    modp->sep[i] = j;
  }
}

static int BM_check_final(BM modp, BMGB char0, mpz_t modulus, int p, mpz_t **points)
{
  int i, j, k, last_deg;
  double tmp, max, bound, log_zero = -1.0e20; /* a large negative number */
  double *pp_value, *log_coeff, *log_point;
  mpz_t D;
  int npoints = modp->npoints;
  int GBsize = modp->GBsize;
  int **pp = modp->pp;
  int *GB = modp->GB;
  int *sep = modp->sep;
  int *reducer = modp->sep;
  int npp = modp->npp;
  int nvars = modp->nvars;

  if (projective_case) sep = modp->projective_sep;
  log_point = (double*)MALLOC(nvars*sizeof(double));
  log_coeff = (double*)MALLOC(npoints*sizeof(double));
  pp_value = (double*)MALLOC(npp*sizeof(double));
  for (i=0; i < npp; i++) pp_value[i] = log_zero;
  for (i=0; i < npoints; i++)
  {
    for (j=0; j < nvars; j++)
      if (mpz_sgn(points[i][j]) == 0) log_point[j] = log_zero;
      else log_point[j] = mpz_log(points[i][j]);
    for (j=0; j < npp; j++)
    {
      tmp = log_eval(nvars, pp[j], log_point);
      if (tmp > pp_value[j]) pp_value[j] = tmp;
    }
  }
  FREE(log_point);

  last_deg = 0;
  bound = 0;
  mpz_init(D);
  if (char0->compute_GB)
  {
    for (i=0; i < GBsize; i++)
    {
      if (projective_case && last_deg != deg(nvars, pp[GB[i]]))
      {
        last_deg = deg(nvars, pp[GB[i]]);
	hack_sep_table(modp, last_deg);
      }
      
      mpz_set_ui(D, 1);
      for (j=0; j < npoints; j++)
        if (mpq_sgn(char0->GB[i][j]) == 0) log_coeff[j] = log_zero;
        else { log_coeff[j] = mpq_log(char0->GB[i][j]); mpz_lcm(D, D, mpq_denref(char0->GB[i][j])); }
      max = pp_value[GB[i]];
      for (k=0; k < npoints; k++)
      {
	if (log_coeff[k] == log_zero) continue;
        tmp = log_coeff[k] + pp_value[reducer[k]];
	if (tmp > max) max = tmp;
      }
      max += mpz_log(D) + logi(1+npoints);
      if (max > bound) bound = max;
    }
  }
  if (char0->compute_separators)
  {
    for (i=0; i < npoints; i++)
    {
      if (projective_case && last_deg != deg(nvars, pp[sep[i]]))
      {
        last_deg = deg(nvars, pp[sep[i]]);
	hack_sep_table(modp, last_deg);
      }
      
      mpz_set_ui(D, 1);
      for (j=0; j < npoints; j++)
        if (mpq_sgn(char0->sep[i][j]) == 0) log_coeff[j] = log_zero;
        else { log_coeff[j] = mpq_log(char0->sep[i][j]); mpz_lcm(D, D, mpq_denref(char0->sep[i][j])); }
      max = log_zero;
      for (k=0; k < npoints; k++)
      {
	if (log_coeff[k] == log_zero) continue;
	tmp = log_coeff[k] + pp_value[reducer[k]];
	if (tmp > max) max = tmp;
      }
      max += mpz_log(D) + logi(npoints);
      if (max > bound) bound = max;
    }
  }
/*printf("--actual precision: %g\n", mpz_log(modulus) + logi(p)); */
/*printf("--predicted bound: %g\n", bound);                       */
  mpz_clear(D);
  FREE(log_coeff);
  FREE(pp_value);
  char0->target_height = (int)(bound/logi(2));
  if (mpz_log(modulus) + logi(p) > bound) return 1;
  return 0;
}


/***************************************************************************/
/* Check to see whether the GBasis elements are "stable".                  */
/* Return value is non-zero if the result in char0 is stable, otherwise 0. */
/* "Stable" means that the rational reconstruction of the previous iteration*/
/* agrees with the new modular result (for every coefficient).             */
/* See also BM_check_with_det which is used if the determinant of the      */
/* separator rows is "stable".                                             */

static int BM_check(BMGB char0, mpz_t modulus, BM modp, mpz_t **points)
{
  FFelem tmp;
  int i, j, p, OK, failed;
  int count, pass;
  int last_fail = char0->last_fail;
  int npoints = modp->npoints;
  int GBsize = modp->GBsize;
  int *sep;
  mpz_t B;
  mpz_alias num, den;

  if (projective_case) sep = modp->projective_sep; else sep = modp->sep;
  mpz_init_set_ui(B, 1);
  i = mpz_sizeinbase(modulus, 2);
  mpz_mul_2exp(B, B, i/2);
  p = CurrentFF.prime;
  OK = 0; /* set to one(=true) only if all tests are passed */
  pass = 1;
  start_pass:
  count = 0;
  if (char0->compute_separators)
  {
    for (i=0; i < npoints; i++)
    {
      for (j=0; j < npoints; j++)
      {
        count++;
        if (pass == 1 && count < last_fail) continue;
        if (pass == 2 && count >= last_fail) continue;
        num = mpq_numref(char0->sep[i][j]);
        den = mpq_denref(char0->sep[i][j]);
        if (pass == 2)
        {
          tmp = FFmul(modp->M[sep[i]][j+npoints], mpz_fdiv_ui(den, p));
          if (tmp == mpz_fdiv_ui(num, p)) continue;
        }
        failed = modular_to_rational(num, den, B, char0->sep_mod[i][j], modulus);
        if (failed) goto tidy_up;
        tmp = FFmul(modp->M[sep[i]][j+npoints], mpz_fdiv_ui(den, p));
        if (tmp != mpz_fdiv_ui(num, p))
          goto tidy_up;
      }
    }
  }

  if (char0->compute_GB)
  {
    for (i=0; i < GBsize; i++)
    {
      for (j=0; j < npoints; j++)
      {
        count++;
        if (pass == 1 && count < last_fail) continue;
        if (pass == 2 && count >= last_fail) continue;
        num = mpq_numref(char0->GB[i][j]);
        den = mpq_denref(char0->GB[i][j]);
        if (pass == 2)
        {
          tmp = FFmul(modp->M[modp->GB[i]][j+npoints], mpz_fdiv_ui(den, p));
          if (tmp == mpz_fdiv_ui(num, p)) continue;
        }
        failed = modular_to_rational(num, den, B, char0->GB_mod[i][j], modulus);
        if (failed) goto tidy_up;
        tmp = FFmul(modp->M[modp->GB[i]][j+npoints], mpz_fdiv_ui(den, p));
        if (tmp != mpz_fdiv_ui(num, p))
          goto tidy_up;
      }
    }
  }
  if (pass == 1) { pass = 2; goto start_pass; }
  OK = 1;
  
  tidy_up:
  mpz_clear(B);
  if (!OK) { char0->last_fail = count; return 0; }
  return BM_check_final(modp, char0, modulus, p, points);
}


/***************************************************************************/
/* Check to see whether the GBasis elements are "stable".                  */
/* Return value is non-zero if they are stable, otherwise zero.            */
/* "Stable" means that the rational reconstruction of the previous iteration*/
/* agrees with the new modular result (for every coefficient).             */
/* See also BM_check which is used if the determinant of the separator rows*/
/* is not "stable".                                                        */


static int BM_check_with_det(BMGB char0, mpz_t modulus, BM modp, mpz_t **points)
{
  int i, j, OK, d, Dp;
  mpz_t gcd, modulus2;
  mpz_alias num, den, D;
  int npoints = modp->npoints;
  int nvars = modp->nvars;
  int GBsize = modp->GBsize;
  int p = CurrentFF.prime;
  int count, pass;
  int last_fail = char0->last_fail;
  int *sep;
    
  if (projective_case) sep = modp->projective_sep; else sep = modp->sep;
  D = char0->det[0];
  Dp = modp->det[0];
  mpz_init(modulus2); mpz_fdiv_q_ui(modulus2, modulus, 2);
  mpz_init(gcd);
  OK = 0; /* set to true(=1) only when all tests have been passed */
  pass = 1;
  start_pass:
  count = 0;
  if (char0->compute_separators)
  {
    for (i=0; i < npoints; i++)
    {
      d = deg(nvars, modp->pp[sep[i]]);
      for (j=0; j < npoints; j++)
      {
        count++;
        if (pass == 1 && count < last_fail) continue;
        if (pass == 2 && count >= last_fail) continue;
        num = mpq_numref(char0->sep[i][j]);
        den = mpq_denref(char0->sep[i][j]);
        if (projective_case) { D = char0->det[d]; Dp = modp->det[d]; }
        mpz_mul(num, char0->sep_mod[i][j], D);
        mpz_fdiv_r(num, num, modulus);
        if (mpz_cmp(num, modulus2) > 0) mpz_sub(num, num, modulus);
        if (FFmul(modp->M[sep[i]][j+npoints], Dp) != mpz_fdiv_ui(num, p))
          goto tidy_up;
        mpz_gcd(gcd, num, D);
        mpz_divexact(num, num, gcd);
        mpz_divexact(den, D, gcd);
      }
    }
  }
  
  if (char0->compute_GB)
  {
    for (i=0; i < GBsize; i++)
    {
      d = deg(nvars, modp->pp[modp->GB[i]]);
      for (j=0; j < npoints; j++)
      {
        count++;
        if (pass == 1 && count < last_fail) continue;
        if (pass == 2 && count >= last_fail) continue;
        num = mpq_numref(char0->GB[i][j]);
        den = mpq_denref(char0->GB[i][j]);
        if (projective_case) { D = char0->det[d]; Dp = modp->det[d]; }
        mpz_mul(num, char0->GB_mod[i][j], D);
        mpz_fdiv_r(num, num, modulus);
        if (mpz_cmp(num, modulus2) > 0) mpz_sub(num, num, modulus);
        if (FFmul(modp->M[modp->GB[i]][j+npoints], Dp) != mpz_fdiv_ui(num, p))
          goto tidy_up;
        mpz_gcd(gcd, num, D);
        mpz_divexact(num, num, gcd);
        mpz_divexact(den, D, gcd);
        if (mpz_sgn(den) < 0) { mpz_neg(num, num); mpz_neg(den, den); }
      }
    }
  }
  if (pass == 1) { pass = 2; goto start_pass; }  
  OK = 1;
  
  tidy_up:
  mpz_clear(modulus2);
  mpz_clear(gcd);
  if (!OK) { char0->last_fail = count; return 0; }
  return BM_check_final(modp, char0, modulus, p, points);
}

/***************************************************************************/
/* Chinese Remainder lifting routines for the separators and GBasis elems. */
/* cra_tmp is just used as workspace to avoid repeated malloc/free calls.  */

static void cra_GBasis(BMGB char0, BM modp, mpz_t modulus, mpz_t cra_tmp)
{
  int i, j, GBji_p;
  int p = CurrentFF.prime;
  int modulus_p = mpz_fdiv_ui(modulus, p); /* used by mpz_cra_ui_raw */
  int GBsize = modp->GBsize;
  int npoints = modp->npoints;
  
  for (j=0; j < GBsize; j++)
  {
    for (i=0; i < npoints; i++)
    {
      GBji_p = modp->M[modp->GB[j]][i+npoints];
      mpz_cra_ui_raw(char0->GB_mod[j][i], modulus, GBji_p, p, cra_tmp, modulus_p);
    }
  }
}

static void cra_separators(BMGB char0, BM modp, mpz_t modulus, mpz_t cra_tmp)
{
  int i, j, sepji_p;
  int *sep;
  int p = CurrentFF.prime;
  int modulus_p = mpz_fdiv_ui(modulus, p); /* used by mpz_cra_ui_raw */
  int npoints = modp->npoints;
  
  if (projective_case) sep = modp->projective_sep; else sep = modp->sep;
  for (j=0; j < npoints; j++)
  {
    for (i=0; i < npoints; i++)
    {
      sepji_p = modp->M[sep[j]][i+npoints];
      mpz_cra_ui_raw(char0->sep_mod[j][i], modulus, sepji_p, p, cra_tmp, modulus_p);
    }
  }
}


static void BMalg(BMGB *char0_ptr, BM *modp_ptr, int nvars, int npoints, mpz_t **points, int (*cmp)(const void *,const void*))
{
  int i, p, OK, tmp, GBsize, changed;
  FFelem **points_modp;
  FF Fp;
  BM modp;
  BMGB char0;
  mpz_t modulus, cra_tmp;

  points_modp = (FFelem**)MALLOC(npoints*sizeof(FFelem*));
  for (i=0; i < npoints; i++)
    points_modp[i] = (FFelem*)MALLOC(nvars*sizeof(FFelem));
  modp = BM_ctor(npoints, nvars, cmp);
  p = 99 + npoints*npoints;
  if (p > 20000) p = 20000;
  GBsize = 0; /* just to keep the compiler quiet */

  mpz_init_set_ui(modulus, 1);
  char0 = NULL;

  while (1)
  {
    p = nextprime(p);
    Fp = FFctor(p);
    FFselect(Fp);
    OK = reduce_mod_p(points_modp, nvars, npoints, points);
    if (!OK) { FFdtor(Fp); continue; }
    tmp = BM_mod_p(modp, points_modp);
/*printf("p=%d, tmp=%d\n", p, tmp);*/
    if (tmp < 0) { FFdtor(Fp); continue; }
    if (tmp > 0 || char0 == NULL)
    {
      GBsize = modp->GBsize;
      BMGB_dtor(char0);
      char0 = BMGB_ctor(npoints, GBsize);
      mpz_set_ui(modulus,1);
    }

    /* Compute determinant modulo p (up to sign). */
    if (projective_case)
    {
      changed = 0;
      for (i=0; i <= modp->maxdeg; i++)
      {
        if (mpz_sgn(char0->det[i]) == 0 && mpz_cmp_ui(modulus, 1)) continue;
        if (modp->det[i] == 0) { mpz_set_ui(char0->det[i], 0); continue; }
        changed |= (mpz_cra_ui(char0->det[i], modulus, modp->det[i], p) != 0);
      }
    }
    else /* affine_case */
    {
      changed = mpz_cra_ui(char0->det[0], modulus, modp->det[0], p);
    }
    const long log_modulus = mpz_sizeinbase(modulus, 2); // really just to keep compiler quiet
    if (log_modulus > char0->target_height)
    {
      if (changed)
      {
        if (BM_check(char0, modulus, modp, points)) { FFdtor(Fp); break; }
      }
      else if (BM_check_with_det(char0, modulus, modp, points)) { FFdtor(Fp); break; }
    }

    mpz_init(cra_tmp); /* workspace for multiple calls to mpz_cra_ui_raw */
    if (char0->compute_GB) cra_GBasis(char0, modp, modulus, cra_tmp);
    if (char0->compute_separators) cra_separators(char0, modp, modulus, cra_tmp);
    mpz_clear(cra_tmp);

    mpz_mul_ui(modulus, modulus, p);
    FFdtor(Fp);
  }
  mpz_clear(modulus);
  for (i=0; i < npoints; i++) FREE(points_modp[i]);
  FREE(points_modp);

  *char0_ptr = char0;
  *modp_ptr = modp;
/*  BM_dtor(modp);    it is the caller's responsibility to destroy these */
/*  BMGB_dtor(char0); .................................................. */
}

/***************************************************************************/
/*  T H I S   I S   T H E   P U B L I C   I N T E R F A C E                */
/***************************************************************************/

BM BM_affine_mod_p(int nvars, int npoints, FFelem **points, int (*cmp)(const void*, const void*))
{
  BM ans;

  if (affine_check_mod_p(nvars, npoints, points) == 0) return NULL;
  affine_case = 1;
  projective_case = 0;
  ans = BM_ctor(npoints, nvars, cmp);
  BM_mod_p(ans, points);
  return ans;
}

BM BM_projective_mod_p(int nvars, int npoints, FFelem **points, int (*cmp)(const void*, const void*))
{
  BM ans;

  if (projective_check_mod_p(nvars, npoints, points) == 0) return NULL;
  affine_case = 0;
  projective_case = 1;
  ans = BM_ctor(npoints, nvars, cmp);
  BM_mod_p(ans, points);
  return ans;
}

static int affine_check(int nvars, int npoints, mpz_t **points)
{
  int i, j, k;

  for (i=1; i < npoints; i++)
    for (j=0; j < i; j++)
    {
      for (k=0; k < nvars; k++)
        if (mpz_cmp(points[i][k], points[j][k])) break;
      if (k == nvars) return 0; /* two points are equal */
    }
  return 1; /* OK, all points are distinct */
}

static int projective_check(int nvars, int npoints, mpz_t **points)
{
  int i, j, k, ans, pivot;
  mpz_t tmp1, tmp2;
  
  ans = 0;
  mpz_init(tmp1); mpz_init(tmp2);
  for (i=0; i < npoints; i++)
  {
    for (j=0; j < nvars; j++)
      if (mpz_sgn(points[i][j]) != 0) break;
    if (j == nvars) goto end; /* disallow (0, 0, ..., 0) in projective case */
    for (j=0; j < i; j++)
    {
      for (pivot=0; pivot < nvars; pivot++)
        if (mpz_sgn(points[i][pivot]) != 0 &&
            mpz_sgn(points[j][pivot]) != 0)
          break;
      if (pivot == nvars) continue;
      for (k=0; k < nvars; k++)
      {
        mpz_mul(tmp1, points[i][k], points[j][pivot]);
        mpz_mul(tmp2, points[j][k], points[i][pivot]);
        if (mpz_cmp(tmp1, tmp2) != 0) break;
      }
      if (k == nvars) goto end; /* two projective points are equal */
    }
  }
  ans = 1; /* OK, all projective points are distinct */
  end:
  mpz_clear(tmp2); mpz_clear(tmp1);
  return ans;
}

void BM_affine(BMGB *char0_ptr, BM *modp_ptr, int nvars, int npoints, mpz_t **points, int (*cmp)(const void *,const void*))
{
  if (affine_check(nvars, npoints, points) == 0) { *char0_ptr = NULL; *modp_ptr = NULL; return; }
  affine_case = 1;
  projective_case = 0;
  BMalg(char0_ptr, modp_ptr, nvars, npoints, points, cmp);
}

void BM_projective(BMGB *char0_ptr, BM *modp_ptr, int nvars, int npoints, mpz_t **points, int (*cmp)(const void *,const void*))
{
  if (projective_check(nvars, npoints, points) == 0) { *char0_ptr = NULL; *modp_ptr = NULL; return; }
  affine_case = 0;
  projective_case = 1;
  BMalg(char0_ptr, modp_ptr, nvars, npoints, points, cmp);
}

