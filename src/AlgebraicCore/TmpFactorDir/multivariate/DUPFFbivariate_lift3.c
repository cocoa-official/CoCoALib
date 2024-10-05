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


// ************************************
// ***** EVERYTHING COMMENTED OUT *****
// ************************************
#if 0

#include <stdio.h>
#include <math.h>
// the line below for debugging only...
#include "DUPFFprint.h"
#include "DUPFFbivariate_lift3.h"
#include "jalloc.h"
#include "jaaerror.h"
#include "jtime.h"



/***************************************************************************/
/* Conversion of poly from a y-adic expansion to a (y-alpha)-adic expansion*/
/* Input is poly  = F[0] +        y *F[1] + ... +         y^deg * F[deg].  */
/* On return poly = F[0] + (y-alpha)*F[1] + ... + (y-alpha)^deg * F[deg].  */


void DUPFFbiv_shift(const DUPFF *F, int deg, int alpha)
{
  int i, j;

  for (i=0; i < deg; i++)
    for (j=deg-1; j >= i; j--)
      DUPFFshift_add(F[j],F[j+1],0,alpha);/* F[j] += alpha * F[j+1] */
}


DUPFF *DUPFFbivariate_transpose(DUPFF *f, int dy)
{
  DUPFF *ans;
  int dx, i, j;

  dx = 0;
  for (i=0; i <= dy; i++) if (DUPFFdeg(f[i]) > dx) dx = DUPFFdeg(f[i]);
  ans = (DUPFF*)MALLOC((1+dx)*sizeof(DUPFF));
  for (i=0; i <= dx; i++) ans[i] = DUPFFnew(dy);
  for (i=0; i <= dx; i++)
  {
    for (j=0; j <= dy; j++)
    {
      if (DUPFFdeg(f[j]) < i) ans[i]->coeffs[j] = 0;
      else ans[i]->coeffs[j] = f[j]->coeffs[i];
    }
    for (j=dy; j >= 0; j--) if (ans[i]->coeffs[j] != 0) break;
    ans[i]->deg = j;
  }
  return ans;
}

DUPFF DUPFFbivariate_content(DUPFF *f, int dy)
{
  int i;
  DUPFF content;

printf("content called on:\n");
for (i=0; i <= dy; i++) DUPFFprint(f[i]);
content=DUPFFnew(0);
  for (i=0; i <= dy; i++)
  {
    content = DUPFFgcd(content, f[i]); /* LEAK LEAK BUG BUG BUG */
    if (DUPFFdeg(content) == 0) break;
  }
DUPFFmake_monic(content);
printf("found content=");DUPFFprint(content);
return content;
}

void DUPFFbivariate_primitive(DUPFF *f, int dy)
{
  DUPFF *tr, *tr2, content, zero;
  int dx, i;

  zero = DUPFFnew(0);
  dx = DUPFFdeg(f[0]); /* CHEAT */
  tr = DUPFFbivariate_transpose(f, dy);
  content = DUPFFbivariate_content(tr, dx);
  for (i=0; i <= dx; i++) DUPFFdiv4(tr[i], zero, tr[i], content);
  tr2 = DUPFFbivariate_transpose(tr, dx);
  for (i=0; i <= dy-DUPFFdeg(content); i++) DUPFFcopy2(f[i], tr2[i]);
}

/***************************************************************************/
/* Lift a univariate factorization to a bivariate one.                     */
/* Input is f a y-adic expansion of the bivariate polynomial,              */
/* and the univariate factorization of f at y=alpha.                       */
/* Output is array such that lifted_factors[i] is y-adic expn of lifted factor*/

DUPFF** DUPFFbivariate_lift(DUPFF *f, int dy, DUPFFlist factors, unsigned int alpha)
{
  DUPFF **lifted_factors, *newf;
  DUPFFlist lcs, iter;
  DUPFF scale, lcf;
  unsigned int lc_image;
  int dlcf;
  DUPFFbiv_lifter lifter;
  int i, j, nfactors, deg, dx, newdy;
  int p = CurrentFF.prime;

printf("shifting...\n");print_time();
  if (alpha != 0)
    DUPFFbiv_shift(f, dy, alpha); /* convert to (y-alpha)-adic expansion */
printf("shifted.\n");print_time();
nfactors = DUPFFlist_length(factors);
dx=DUPFFdeg(f[0]);
lcf=DUPFFnew(dy);
for(i=0;i<=dy;i++)if(DUPFFdeg(f[i])==dx)lcf->coeffs[i]=f[i]->coeffs[dx];else lcf->coeffs[i]=0;
for(i=dy;i>=0;i--)if(lcf->coeffs[i]!=0)break;lcf->deg=i;
printf("leading coeff is: ");DUPFFprint(lcf);
dlcf=DUPFFdeg(lcf);
scale=DUPFFexpt(lcf, nfactors - 1);
newdy=dy+dlcf*(nfactors-1);
printf("newf has degree %d.\n", newdy);
newf=(DUPFF*)MALLOC((1+newdy)*sizeof(DUPFF));
for(i=0;i<=newdy;i++)newf[i]=DUPFFnew(dx);
for(i=0;i<=dy;i++)for(j=0;j<=DUPFFdeg(scale);j++)DUPFFshift_add(newf[i+j],f[i], 0, lcf->coeffs[j]);
printf("newf is:\n");for(i=0;i<=newdy;i++){printf("y^%d*",i);DUPFFprint(newf[i]);}
lcs=NULL;
for(i=0;i<nfactors;i++) lcs = DUPFFlist_append(DUPFFlist_ctor(lcf, 1), lcs);
printf("lcs are:\n");
for(iter=lcs;iter;iter=iter->next)DUPFFprint(iter->poly);
lc_image=DUPFFeval(lcf,alpha);
for(iter=factors;iter;iter=iter->next) DUPFFmul2ff(iter->poly,FFdiv(lc_image,DUPFFlc(iter->poly)));
  lifter = DUPFFbiv_lift_init(newf, newdy, factors, lcs);
  for (i=1; i <= newdy/2; i++)
  {if(i%1==0){printf("level=%d\n",i);print_time();}
    DUPFFbiv_lift_step(lifter, newf[i]);
  }
printf("lifting completed (there may be 1 bogus factor)\n");print_time();
  lifted_factors = (DUPFF **)MALLOC(2*sizeof(DUPFF *));
  nfactors = DUPFFbiv_lift_output(lifted_factors, lifter, 0);
  DUPFFbiv_lifter_dtor(lifter);
printf("first factor:\n");
  for (i=0; i <= newdy/2; i++)
  {
    printf("part of degree %d in y:", i);
    DUPFFprint(lifted_factors[0][i]);
  }
for(i=0;i<nfactors;i++)DUPFFbivariate_primitive(lifted_factors[i], newdy/2);  

printf("factors output\n");print_time();  
  if (alpha != 0)
    for (i=0; i < nfactors; i++)
    {
      for (deg=newdy/2; deg > 0; deg--)
        if (DUPFFdeg(lifted_factors[i][deg]) >= 0) break;
      DUPFFbiv_shift(lifted_factors[i], deg, p-alpha);
    }
printf("first factor:\n");
  for (i=0; i <= newdy/2; i++)
  {
    printf("part of degree %d in y:", i);
    DUPFFprint(lifted_factors[0][i]);
  }
printf("factors shifted\n");print_time();  
  return lifted_factors;
}



/***************************************************************************/

void solve_for_correction(DUPFF *gdelta, const DUPFF E, const DUPFF g0, const DUPFF grecip)
{
  DUPFF junk, tmp;
  int dg0, dE;

  dg0 = DUPFFdeg(g0);
  dE = DUPFFdeg(E);
  junk = DUPFFnew(dE + dg0);  /* these are all a bit too big */
  tmp = DUPFFnew(2*dg0);      /*                             */
  *gdelta = DUPFFnew(dg0);    /* extra space necessary for lc fix */

  DUPFFdiv4(junk, tmp, E, g0);
  DUPFFmul3(tmp, tmp, grecip);
  DUPFFdiv4(junk, *gdelta, tmp, g0);

  DUPFFfree(tmp);
  DUPFFfree(junk);
}

/***************************************************************************/


static DUPFFbiv_lifter DUPFFbiv_lifter_ctor(DUPFF g, DUPFF h, DUPFFbiv_lifter g_lifter, DUPFFbiv_lifter h_lifter, DUPFF lcg, DUPFF lch, int rmax)
{
  DUPFFbiv_lifter ans;
  DUPFF grecip, hrecip, tmp;
  int df;
  
//  if (DUPFFdeg(g) > DUPFFdeg(h)) return DUPZfactor_lifter_ctor(h, g, h_lifter, g_lifter, lch, lcg, r);
  ans = (DUPFFbiv_lifter)MALLOC(sizeof(struct DUPFFbiv_lift_struct));
  ans->r = 1;
  ans->g = (DUPFF*)MALLOC(rmax*sizeof(DUPFF));
  ans->h = (DUPFF*)MALLOC(rmax*sizeof(DUPFF));
  ans->g[0] = DUPFFcopy(g);
  ans->lcg = lcg;
  ans->h[0] = DUPFFcopy(h);
  ans->lch = lch;
  df = DUPFFdeg(g) + DUPFFdeg(h);
  ans->E = DUPFFnew(df);
  
  tmp = DUPFFexgcd(&grecip, &hrecip, g, h);
  if (DUPFFdeg(tmp) != 0) { JERROR(JERROR_HENSEL); exit(1); }
  /* since the gcd is not monic we have to do this... */
  DUPFFdiv2ff(grecip, DUPFFlc(tmp));
  DUPFFdiv2ff(hrecip, DUPFFlc(tmp));
  ans->grecip = grecip;
  ans->hrecip = hrecip;
  ans->g_lifter = g_lifter;
  ans->h_lifter = h_lifter;
  DUPFFfree(tmp);
  return ans;
}



void DUPFFbiv_lifter_dtor(DUPFFbiv_lifter THIS)
{
//  int i;
  int r;

  if (THIS == NULL) return;
  r = THIS->r;
  DUPFFfree(THIS->E);
  if (THIS->g_lifter) DUPFFbiv_lifter_dtor(THIS->g_lifter);
  if (THIS->h_lifter) DUPFFbiv_lifter_dtor(THIS->h_lifter);
//  for (i=0; i < r; i++)
//  {
//    DUPFFfree(THIS->g[i]);
//    DUPFFfree(THIS->h[i]);
//  }
//  FREE(THIS->g);
//  FREE(THIS->h);
  DUPFFfree(THIS->grecip);
  DUPFFfree(THIS->hrecip);
  FREE(THIS);
}

/***************************************************************************/

void DUPFFbiv_lift_step(DUPFFbiv_lifter THIS, DUPFF f_cmpt_r)
{
  int i;
  /* the rest are just aliases to make the code easier to read (!?!) */
  DUPFF       *g     = THIS->g;
  DUPFF       *h     = THIS->h;
  DUPFF       lcg    = THIS->lcg;
  DUPFF       lch    = THIS->lch;
  int         r      = THIS->r;
  const DUPFF grecip = THIS->grecip;
  const DUPFF hrecip = THIS->hrecip;
  DUPFF       E      = THIS->E;
DUPFF tmp;

  tmp = DUPFFnew(DUPFFdeg(g[0]) + DUPFFdeg(h[0]));
  DUPFFcopy2(E, f_cmpt_r);
  for (i=1; i < r; i++)
  {
    DUPFFmul3(tmp, g[i], h[r-i]);
    DUPFFsub3(E, E, tmp);
  }
  DUPFFfree(tmp);

//printf("before solve g[0]=");DUPFFprint(g[0]);
  solve_for_correction(&g[r], E, g[0], hrecip);
  if (r <= DUPFFdeg(lcg) && lcg->coeffs[r] != 0)
    DUPFFshift_add(g[r], g[0], 0, lcg->coeffs[r]);
    /* g[r] += lcg->coeffs[r]*g[0]; */
//printf("after solve g[r]=");DUPFFprint(g[r]);
//printf("before solve h[0]=");DUPFFprint(h[0]);
  solve_for_correction(&h[r], E, h[0], grecip);
  if (r <= DUPFFdeg(lch) && lch->coeffs[r] != 0)
    DUPFFshift_add(h[r], h[0], 0, lch->coeffs[r]);
    /* h[r] += lch->coeffs[r]*h[0]; */
//printf("after solve h[r]=");DUPFFprint(h[r]);
  THIS->r++;

  if (THIS->g_lifter) DUPFFbiv_lift_step(THIS->g_lifter, g[r]);
  if (THIS->h_lifter) DUPFFbiv_lift_step(THIS->h_lifter, h[r]);
}


/***************************************************************************/
/* Make the lifted factors available to the caller.                        */

int DUPFFbiv_lift_output(DUPFF **factors, const DUPFFbiv_lifter THIS, int i)
{
  int p, k;

  p = CurrentFF.prime;
  k = THIS->r;
  if (THIS->g_lifter) i = DUPFFbiv_lift_output(factors, THIS->g_lifter, i);
  else factors[i++] = THIS->g;
  if (THIS->h_lifter) i = DUPFFbiv_lift_output(factors, THIS->h_lifter, i);
  else factors[i++] = THIS->h;
  return i;
}



/***************************************************************************/

static DUPFFbiv_lifter DUPFFbiv_lift_init_stupid(const DUPFFlist pfactors, const DUPFFlist lcs, int rmax)
{
  int df;
  DUPFFbiv_lifter ans;
  DUPFF tmp, tmp_lc;
  DUPFFlist iter, lc;

  df = 0;
printf("factors are:\n");
  for (iter = pfactors; iter; iter = iter->next)
  {
    df += DUPFFdeg(iter->poly);
DUPFFprint(iter->poly);
  }
  ans = NULL;
  tmp = DUPFFnew(df);
  DUPFFcopy2(tmp, pfactors->poly);
  tmp_lc = lcs->poly;
  for (iter = pfactors->next, lc = lcs->next; iter; iter = iter->next, lc = lc->next)
  {
    ans = DUPFFbiv_lifter_ctor(iter->poly, tmp, NULL, ans, lc->poly, tmp_lc, rmax);
    DUPFFmul3(tmp, tmp, iter->poly);
    tmp_lc = DUPFFmul(tmp_lc, lc->poly);
  }
  DUPFFfree(tmp);
  DUPFFfree(tmp_lc);
  return ans;
}


DUPFFbiv_lifter DUPFFbiv_lift_init(DUPFF *f, int dy, const DUPFFlist pfactors, const DUPFFlist lcs)
{
  return DUPFFbiv_lift_init_stupid(pfactors, lcs, dy);
}



#if 0
///****************************************************************************/
///* This function is called after a speculative early search which has found */
///* at least one factor but further lifting is still required.               */
///* It creates a new lifting tree for the remaining factors.                 */
//
//void DUPFFbiv_lift_revise1(DUPFFfactor_info info, int height)
//{
//  int i;
//
//  DUPFFbiv_lift_init(info);
//  DUPFFpadic_expand_end(info->f_expansion);
//  info->f_expansion = DUPFFpadic_expand_init(info->f, info->p, info->rmax);
//  DUPFFpadic_expand_step(info->f_expansion);
//  if (height > info->rmax) height = info->rmax;
//  for(i=1; i < height; i++)
//  {
//    DUPFFpadic_expand_step(info->f_expansion);
//    DUPFFbiv_lift_step(info->lifter, info->f_expansion->monic_cmpts[i]);
//  }
//}
//
//
//void DUPFFbiv_lift_revise(DUPFFfactor_info info)
//{
//  int height = info->lifter->r;
//
//  DUPFFbiv_lift_dtor(info->lifter); info->lifter = NULL;
//  DUPFFpadic_expand_end(info->f_expansion); info->f_expansion = NULL;
//  DUPFFbiv_lift_revise1(info, height);
//}
#endif

/***************************************************************************/
void liftprint(DUPFFbiv_lifter THIS, int indent)
{
int i;
if (THIS == NULL) return;
putchar('>');for(i=0;i<indent;i++)putchar(' ');
for(i=0;i<THIS->r;i++)printf("%d,",DUPFFdeg(THIS->g[i]));printf("\n");
for(i=0;i<=indent;i++)putchar(' ');
for(i=0;i<THIS->r;i++)printf("%d,",DUPFFdeg(THIS->h[i]));printf("\n");
liftprint(THIS->g_lifter, indent+2);
liftprint(THIS->h_lifter, indent+2);
}

#endif  // end of commented out section
