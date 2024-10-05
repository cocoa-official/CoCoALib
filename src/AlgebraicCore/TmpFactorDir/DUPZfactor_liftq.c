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
#include <stdio.h>
#include "DUPZfactor_liftq.h"
#include "jalloc.h"
#include "jaaerror.h"
#include "DUPZ_DUPFF.h"
#include "DUPZfactor_info.h"


/***************************************************************************/


DUPZfactor_lifter DUPZfactor_lifter_ctor(DUPFF g, DUPFF h, DUPZfactor_lifter g_lifter, DUPZfactor_lifter h_lifter)
{
  DUPZfactor_lifter ans;
  DUPFF corr_factor, junk, tmp;
  int df;
  
  if (DUPFFdeg(g) > DUPFFdeg(h)) return DUPZfactor_lifter_ctor(h, g, h_lifter, g_lifter);
  ans = (DUPZfactor_lifter)MALLOC(sizeof(struct DUPZfactor_lift_struct));
  ans->g = DUPFF_to_DUPZ(g);
  ans->h = DUPFF_to_DUPZ(h);
  df = DUPFFdeg(g) + DUPFFdeg(h);
  tmp = DUPFFexgcd(&junk, &corr_factor, g, h);
  if (DUPFFdeg(tmp) != 0) { JERROR(JERROR_HENSEL); exit(1); }
  /* since the gcd is not monic we have to do this... */
  DUPFFdiv2ff(corr_factor, DUPFFlc(tmp));
  ans->corr_factor = DUPZnew(df-1); /* need spare space for revise_lift_tree */
  DUPFF_to_DUPZ2(ans->corr_factor, corr_factor);
  ans->g_lifter = g_lifter;
  ans->h_lifter = h_lifter;
  DUPFFfree(tmp);
  DUPFFfree(corr_factor);
  DUPFFfree(junk);

  return ans;
}


void DUPZfactor_lift_dtor(DUPZfactor_lifter THIS)
{
  if (THIS == NULL) return;
  if (THIS->g_lifter) DUPZfactor_lift_dtor(THIS->g_lifter);
  if (THIS->h_lifter) DUPZfactor_lift_dtor(THIS->h_lifter);
  DUPZfree(THIS->g);
  DUPZfree(THIS->h);
  DUPZfree(THIS->corr_factor);
  FREE(THIS);
}

/***************************************************************************/

void DUPZfactor_lift_step(DUPZfactor_lifter THIS, DUPZ f, mpz_t Q)
{
  /* local workspace -- maybe put these in lifter struct??? */
  DUPZ Delta, delta_g, tmp;
  /* just aliases to make the code easier to read (!?!) */
  DUPZ       g     = THIS->g;
  DUPZ       h     = THIS->h;
  DUPZ       alpha = THIS->corr_factor;

  Delta = DUPZmul(g, h);
  DUPZsub3(Delta, f, Delta);
  DUPZdiv2z(Delta, Q);
  DUPZmod2z(Delta, Q);

  tmp = DUPZnew(DUPZdeg(f));
  DUPZmul3(tmp, alpha, h);
  mpz_sub_ui(tmp->coeffs[0], tmp->coeffs[0], 2);
  DUPZmod2(tmp, g, Q);
  DUPZmul3(tmp, alpha, tmp);
  DUPZneg1(tmp);
  DUPZmod2(tmp, g, Q);
  DUPZcopy2(alpha, tmp);

  delta_g = DUPZnew(DUPZdeg(f));
  DUPZcopy2(delta_g, Delta);
  DUPZmod2(delta_g, g, Q);
  DUPZmul3(delta_g, alpha, delta_g);
  DUPZmod2(delta_g, g, Q);
  DUPZmul3(tmp, h, delta_g);
  DUPZsub3(tmp, Delta, tmp);
  DUPZmdiv2(tmp, g, Q);
  DUPZmul2z(delta_g, Q);
  DUPZmul2z(tmp, Q);
  DUPZadd3(g, g, delta_g);
  DUPZadd3(h, h, tmp);
  DUPZfree(Delta);DUPZfree(delta_g);DUPZfree(tmp);
  
  if (THIS->g_lifter) DUPZfactor_lift_step(THIS->g_lifter, g, Q);
  if (THIS->h_lifter) DUPZfactor_lift_step(THIS->h_lifter, h, Q);
}


/***************************************************************************/

static DUPZfactor_lifter DUPZfactor_lift_init_stupid(const DUPFFlist pfactors)
{
  DUPZfactor_lifter ans;
  DUPFFlist iter;
  DUPFF tmp;
  int df;

  df = 0;
  for (iter = pfactors; iter; iter = iter->next)
    df += DUPFFdeg(iter->poly);
  ans = NULL;
  tmp = DUPFFnew(df);
  DUPFFcopy2(tmp, pfactors->poly);
  for (iter = pfactors->next; iter; iter = iter->next)
  {
    ans = DUPZfactor_lifter_ctor(iter->poly, tmp, NULL, ans);
    DUPFFmul3(tmp, tmp, iter->poly);
  }
  DUPFFfree(tmp);

  return ans;
}


void DUPZfactor_lift_init(DUPZfactor_info THIS)
{
  THIS->lifter = DUPZfactor_lift_init_stupid(THIS->pfactors);
}


#if 0
/****************************************************************************/
/* This function is called after a speculative early search which has found */
/* at least one factor but further lifting is still required.               */
/* It creates a new lifting tree for the remaining factors.                 */

static DUPZfactor_lifter revise_lifter(DUPZfactor_lifter old_tree, mpz_t Q)
{
?????????????????????????????????????????????????????????????????????????????  
  if (old_tree->g_lifter != NULL) g1 = revise_lifter(old_tree->g_lifter, Q);
  if (old_tree->h_lifter != NULL) h1 = revise_lifter(old_tree->h_lifter, Q);
  ans->g = DUPZmul(g1->g, h1->g);
  ans->h = DUPZmul(g1->h, h1->h);
  ans->g_lifter = NULL;
  
  if (DUPZdeg(g1->g) == 0 || DUPZdeg(h1->g) == 0) goto tidy_up;
  if (DUPZdeg(g1->g) > DUPZdeg(h1->g))
  {
    new_tree->g = g1->g; g1->g = NULL;
    new_tree->h = h1->g; h1->g = NULL;
    new_tree->g_lifter = g1->g_lifter;
    new_tree->h_lifter = h1->g_lifter;
    new_tree->corr_factor = DUPZnew(DUPZdeg(ans->g));
    DUPZmul3(new_tree->corr_factor, old_tree->corr_factor, h1->h);
    DUPZmod3(new_tree->corr_factor, new_tree->g, Q);
  }
  else
  {
    new_tree->g = h1->g; h1->g = NULL;
    new_tree->h = g1->g; g1->g = NULL;
    new_tree->g_lifter = h1->g_lifter;
    new_tree->h_lifter = g1->g_lifter;
    new_tree->corr_factor = DUPZnew(df-1);
    DUPZmul3(new_tree->corr_factor, old_tree->corr_factor, g1->h);
    DUPZmod3(new_tree->corr_factor, new_tree->h, Q);
    DUPZmul3(new_tree->corr_factor, new_tree->corr_factor, new_tree->g);
    DUPZmod3(new_tree->corr_factor, new_tree->h, Q);
    DUPZneg(new_tree->corr_factor);
    DUPZadd3(new_tree->corr_factor, new_tree->corr_factor, 1);
    DUPZdivmodQ(new_tree->corr_factor, new_tree->h);
  }
  tidy_up:
  DUPZfactor_lift_dtor(g1);
  DUPZfactor_lift_dtor(h1);
  return ans;
}
#endif

static DUPZ revise_lift_tree(DUPZfactor_lifter *tree, mpz_t Q)
{
  DUPZ G, H, ans;
  int g_prune, h_prune;
  DUPZfactor_lifter THIS = *tree;     /* alias */
  DUPZ g = THIS->g;                /* alias */
  DUPZ h = THIS->h;                /* alias */
  DUPZ alpha = THIS->corr_factor;  /* alias */

  if (THIS->g_lifter != NULL)
    G = revise_lift_tree(&THIS->g_lifter, Q);
  else
  {
    if (DUPZdeg(g) > 0) G = int_to_DUPZ(1);
    else
    {
      g->deg = - g->deg;
      G = DUPZcopy(g); /* wasteful */
    }
  }
  g_prune = (DUPZdeg(G) == DUPZdeg(g));

  if (THIS->h_lifter != NULL)
    H = revise_lift_tree(&THIS->h_lifter, Q);
  else
  {
    if (DUPZdeg(h) > 0) H = int_to_DUPZ(1);
    else
    {
      h->deg = - h->deg;
      H = DUPZcopy(h); /* wasteful */
    }
  }
  h_prune = (DUPZdeg(H) == DUPZdeg(h));

  if (g_prune && h_prune)
  {
    DUPZfactor_lift_dtor(THIS);
    *tree = NULL;
    goto end;
  }
  if (g_prune)
  {
    *tree = THIS->h_lifter;
    THIS->h_lifter = NULL;
    DUPZfactor_lift_dtor(THIS);
    goto end;
  }
  if (h_prune)
  {
    *tree = THIS->g_lifter;
    THIS->g_lifter = NULL;
    DUPZfactor_lift_dtor(THIS);
    goto end;
  }
  DUPZmdiv2(g, G, Q);
  DUPZmdiv2(h, H, Q);
  DUPZmod2(alpha, g, Q);
  DUPZmul3(alpha, alpha, H);
  DUPZmod2(alpha, g, Q);

 end:
  ans = DUPZmul(G, H);
  DUPZfree(G); DUPZfree(H);
  return ans;
}


void DUPZfactor_lift_revise(DUPZfactor_info info)
{
  DUPZfree(revise_lift_tree(&info->lifter, info->Q));
}


/***************************************************************************/
/* Make the lifted factors available to outsiders.                         */
/* We actually return pointers to our own copies so that outsiders can     */
/* modify them: they should negate the degree if the factor was used to    */
/* form a true factor (see DUPZfactor_lift_revise above).                  */

int DUPZfactor_lift_output(DUPZ **factors, const DUPZfactor_lifter THIS, int i)
{
  if (THIS->g_lifter) i = DUPZfactor_lift_output(factors, THIS->g_lifter, i);
  else factors[i++] = &THIS->g;
  if (THIS->h_lifter) i = DUPZfactor_lift_output(factors, THIS->h_lifter, i);
  else factors[i++] = &THIS->h;
  return i;
}

/***************************************************************************/
