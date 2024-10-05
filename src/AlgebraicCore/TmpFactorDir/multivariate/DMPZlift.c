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

#include "DMPZlift.h"
#include "DMPZmap_to_univariate.h"
#include "DMPZeval.h"
#include "DMPZgcd.h"
#include "DMPZ_to_DUPZ.h"
#include "DUPZ.h"
#include "DUPZgcd.h"
#include "DUPZexgcd.h"
#include "jalloc.h"
#include "jaaerror.h"
#include "Zmat.h"
#include "Qmat.h"
#include "Qsolve.h"

DMPZlifter DMPZlifter_ctor(DMPZ f, DUPZ g1, DUPZ h1, int *substitution, int x)
{
  int i;
  
  DMPZlifter THIS = (DMPZlifter)MALLOC(sizeof(struct DMPZlifter_struct));
  THIS->f = f; /* RISKY */
  THIS->flcf = NULL;
  THIS->x = x;
  THIS->degs = (int*)MALLOC(NVARS*sizeof(int));
  DMPZdegs(THIS->degs, f);
  THIS->substitution = (int*)MALLOC(NVARS*sizeof(int));
  for (i=0; i < NVARS; i++) THIS->substitution[i] = substitution[i];
  THIS->lifted = (int*)MALLOC(NVARS*sizeof(int));
  /* Mark absent variables as having been lifted... */
  for (i=0; i < NVARS; i++) THIS->lifted[i] = (THIS->degs[i] == 0);
  THIS->lifted[x] = 1;

  THIS->g1 = DUPZnew(0);
  THIS->g = DUPZ_to_DMPZ(g1, x);
  THIS->g2 = NULL;
  THIS->lcg = DMPZlc(f, x);
  THIS->g_skeleton = NULL;

  THIS->h1 = DUPZnew(0);
  THIS->h = DUPZ_to_DMPZ(h1, x);
  THIS->h2 = NULL;
  THIS->lch = DMPZlc(f, x);
  THIS->h_skeleton = NULL;

  THIS->FAILED = -1; /* used to trigger the call to DMPZlift_check_initial_subst */
                     /* see the first few lines of DMPZlift(...).                */

  return THIS;
}

void DMPZlifter_dtor(DMPZlifter THIS)
{
  DMPZdtor(THIS->flcf);/*??????????????*/
  
  DUPZfree(THIS->g1);
  DMPZdtor(THIS->g);
  DMPZdtor(THIS->g2);
  DMPZdtor(THIS->lcg);
  DMPZdtor(THIS->g_skeleton);

  DUPZfree(THIS->h1);
  DMPZdtor(THIS->h);
  DMPZdtor(THIS->h2);
  DMPZdtor(THIS->lch);
  DMPZdtor(THIS->h_skeleton);

  FREE(THIS->lifted);
  FREE(THIS->substitution);
  FREE(THIS->degs);

  FREE(THIS);
}

/***************************************************************************/

static DMPZ DMPZbivariate_to_DMPZ(DUPZ *f, int dy, int x, int y)
{
  DMPZ ans = NULL;
  int i, j, k;
  int *exps;

  for (j=0; j <= dy; j++)
  {
    if (DUPZdeg(f[j]) < 0) continue;
    for (i=0; i <= DUPZdeg(f[j]); i++)
    {
      if (mpz_sgn(f[j]->coeffs[i]) == 0) continue;
      exps = (int*)MALLOC(NVARS*sizeof(int));
      for (k=0; k < NVARS; k++) exps[k] = 0;
      exps[x] = i;
      exps[y] = j;
      ans = DMPZprepend(f[j]->coeffs[i], exps, ans);
    }
  }
  return ans;
}


/***************************************************************************/
/* Lex comparison of two power products but ignoring the exponent of the   */
/* indeterminate x.                                                        */

static int cmp_except(int *PP1, int *PP2, int x)
{
  int i;
  for (i=0; i < NVARS; i++)
  {
    if (i == x) continue;
    if (PP1[i] != PP2[i]) return PP1[i]-PP2[i];
  }
  return 0;
}

/***************************************************************************/
/* Compute the "skeleton" of a multivariate polynomial with respect to x   */
/* The skeleton is effectively the set of power products which remain after*/
/* the variable x is replaced by a generic value (i.e. so that no coeffs   */
/* cancel out).  Here we represent it as a polynomial being the sum of     */
/* these power products (with coeff = 1).                                  */

static DMPZ make_skeleton(const DMPZ f, int x)
{
  DMPZ iter, skel_iter;
  DMPZ skeleton = NULL;
  int i, *exps;
  mpz_t one;
  mpz_init_set_ui(one, 1);

  for (iter=f; iter; iter = iter->next)
  {
    for (skel_iter=skeleton; skel_iter; skel_iter = skel_iter->next)
      if (cmp_except(iter->exps, skel_iter->exps, x) == 0) break;
    if (skel_iter != NULL) continue;
    exps = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++) exps[i] = iter->exps[i];
    exps[x] = 0;
    skeleton = DMPZprepend(one, exps, skeleton);
  }
  mpz_clear(one);
  return skeleton;
}


/***************************************************************************/
/* Map a multivariate polynomial to a bivariate one by evaluation.         */
/* The given indeterminates (x and y) are retained in the bivariate image, */
/* all others are mapped to integer values given by substitution.          */
/* The result is a DMPZ for the image of f.                                */

static DMPZ DMPZmap_to_bivariate(const DMPZ f, int x, int y, int *substitution)
{
  int *flag, i;
  DMPZ ans;
  
  flag = (int*)MALLOC(NVARS*sizeof(int));
  for (i=0; i < NVARS; i++)
    flag[i] = 1;
  flag[x] = flag[y] = 0;
  ans = DMPZeval_partial(f, substitution, flag);
  FREE(flag);
  return ans;
}


/***************************************************************************/
/* Convert a bivariate DMPZ into a C vector of DUPZs.                      */
/* The variable in the DUPZs corresponds to x, while the i-th entry in the */
/* vector corresponds to those terms divisible exactly by y^i.             */
/* Do not pass in a DMPZ involving variable other than x and y!            */

static DUPZ* DMPZbivariate_padic_expand1(const DMPZ f, int x, int y)
{
  int i, dx, dy;
  DUPZ *ans;
  DMPZ iter;

  dx = DMPZdeg(f, x);
  dy = DMPZdeg(f, y);
  ans = (DUPZ*)MALLOC((1+dy)*sizeof(DUPZ));
  for (i=0; i <= dy; i++)
    ans[i] = DUPZnew(dx);
  for (iter=f; iter; iter = iter->next)
  {
    int expx = iter->exps[x];
    int expy = iter->exps[y];
    mpz_add(ans[expy]->coeffs[expx], ans[expy]->coeffs[expx], iter->coeff);
  }
  /* Now fix the degree of each entry in the vector. */
  for (i=0; i <= dy; i++)
  {
    int deg = dx;
    for (deg = dx; deg >= 0; deg--)
      if (mpz_sgn(ans[i]->coeffs[deg])) break;
    ans[i]->deg = deg;
  }
  return ans;
}


/***************************************************************************/
/* Conversion of poly from a y-adic expansion to a (y-alpha)-adic expansion*/
/* Input is poly  = F[0] +        y *F[1] + ... +         y^deg * F[deg].  */
/* On return poly = F[0] + (y-alpha)*F[1] + ... + (y-alpha)^deg * F[deg].  */

static void DMPZbivariate_padic_shift(const DUPZ *F, int deg, mpz_t alpha)
{
  int i, j;

  for (i=deg-1; i >= 0; i--)
    for (j=i; j < deg; j++)
      DUPZshift_add(F[j], F[j+1], 0, alpha);/* F[j] += alpha * F[j+1] */
}


/***************************************************************************/
/* Compute the (y-alpha)-adic expansion of the bivariate polynomial f.     */
/* x is the other variable involved in f  -- f must not involve variables  */
/* other than x and y!!                                                    */
/* Result is a C vector of DUPZs.  Entry j is the (y-alpha)^j component of */
/* expansion (as a DUPZ representing a polynomial in x).                   */

static DUPZ* DMPZbivariate_padic_expand(const DMPZ f, int x, int y, mpz_t alpha)
{
  DUPZ *ans;

  ans = DMPZbivariate_padic_expand1(f, x, y);
  DMPZbivariate_padic_shift(ans, DMPZdeg(f, y), alpha); /* or -alpha ??? */
  return ans;
}


/***************************************************************************/
/* This lifts two univariate factors to bivariate factors.                 */
/* The bivariate polynomial is mapped to univariate by y |-> a             */
/* Using linear lifting...                                                 */
/* The return value is non-zero if the lifting succeeded, and 0 otherwise. */

static int DMPZlift_to_bivariate(DMPZlifter THIS)
{
  DUPZ g_cofac, h_cofac, tmp, tmp2, gcd, corr, error;
  DUPZ lcg, lch;
  int dx, dy, level, j, max_height;
  DUPZ *f2, *glift, *hlift;
  mpz_t alpha;

  {
    DMPZ lcf = DMPZlc(THIS->f, THIS->x);
    mpz_t lcf1, scale;
    /* Ensure that the leading coeffs of the univariate images are "right" */

    mpz_init(lcf1);
    mpz_init(scale);
    DMPZeval(lcf1, lcf, THIS->substitution);
    mpz_divexact(scale, lcf1, DUPZlc(THIS->g1));
    DUPZmul2z(THIS->g1, scale);
    mpz_divexact(scale, lcf1, DUPZlc(THIS->h1));
    DUPZmul2z(THIS->h1, scale);
    mpz_clear(scale);
    mpz_clear(lcf1);
    DMPZdtor(lcf);
  }
#ifdef FACTOR_DEBUG
printf("Starting lift to bivariate.\n");
printf("Initial univariate factors are:\n");
DUPZprint(THIS->g1);
DUPZprint(THIS->h1);
printf("Variable being lifted is x[%d]\n", THIS->y);
printf("Evaluation of x[%d] is %d.\n", THIS->y, THIS->substitution[THIS->y]);
#endif

  dx = DMPZdeg(THIS->f, THIS->x);
  if (THIS->flcf == NULL)
  {
    DMPZ lcf;
    lcf = DMPZlc(THIS->f, THIS->x);
    THIS->flcf = DMPZmul(THIS->f, lcf);
    DMPZdtor(lcf);
  }
  {
    DMPZ flcf2;
    flcf2 = DMPZmap_to_bivariate(THIS->flcf, THIS->x, THIS->y, THIS->substitution);
    dy = DMPZdeg(THIS->flcf, THIS->y);
    if (DMPZdeg(flcf2, THIS->y) != dy) { DMPZdtor(flcf2); return 0; }
    mpz_init_set_si(alpha, THIS->substitution[THIS->y]);
    f2 = DMPZbivariate_padic_expand(flcf2, THIS->x, THIS->y, alpha);
    DMPZdtor(flcf2);
  }

  lcg = DMPZmap_to_univariate(THIS->lcg, THIS->y, THIS->substitution);
#ifdef FACTOR_DEBUG
printf("BIVARIATE: forcing ldcf of g (wrt x[%d]) to be ", THIS->x);DUPZprint(lcg);
#endif
  DUPZlinear_shift(lcg, alpha);
//printf("Same thing but as a poly in (x[%d]-%d): ",THIS->y, THIS->substitution[THIS->y]);DUPZprint(lcg);  
  lch = DMPZmap_to_univariate(THIS->lch, THIS->y, THIS->substitution);
  DUPZlinear_shift(lch, alpha);
//printf("lch as a poly in (x[%d]-%d): ",THIS->y, THIS->substitution[THIS->y]);DUPZprint(lch);

  max_height = dy;
  glift = (DUPZ*)MALLOC((1+max_height)*sizeof(DUPZ));
  hlift = (DUPZ*)MALLOC((1+max_height)*sizeof(DUPZ));
  glift[0] = DUPZcopy(THIS->g1);
  hlift[0] = DUPZcopy(THIS->h1);
  
  g_cofac = DUPZnew(DUPZdeg(THIS->h1)-1);
  h_cofac = DUPZnew(DUPZdeg(THIS->g1)-1);
  gcd = DUPZexgcd(g_cofac, h_cofac, THIS->g1, THIS->h1);
  if (DUPZdeg(gcd) != 0) JERROR(JERROR_HENSEL);
//printf("cofactor is ");DUPZprint(g_cofac);

#ifdef FACTOR_DEBUG
  printf("Entering bivariate lifting loop\nLevel=");
#endif  
  corr = DUPZnew(2*dx); /* too big really */
  /* tmp and tmp2 are used as workspace */
  tmp = DUPZnew(dx);
  tmp2 = DUPZnew(dx);
  for (level=1; level <= max_height; level++)
  {
#ifdef FACTOR_DEBUG
printf("%d..",level);fflush(stdout);
#endif
    tmp->deg = -1;
    for (j=1; j < level; j++)
    {
      DUPZmul3(tmp2, glift[j], hlift[level-j]);
      DUPZadd3(tmp, tmp, tmp2);
    }
    error = DUPZsub(f2[level], tmp);
//printf("Error is ");DUPZprint(error);
    DUPZmul3(corr, error, h_cofac);
    DUPZmul2z(corr, DUPZlc(THIS->g1));
//printf("Correction is ");DUPZprint(corr);
//printf("g1 is ");DUPZprint(THIS->g1);
    glift[level] = DUPZnew(DUPZdeg(glift[0]));
    DUPZdiv4(tmp, glift[level], corr, THIS->g1); /* tmp is discarded */
//    glift[level] = DUPZrem(corr, THIS->g1);
//printf("corr to g is ");DUPZprint(glift[level]);    
//    DUPZcontent(content, glift[level]);
    DUPZdiv2z(glift[level], gcd->coeffs[0]);/*CHECK IF EXACT!!!*/
    DUPZmul3(corr, error, g_cofac);
    DUPZmul2z(corr, DUPZlc(THIS->h1));
//printf("Correction is ");DUPZprint(corr);
//printf("h1 is ");DUPZprint(THIS->h1);
    hlift[level] = DUPZnew(DUPZdeg(hlift[0]));
    DUPZdiv4(tmp, hlift[level], corr, THIS->h1); /* tmp is discarded */
//    hlift[level] = DUPZrem(corr, THIS->h1);
//printf("corr to h is ");DUPZprint(hlift[level]);
    DUPZdiv2z(hlift[level], gcd->coeffs[0]);/*CHECK IF EXACT!!!*/
    DUPZfree(error);
//    glift[level] += lcg->coeffs[level]*glift[0];
    if (level <= DUPZdeg(lcg))
    {
      DUPZshift_add(glift[level], glift[0], 0, lcg->coeffs[level]);
    }
    DUPZdiv2z(glift[level], DUPZlc(THIS->g1));
    if (level <= DUPZdeg(lch))
    {
      DUPZshift_add(hlift[level], hlift[0], 0, lch->coeffs[level]);
    }
    DUPZdiv2z(hlift[level], DUPZlc(THIS->h1));
//printf("Finally:\n");
//printf("glift[%d]=",level);DUPZprint(glift[level]);
//printf("hlift[%d]=",level);DUPZprint(hlift[level]);
  }
#ifdef FACTOR_DEBUG
  printf("\n");
#endif
  
  /* At this point at least one factor is fully lifted */
  DUPZfree(lcg); DUPZfree(lch);
  DUPZfree(g_cofac); DUPZfree(h_cofac);
  DUPZfree(corr);
  DUPZfree(gcd); DUPZfree(tmp); DUPZfree(tmp2);
  for (j=0; j <= dy; j++) DUPZfree(f2[j]);
  FREE(f2);

  DMPZdtor(THIS->g2);/*get rid of any previous value*/
  THIS->g2 = NULL;
  {
    /* In this block we negate alpha to undo the map y |-> y-alpha */
    mpz_neg(alpha, alpha);
    DMPZbivariate_padic_shift(glift, level-1, alpha);
    DMPZbivariate_padic_shift(hlift, level-1, alpha);
    mpz_clear(alpha);
  }
  DMPZdtor(THIS->g2);
  DMPZdtor(THIS->h2);

THIS->g2 = DMPZbivariate_to_DMPZ(glift, level-1, THIS->x, THIS->y);
THIS->h2 = DMPZbivariate_to_DMPZ(hlift, level-1, THIS->x, THIS->y);
#if 0
  {
    /* Get rid of spurious content due to forcing the lc */
    DMPZ content, tmp;
    tmp = DMPZbivariate_to_DMPZ(glift, level-1, THIS->x, THIS->y);    
    content = DMPZcontent_var(tmp, THIS->x);
    THIS->g2 = DMPZdiv_exact(tmp, content);
    DMPZdtor(content); DMPZdtor(tmp);
    tmp = DMPZbivariate_to_DMPZ(hlift, level-1, THIS->x, THIS->y);    
    content = DMPZcontent_var(tmp, THIS->x);
    THIS->h2 = DMPZdiv_exact(tmp, content);
    DMPZdtor(content); DMPZdtor(tmp);
  }
#endif  
  for (j=0; j < level; j++) DUPZfree(glift[j]);
  FREE(glift);
  for (j=0; j < level; j++) DUPZfree(hlift[j]);
  FREE(hlift);
  return 1;
}



/***************************************************************************/
/* Pick a (new) random substitution during the Zippel sparse interpolation */
/* of DMPZlift.  Must avoid substitutions which kill the leading term wrt  */
/* x or y.                                                                 */

static void pick_random_substitution(DMPZlifter THIS)
{
  int j;
  mpz_t tmp;
  DMPZ lcx, lcy;

  lcx = DMPZlc(THIS->f, THIS->x);
  lcy = DMPZlc(THIS->f, THIS->y);
  mpz_init(tmp);
  while(1)
  {
    for (j=0; j < NVARS; j++)
      if (THIS->lifted[j]) THIS->substitution[j] = (rand()%199); /* range is -99..99 */
    DMPZeval(tmp, lcx, THIS->substitution);
    if (mpz_sgn(tmp) == 0) continue;
    DMPZeval(tmp, lcy, THIS->substitution);
    if (mpz_sgn(tmp) == 0) continue;
    break;
  }
  mpz_clear(tmp);
  DMPZdtor(lcx);
  DMPZdtor(lcy);
}


/***************************************************************************/
#if FACTOR_DEBUG
static void Qmat_print(Qmat M)
{
  int i, j;

  for (i=0; i < M->nrows; i++)
  {
    for (j=0; j < M->ncols; j++)
    {
      mpz_out_str(stdout,10, mpq_numref(M->entry[i][j]));
      printf("  ");
    }
    printf("\n");
  }
  printf("\n");
}
#endif


/***************************************************************************/

static void DMPZlift_check_initial_subst(DMPZlifter THIS)
{
  int var;
  int *df = (int*)MALLOC(NVARS*sizeof(int));
  int *flag = (int*)MALLOC(NVARS*sizeof(int));
  DMPZ lcf;

  THIS->FAILED = 0;
  DMPZdegs(df, THIS->f);
  for (var = 0; var < NVARS; var++) { flag[var] = 1; }
  lcf = DMPZlc(THIS->f, THIS->x);

  for (var = 0; var < NVARS; var++)
  {
    DMPZ f1,lcf1;
    int d, dlcf;
    if (var == THIS->x) continue;
    if (df[var] == 0) continue;
    flag[var] = 0;
    f1 = DMPZeval_partial(THIS->f, THIS->substitution, flag);
    flag[var] = 1;
    d = DMPZdeg(f1, var);
    DMPZdtor(f1);
    if (d != df[var]) break;
    dlcf = DMPZdeg(lcf, var);
    flag[var] = 0;
    lcf1 = DMPZeval_partial(lcf, THIS->substitution, flag);
    flag[var] = 1;
    d = DMPZdeg(lcf1, var);
    DMPZdtor(lcf1);
    if (d != dlcf) break;
  }
  DMPZdtor(lcf);
  FREE(df);
  FREE(flag);
  if (var != NVARS) THIS->FAILED = 1;
}

/***************************************************************************/
/* At least 2 factors, and at least 2 variables.                           */

void DMPZlift(DMPZlifter THIS)
{
  DMPZ iter, newg;
  const int x = THIS->x;
  int i, j, k, n, y, dy, gx, gy, rank;
  Qmat M, rhs, soln;

if(THIS->FAILED<0)DMPZlift_check_initial_subst(THIS);
if(THIS->FAILED)return;
#ifdef FACTOR_DEBUG
  printf("Starting to lift\n");
  printf("DMPZlift: Initial factors are g=\n");DMPZprint(THIS->g);
  printf("DMPZlift: ................... h=\n");DMPZprint(THIS->h);
  printf("DMPZlift: Initial subst is:\n");for(i=0;i<NVARS;i++)printf("%d ",THIS->substitution[i]);printf("\n");
#endif
/*
 *  while (p = NextPrime(PS))
 *  {
 *    for (iter=F; iter; iter=iter->next)
 *      if (mpz_fdiv_ui(iter->coeff, p) == 0) break;
 *    if (iter == NULL) break;
 *  }
 *  if (p == 0) { printf("DMPZlift: ran out of primes.\n"); exit(1); }
 *  Fp = DMPZ_to_DMPFF(F, p);
 */

  /* Pick which variable to lift next (just take the first we find...) */
  for (y=0; y < NVARS; y++)
    if (THIS->degs[y] > 0 && THIS->lifted[y] == 0) break;
  if (y == NVARS) /* no more variables to lift */
  {
    /* Before returning remove any excess content arising from forcing the lc */
    DMPZ c, tmp;
//printf("Removing excess content from ");DMPZprint(THIS->g);
    c = DMPZcontent_var(THIS->g, THIS->x);
//printf("content = ");DMPZprint(c);
    tmp = DMPZdiv_exact(THIS->g, c);
    DMPZdtor(c);
    DMPZdtor(THIS->g); THIS->g = tmp;
//printf("Removed excess content: ");DMPZprint(THIS->g);
    c = DMPZcontent_var(THIS->h, THIS->x);
    tmp = DMPZdiv_exact(THIS->h, c);
    DMPZdtor(c);
    DMPZdtor(THIS->h); THIS->h = tmp;
//printf("Lifting completed: all variables lifted.\n");
//printf("Resulting factor is:\n");
//DMPZprint(THIS->g);
    return;
  }
  THIS->y = y;
#ifdef FACTOR_DEBUG
printf("We shall be lifting the variable x[%d]\n", y);
#endif

  DMPZdtor(THIS->g_skeleton);
  THIS->g_skeleton = make_skeleton(THIS->g, x);
#ifdef FACTOR_DEBUG
printf("g_skeleton is "); DMPZprint(THIS->g_skeleton);
#endif
  DMPZdtor(THIS->h_skeleton);
  THIS->h_skeleton = make_skeleton(THIS->h, x);
#ifdef FACTOR_DEBUG
printf("h_skeleton is "); DMPZprint(THIS->h_skeleton);
#endif
/*
 *  {
 *    int leng = DMPZ_nterms(THIS->g_skeleton);
 *    int lenh = DMPZ_nterms(THIS->h_skeleton);
 *    if (leng > lenh) n = leng; else n = lenh;
 *  }
 */
  n = DMPZ_nterms(THIS->g_skeleton); /* Should be smarter about this */
  dy = DMPZdeg(THIS->f, y);
  gx = DMPZdeg(THIS->g, THIS->x);
  gy = 2*dy; /* upper bound for deg in y of g*lc(f) */
#ifdef FACTOR_DEBUG
printf("Height to lift to is %d\n", gy);
printf("Size of linear system is %dx%d\n", n, n);
printf("RHS will have %d*%d=%d vectors.\n", 1+gx, 1+gy, (1+gx)*(1+gy));
#endif
  M = Qmat_ctor(n, n);
  rhs = Qmat_ctor(n, (1+gx)*(1+gy));
  soln = Qmat_ctor(n, (1+gx)*(1+gy));
  for (i=0; i < n; i++)
  {
#ifdef FACTOR_DEBUG
    printf("%d/%d..",i,n);fflush(stdout);
#endif
  pick_new_rnd_subst:
    /* Find a new substitution which preserves degree in x and y */
    pick_random_substitution(THIS);
#ifdef FACTOR_DEBUG
//printf("Creating new univariate images...\n");
//printf("Univariate variable is x[%d].\n", THIS->x);
printf("Substitution is: ");for (j=0; j < NVARS; j++)printf("%3d ",THIS->substitution[j]);printf("\n");
#endif
    DUPZfree(THIS->g1);
    DUPZfree(THIS->h1);
    THIS->g1 = DMPZmap_to_univariate(THIS->g, THIS->x, THIS->substitution);
    THIS->h1 = DMPZmap_to_univariate(THIS->h, THIS->x, THIS->substitution);
    {
      DUPZ gcd = DUPZgcd(THIS->g1, THIS->h1);
      int d;
      d = DUPZdeg(gcd);
      DUPZfree(gcd);
      if (d > 0) goto pick_new_rnd_subst;
    }
    {
      int OK = DMPZlift_to_bivariate(THIS);
      if (!OK) goto pick_new_rnd_subst;
    }

#ifdef FACTOR_DEBUG
printf("Bivariate lifted is ");DMPZprint(THIS->g2);
#endif
    /* Fill row of M using skel */
    iter = THIS->g_skeleton;
    for (j=0; j < n; j++)
    {
//printf("Filling row of matrix with values of PP: ");
//{int j;for (j=0;j<NVARS;j++)printf("%3d ", iter->exps[j]);printf("\n");}
      DMPZeval_power_product(mpq_numref(M->entry[i][j]), iter->exps, THIS->substitution);
      iter = iter->next;
    }
    /* Fill row of rhs using lifted bivariate factor */
    for (iter = THIS->g2; iter; iter = iter->next)
    {
      int dx, dy;
      dx = iter->exps[x];
      dy = iter->exps[y];
//printf("Found a term with exponent (%d, %d) <--> %d\n", dx, dy, dx*(1+gy)+dy);
      mpz_add(mpq_numref(rhs->entry[i][dx*(1+gy)+dy]),
              mpq_numref(rhs->entry[i][dx*(1+gy)+dy]),
              iter->coeff);
    }
  }
#ifdef FACTOR_DEBUG
printf("\nSolving linear system.\n");
//Qmat_print(M);
//Qmat_print(rhs);
#endif
  rank = Qsolve(soln, M, rhs);
  if (rank != n)
  {
#ifdef FACTOR_DEBUG
    printf("UNLUCKY, system was singular.\n");
#endif
    Qmat_dtor(soln);
    Qmat_dtor(rhs);
    Qmat_dtor(M);
          
    DMPZlift(THIS);
    return;
  }
#ifdef FACTOR_DEBUG
  printf("solved linear system\n");
#endif
  /* --> new skeleton */
  THIS->lifted[y] = 1;
  newg = NULL;
  i = 0;
  for (iter = THIS->g_skeleton; iter; iter = iter->next)
  {
    for (j=0; j < (1+gx)*(1+gy); j++)
      if (mpz_sgn(mpq_numref(soln->entry[i][j])))
      {
        int *exps = (int*)MALLOC(NVARS*sizeof(int));
if (mpz_cmp_ui(mpq_denref(soln->entry[i][j]), 1) != 0) THIS->FAILED = 1;
        for (k=0; k < NVARS; k++) exps[k] = iter->exps[k];
//printf("i=%d\n",j);
        exps[x] = j/(1+gy);
        exps[y] = j%(1+gy);
//printf("extracting a term with exponent (%d, %d) <--> %d\n", exps[x], exps[y], j);
        newg = DMPZprepend(mpq_numref(soln->entry[i][j]), exps, newg);
      }
    i++;
  }
  Qmat_dtor(soln);
  Qmat_dtor(rhs);
  Qmat_dtor(M);
  DMPZdtor(THIS->g); THIS->g = newg;
#ifdef FACTOR_DEBUG
  if (THIS->FAILED) printf("Linear system was bad because it had a rational solution...\n");
#endif  
  if (THIS->FAILED) return;
  THIS->lifted[y] = 1;
#ifdef FACTOR_DEBUG
printf("DMPZlift: new g is "); DMPZprint(newg);
printf("DMPZlift: Finding h by division...\n");
#endif
  {
    DMPZ newh, lcnewh, lch, scale;
    int *flag;
    flag = (int*)MALLOC(NVARS*sizeof(int));
    for (i=0; i < NVARS; i++) flag[i] = (THIS->lifted[i] == 0);
    {
      DMPZ flcf = DMPZeval_partial(THIS->flcf, THIS->substitution, flag);
      newh = DMPZdiv_exact(flcf, THIS->g);
      DMPZdtor(flcf);
    }
    if (newh == NULL) { THIS->FAILED = 1; goto tidy_up; }
    {
      DMPZ contentx = DMPZcontent_var(newh, x);
      DMPZ tmp = DMPZdiv_exact(newh, contentx);
      DMPZdtor(newh);
      DMPZdtor(contentx);
      newh = tmp;
    }
    lcnewh = DMPZlc(newh, x);
    lch = DMPZeval_partial(THIS->lch, THIS->substitution, flag);
//printf("lch=");DMPZprint(THIS->lch);
//printf("lcnewh=");DMPZprint(lcnewh);
    scale = DMPZdiv_exact(lch, lcnewh);
//printf("scale=");DMPZprint(scale);
    DMPZdtor(THIS->h); THIS->h = DMPZmul(newh, scale);
//printf("h computed.\n");
    DMPZdtor(lcnewh);
    DMPZdtor(lch);
    DMPZdtor(scale);
  tidy_up:
    DMPZdtor(newh);
    FREE(flag);
  }
#ifdef FACTOR_DEBUG
printf("DMPZlift: h found by division is ");DMPZprint(THIS->h);
#endif
if (THIS->FAILED) return;
//printf("newg = ");DMPZprint(THIS->g);  
//printf("newh = ");DMPZprint(THIS->h);  
  DMPZlift(THIS);
}
