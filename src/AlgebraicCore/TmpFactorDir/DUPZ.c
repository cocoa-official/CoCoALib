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

/* need stddef for NULL */
#include <stddef.h>
#include "jalloc.h"
#include "DUPZ.h"
#include "jaaerror.h"

int DUPZdeg(const DUPZ f)
{
  return f->deg;
}


DUPZ DUPZnew(int maxdeg)
{
  DUPZ ans = (DUPZ)MALLOC(sizeof(struct DUPZstruct));
  int i;

  ans->coeffs = NULL;
  if (maxdeg >= 0) ans->coeffs = (mpz_t*)MALLOC((maxdeg+1)*sizeof(mpz_t));
  for (i=0; i<=maxdeg; i++) mpz_init(ans->coeffs[i]);
  ans->maxdeg = maxdeg;
  ans->deg = -1;
  return ans;
}

DUPZ int_to_DUPZ(int i)
{
  DUPZ ans = DUPZnew(0);

  mpz_set_si(ans->coeffs[0], i);
  ans->deg = 0;
  return ans;
}

void DUPZinit(DUPZ x, int maxdeg)
{
  int i;

  x->coeffs = NULL;
  if (maxdeg >= 0) x->coeffs = (mpz_t*)MALLOC((maxdeg+1)*sizeof(mpz_t));
  for (i=0; i<=maxdeg; i++) mpz_init(x->coeffs[i]);
  x->maxdeg = maxdeg;
  x->deg = -1;
}

void DUPZfree(DUPZ x)
{
  int i;

  for (i=0; i <= x->maxdeg; i++) mpz_clear(x->coeffs[i]);
  if (x->maxdeg >= 0) FREE(x->coeffs);
  FREE(x);
}

void DUPZswap(DUPZ x, DUPZ y)
{
  mpz_t *swap1;
  int swap2;

  swap1 = x->coeffs; x->coeffs = y->coeffs; y->coeffs = swap1;
  swap2 = x->deg;    x->deg    = y->deg;    y->deg    = swap2;
  swap2 = x->maxdeg; x->maxdeg = y->maxdeg; y->maxdeg = swap2;
}


DUPZ DUPZcopy(const DUPZ x)
{
  int dx, i;
  DUPZ ans = (DUPZ)MALLOC(sizeof(struct DUPZstruct));

  dx = x->deg;
  ans->deg = dx;
  ans->maxdeg = dx;
  if (dx < 0) { ans->coeffs = NULL; return ans; }
  ans->coeffs = (mpz_t*)MALLOC((dx+1)*sizeof(mpz_t));
  for (i=0; i<=dx; i++) mpz_init_set(ans->coeffs[i], x->coeffs[i]);
  return ans;
}


void DUPZcopy2(DUPZ dest, const DUPZ src)
{
  int i;

  if (dest == src) return; /* copying onto itself */
  if (src->deg > dest->maxdeg)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    dest->deg = -1;
    return;
  }
  dest->deg = src->deg;
  for (i=src->deg; i>=0; i--)
    mpz_set(dest->coeffs[i], src->coeffs[i]);
}


int DUPZequal(const DUPZ x, const DUPZ y)
{
  int dx, dy;

  dx = x->deg;
  dy = y->deg;
  if (dx < dy) return DUPZequal(y, x);
  for (; dx > dy; dx--) if (mpz_sgn(x->coeffs[dx]) != 0) return 0;
  for (; dx >= 0; dx--) if (mpz_cmp(x->coeffs[dx], y->coeffs[dx])) return 0;
  return 1;
}


/***************************************************************************/

void DUPZshift(DUPZ f, int xpower)
{
  int i, df;

  df = DUPZdeg(f);
  if (df < 0 || xpower == 0) return;
  if (df + xpower < 0) { f->deg = -1; return; }
  if (df + xpower > f->maxdeg) { JERROR(JERROR_DEG_TOO_LOW); return; }
  f->deg = df+xpower;
  if (xpower < 0)
  {
    for (i=-xpower; i <= df; i++) mpz_set(f->coeffs[i+xpower], f->coeffs[i]);
    return;
  }
  for (i=df; i >= 0; i--) mpz_set(f->coeffs[i+xpower], f->coeffs[i]);
  for (i=xpower-1; i >= 0; i--) mpz_set_ui(f->coeffs[i], 0);
}

void DUPZreverse(DUPZ f)
{
  int i, df;
  mpz_t swap;
  
  if (mpz_sgn(f->coeffs[0]) == 0) { JERROR(JERROR_DUPZ_REVERSE); return; }
  df = DUPZdeg(f);
  mpz_init(swap);
  for (i=(df-1)/2; i >= 0; i--)
  {
    mpz_set(swap, f->coeffs[i]);
    mpz_set(f->coeffs[i], f->coeffs[df-i]);
    mpz_set(f->coeffs[df-i], swap);
  }
  mpz_clear(swap);
}

/***************************************************************************/

void DUPZneg1(DUPZ f)
{
  int i;
  for (i=DUPZdeg(f); i >= 0; i--)
    mpz_neg(f->coeffs[i], f->coeffs[i]);
}

/***************************************************************************/


void DUPZadd3(DUPZ sum, const DUPZ x, const DUPZ y)
{
  mpz_t tmp;
  int d, dx, dy;
  
  dx = DUPZdeg(x);
  dy = DUPZdeg(y);
  if (dx == dy)
  {
    mpz_init(tmp);
    for(d = dx; d >= 0; d--)
    {
      mpz_add(tmp, x->coeffs[d], y->coeffs[d]);
      if (mpz_sgn(tmp) != 0) break;
    }
    if (sum->maxdeg < d) { JERROR(JERROR_DEG_TOO_LOW); return; }
    if (d >= 0) mpz_set(sum->coeffs[d], tmp);
    mpz_clear(tmp);
    sum->deg = d;
    d--;
  }
  else
  {
    d = (dx>dy)?dx:dy;
    if (sum->maxdeg < d) { JERROR(JERROR_DEG_TOO_LOW); return; }
    sum->deg = d;
    for (; d>dx; d--) mpz_set(sum->coeffs[d], y->coeffs[d]);
    for (; d>dy; d--) mpz_set(sum->coeffs[d], x->coeffs[d]);
  }
  for (; d >= 0; d--)
    mpz_add(sum->coeffs[d], x->coeffs[d], y->coeffs[d]);
}


DUPZ DUPZadd(const DUPZ x, const DUPZ y)
{
  DUPZ ans;
  int d, dx, dy;

  dx = DUPZdeg(x);
  dy = DUPZdeg(y);
  d = (dx>dy)?dx:dy;
  ans = DUPZnew(d);
  ans->deg = d;
  for (; d>dx; d--) mpz_set(ans->coeffs[d], y->coeffs[d]);
  for (; d>dy; d--) mpz_set(ans->coeffs[d], x->coeffs[d]);
  for (; d>=0; d--)
    mpz_add(ans->coeffs[d], x->coeffs[d], y->coeffs[d]);
  for (d = ans->deg; d >= 0 && mpz_sgn(ans->coeffs[d]) == 0; d--) {}
  ans->deg = d;
  return ans;
}

/***************************************************************************/

void DUPZsub3(DUPZ diff, const DUPZ x, const DUPZ y)
{
  mpz_t tmp;
  int d, dx, dy;
  
  dx = DUPZdeg(x);
  dy = DUPZdeg(y);
  if (dx == dy)
  {
    mpz_init(tmp);
    for(d = dx; d >= 0; d--)
    {
      mpz_sub(tmp, x->coeffs[d], y->coeffs[d]);
      if (mpz_sgn(tmp) != 0) break;
    }
    if (diff->maxdeg < d) { JERROR(JERROR_DEG_TOO_LOW); return; }
    if (d >= 0) mpz_set(diff->coeffs[d], tmp);
    mpz_clear(tmp);
    diff->deg = d;
    d--;
  }
  else
  {
    d = (dx>dy)?dx:dy;
    if (diff->maxdeg < d) { JERROR(JERROR_DEG_TOO_LOW); return; }
    diff->deg = d;
    for (; d>dx; d--) mpz_neg(diff->coeffs[d], y->coeffs[d]);
    for (; d>dy; d--) mpz_set(diff->coeffs[d], x->coeffs[d]);
  }
  for (; d >= 0; d--)
    mpz_sub(diff->coeffs[d], x->coeffs[d], y->coeffs[d]);
}


DUPZ DUPZsub(const DUPZ x, const DUPZ y)
{
  DUPZ ans;
  int d, dx, dy;

  dx = DUPZdeg(x);
  dy = DUPZdeg(y);
  d = (dx>dy)?dx:dy;
  ans = DUPZnew(d);
  ans->deg = d;
  for (; d>dx; d--) mpz_neg(ans->coeffs[d], y->coeffs[d]);
  for (; d>dy; d--) mpz_set(ans->coeffs[d], x->coeffs[d]);
  for (; d>=0; d--)
    mpz_sub(ans->coeffs[d], x->coeffs[d], y->coeffs[d]);
  for (d = ans->deg; d >= 0 && mpz_sgn(ans->coeffs[d]) == 0; d--) {}
  ans->deg = d;
  return ans;
}

/***************************************************************************/

void DUPZmul2z(DUPZ f, const mpz_t c) /* c must be non-zero */
{
  int i;

  if (mpz_cmp_ui(c, 1) == 0) return;
  for (i=DUPZdeg(f); i >= 0; i--)
    mpz_mul(f->coeffs[i], f->coeffs[i], c);
}


DUPZ DUPZmul(const DUPZ x, const DUPZ y)
{
  int dx, dy, dans;
  DUPZ ans;
  dx = DUPZdeg(x);
  dy = DUPZdeg(y);
  if (dx < 0 || dy < 0) return DUPZnew(-1); /* zero times anything is zero */
  dans = dx+dy;
  ans = DUPZnew(dans);
  DUPZmul3(ans, x, y);
  return ans;
}

/* PLACEHOLDER TRIVIAL IMPLEMENTATION */
void DUPZsquare(DUPZ f)
{
  DUPZ fcopy;
  
  if (DUPZdeg(f) < 0) return;
  fcopy = DUPZcopy(f);
  DUPZmul3(f, fcopy, f);
  DUPZfree(fcopy);
}

void DUPZmul3(DUPZ ans, const DUPZ x, const DUPZ y)
{
  int i, dx, dy, dans, j;
  mpz_t *ac, *xc, *yc, tmp, xci;

  if (ans == y)
  {
    if (ans == x) DUPZsquare(ans);
    else DUPZmul3(ans, y, x);
    return;
  }
  /* we are now sure that ans and y are different */
  dx = DUPZdeg(x);
  dy = DUPZdeg(y);
  if (dx < 0 || dy < 0)
  {
    ans->deg = -1; /* zero times anything is zero */
    return;
  }
  dans = dx+dy;
  if (ans->maxdeg < dans)  /* error case */
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }

  dans = dx+dy;
  ans->deg = dans;
  ac = ans->coeffs;
  xc = x->coeffs;
  yc = y->coeffs;
  mpz_init(tmp);
  mpz_init(xci);

  for (i=1; i<=dy; i++) mpz_set_ui(ac[dx+i], 0);
  for (i=dx; i >= 0; i--)
  {
    mpz_set(xci, xc[i]);
    mpz_set_ui(ac[i], 0);  /*  NB xc[i] and ac[i] might be identical */
    if (mpz_sgn(xci) == 0) continue;
    for (j=dy; j >= 0; j--)
    {
      mpz_mul(tmp, xci, yc[j]);
      mpz_add(ac[i+j], ac[i+j], tmp);
    }
  }

  mpz_clear(xci);
  mpz_clear(tmp);
  return;
}


/***************************************************************************/

void DUPZmod2z(DUPZ f, const mpz_t c) /* c must be non-zero */
{
  int i;

  for (i=DUPZdeg(f); i >= 0; i--)
    mpz_fdiv_r(f->coeffs[i], f->coeffs[i], c);
  for (i=DUPZdeg(f); i >= 0; i--) if (mpz_sgn(f->coeffs[i]) != 0) break;
  f->deg = i;
}

/***************************************************************************/

void DUPZshift_sub(DUPZ f, const DUPZ g, const int deg, const mpz_t coeff)
{
  mpz_t negated;
  mpz_init_set(negated, coeff);
  mpz_neg(negated, negated);
  DUPZshift_add(f, g, deg, negated);
  mpz_clear(negated);
}

void DUPZshift_add(DUPZ f, const DUPZ g, const int deg, const mpz_t coeff)
{
  int i, dg, df, dans;
  mpz_t tmp;

  df = f->deg;
  dg = g->deg;
  if (dg < 0) return;
  if (mpz_sgn(coeff) == 0) return;
  dans = (deg+dg > df)? deg+dg : df;
  if (dans > f->maxdeg) return; /* ERROR */
  for (i=df+1; i <= dans; i++) mpz_set_ui(f->coeffs[i], 0);

  mpz_init(tmp);
  /* THE ORDER OF ITERATION in this loop is important (see e.g. DUPZmod2) */
  for (i=0; i <= dg; i++)
  {
    mpz_mul(tmp, coeff, g->coeffs[i]);
    mpz_add(f->coeffs[i+deg], f->coeffs[i+deg], tmp);
  }

  mpz_clear(tmp);

  for (; (dans >= 0) && (mpz_sgn(f->coeffs[dans]) == 0); dans--){}
  f->deg = dans;
  return;
}

/***************************************************************************/

void DUPZdiv2z(DUPZ f, mpz_t n)
{
  int i;

  if (mpz_sgn(n) == 0)
  {
      JERROR(JERROR_DIV_BY_ZERO);
      return;
  }
  if (mpz_cmp_ui(n, 1) == 0) return;
  for (i=0; i <= DUPZdeg(f); i++)
    mpz_fdiv_q(f->coeffs[i], f->coeffs[i], n);
}


DUPZ DUPZdiv(const DUPZ num, const DUPZ den)
{
  int dnum, dden;
  DUPZ quot, rem;
  dnum = num->deg;
  dden = den->deg;
  if (dden < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return DUPZnew(-1);
  }
  if (dnum < 0) return DUPZnew(-1);
  quot = DUPZnew(dnum-dden);
  rem = DUPZnew(dden-1);
  DUPZdiv4(quot, rem, num, den);
  DUPZfree(rem);
  return quot;
}


DUPZ DUPZrem(const DUPZ num, const DUPZ den)
{
  int dnum, dden;
  DUPZ quot, rem;

  dnum = num->deg;
  dden = den->deg;
  if (dden < 0)
    {
      JERROR(JERROR_DIV_BY_ZERO);
      return DUPZnew(-1);
    }
  if (dnum < 0) return DUPZnew(-1);
  quot = DUPZnew(dnum-dden);
  rem = DUPZnew(dden-1);
  DUPZdiv4(quot, rem, num, den);
  DUPZfree(quot);
  return rem;
}


void DUPZrem2(DUPZ x, const DUPZ m)
{
  int dx, dm;
  mpz_t q, r, *xc;

  dm = DUPZdeg(m);
  dx = DUPZdeg(x);
  if (dm < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return;
  }

  mpz_init(q); mpz_init(r);
  xc = x->coeffs;
  while (dx >= dm)
  {
    mpz_tdiv_qr(q, r, xc[dx], DUPZlc(m));
    DUPZshift_sub(x, m, (dx-dm), q);
    if (mpz_sgn(r) != 0) break;
    dx = DUPZdeg(x);
  }
  mpz_clear(q); mpz_clear(r);
}


void DUPZdiv4(DUPZ quot, DUPZ rem, const DUPZ num, const DUPZ den)
{
  int dden, dnum, dquot, i;
  mpz_t *q, r;
  DUPZ tmp;

  if (quot == rem || quot == den || rem == den)
  {
    JERROR(JERROR_DIV4_ARGS);
    return;
  }
  dden = DUPZdeg(den);
  if (dden < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return;
  }

  dnum = DUPZdeg(num);
  if (dnum < 0)  /* quotient and remainder are 0 */
  {
    quot->deg = -1;
    rem->deg = -1;
    return; 
  }


  if (dnum < dden) /* quotient is zero, remainder = num */
  {
    if (rem->maxdeg < dnum)
    {
      JERROR(JERROR_DEG_TOO_LOW);
      return;
    }
    DUPZcopy2(rem, num);
    quot->deg = -1;
    return;
  }
  dquot = dnum-dden;
  if (quot->maxdeg < dquot)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    return;
  }

  /* we work with a copy of num; do this before zeroing quot! */
  tmp = DUPZcopy(num);

  /* get quot ready for the answer */
  quot->deg = dquot;
  q = quot->coeffs;
  /* fill quot with 0 as the loop below doesn't assign zero coeffs */
  for (i=dquot; i>=0; i--) mpz_set_ui(q[i], 0);

  mpz_init(r); /* used to check that divisions are exact */

  while (dquot >= 0)
  {
    mpz_fdiv_qr(q[dquot], r, DUPZlc(tmp), DUPZlc(den));
    DUPZshift_sub(tmp, den, dquot, q[dquot]);
    if (mpz_sgn(r) != 0) break;
    dquot = DUPZdeg(tmp) - DUPZdeg(den);
  }

  DUPZcopy2(rem, tmp);
  DUPZfree(tmp);
  mpz_clear(r);
}

/***************************************************************************/

/* g MUST BE MONIC!!! */
void DUPZmod2(DUPZ f, DUPZ g, mpz_t Q)
{
  int df, dg;

  df = DUPZdeg(f);
  dg = DUPZdeg(g);
  while (1)
  {
    for (; df >= 0; df--)
    {
      mpz_fdiv_r(f->coeffs[df], f->coeffs[df], Q);
      if (mpz_sgn(f->coeffs[df]) != 0) break;
    }
    f->deg = df;
    if (df < dg) break;
    DUPZshift_sub(f, g, df-dg, f->coeffs[df]);
  }
  /* Reduce remaining coeffs of result mod Q */
  for (; df >= 0; df--)
    mpz_fdiv_r(f->coeffs[df], f->coeffs[df], Q);
}



void DUPZmdiv2(DUPZ f, DUPZ g, mpz_t Q)
{
  int df, dg, df_orig, i;

  df = DUPZdeg(f);
  if (df < 0) return;
  df_orig = df;
  dg = DUPZdeg(g);
  g->deg = dg-1;  /* DANGEROUS HACK */
  while (1)
  {
    for (; df >= 0; df--)
    {
      mpz_fdiv_r(f->coeffs[df], f->coeffs[df], Q);
      if (mpz_sgn(f->coeffs[df]) != 0) break;
    }
    f->deg = df;
    if (df < dg) break;
    DUPZshift_sub(f, g, df-dg, f->coeffs[df]);
    --df;
  }
  g->deg = dg;    /* UNHACK G */
  df = df_orig;
  f->deg = -1;
  for (i = 0; i <= df-dg; i++)
  {
    mpz_set(f->coeffs[i], f->coeffs[i+dg]);
    if (mpz_sgn(f->coeffs[i])) f->deg = i;
  }
}

/***************************************************************************/

void DUPZmmul3z(DUPZ product, DUPZ f, mpz_t c, mpz_t Q)
{
  int i;

  if (DUPZdeg(f) > product->maxdeg) JERROR(JERROR_DEG_TOO_LOW);
  for (i=0; i <= DUPZdeg(f); i++)
  {
    mpz_fdiv_r(product->coeffs[i], f->coeffs[i], Q);
    mpz_mul(product->coeffs[i], product->coeffs[i], c);
    mpz_fdiv_r(product->coeffs[i], product->coeffs[i], Q);
    if (mpz_sgn(product->coeffs[i])) product->deg = i;
  }
}

/***************************************************************************/
/* On calling we have poly = sum f[i]*x^i where f[i] is f->coeffs[i].      */
/* On return we have  poly = sum f[i]*(x-a)^i, i.e. a linear shift in "x". */
/* Thus in particular on output f[0] = poly evaluated at a.                */

void DUPZlinear_shift(DUPZ f, mpz_t a)
{
  int i, j, df;
  mpz_t tmp;

  mpz_init(tmp);
  df = DUPZdeg(f);
  for (i=df-1; i >= 0; i--)
    for (j=i; j < df; j++)
    {
      mpz_mul(tmp, a, f->coeffs[j+1]);
      mpz_add(f->coeffs[j], f->coeffs[j], tmp);
    }
  mpz_clear(tmp);
}
