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
#include "DUPI.h"
#include "jaaerror.h"

int DUPIdeg(const DUPI f)
{
  return f->deg;
}


DUPI DUPInew(const int maxdeg)
{
  DUPI ans = (DUPI)MALLOC(sizeof(struct DUPIstruct));
  ans->coeffs = NULL;
  if (maxdeg >= 0) ans->coeffs = (int32*)MALLOC((maxdeg+1)*sizeof(int32));
  ans->maxdeg = maxdeg;
  ans->deg = -1;
  return ans;
}

void DUPIfree(DUPI x)
{
  if (x->maxdeg >= 0) FREE(x->coeffs);
  FREE(x);
}

void DUPIswap(DUPI x, DUPI y)
{
  int32 *swap1;
  int swap2;

  swap1 = x->coeffs; x->coeffs = y->coeffs; y->coeffs = swap1;
  swap2 = x->deg;    x->deg    = y->deg;    y->deg    = swap2;
  swap2 = x->maxdeg; x->maxdeg = y->maxdeg; y->maxdeg = swap2;
}


DUPI DUPIcopy(const DUPI x)
{
  int dx, i;
  DUPI ans = (DUPI)MALLOC(sizeof(struct DUPIstruct));

  dx = x->deg;
  ans->deg = dx;
  ans->maxdeg = dx;
  if (dx < 0) { ans->coeffs = NULL; return ans; }
  ans->coeffs = (int32*)MALLOC((dx+1)*sizeof(int32));
  for (i=0; i<=dx; i++) ans->coeffs[i] = x->coeffs[i];
  return ans;
}


void DUPIcopy2(DUPI dest, const DUPI src)
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
    dest->coeffs[i] = src->coeffs[i];
}


int DUPIequal(const DUPI x, const DUPI y)
{
  int dx, dy;

  dx = x->deg;
  dy = y->deg;
  if (dx < dy) return DUPIequal(y, x);
  for (; dx > dy; dx--) if (x->coeffs[dx] != 0) return 0;
  for (; dx >= 0; dx--) if (x->coeffs[dx] != y->coeffs[dx]) return 0;
  return 1;
}


/***************************************************************************/


DUPI DUPIadd(const DUPI x, const DUPI y)
{
  DUPI ans;
  int d, dx, dy;

  dx = DUPIdeg(x);
  dy = DUPIdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (x->coeffs[d] + y->coeffs[d] == 0)) d--;
  ans = DUPInew(d);
  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] = y->coeffs[d];
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = x->coeffs[d] + y->coeffs[d];
  return ans;
}

void DUPIadd3(DUPI ans, const DUPI x, const DUPI y)
{
  int d, dx, dy;

  dx = DUPIdeg(x);
  dy = DUPIdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (x->coeffs[d] + y->coeffs[d] == 0)) d--;
  if (ans->maxdeg < d)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }

  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] = y->coeffs[d];
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = x->coeffs[d] + y->coeffs[d];
  return;
}

/***************************************************************************/


DUPI DUPIsub(const DUPI x, const DUPI y)
{
  DUPI ans;
  int d, dx, dy;

  dx = DUPIdeg(x);
  dy = DUPIdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (x->coeffs[d] - y->coeffs[d] == 0)) d--;
  ans = DUPInew(d);
  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] =  - y->coeffs[d];
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = x->coeffs[d] - y->coeffs[d];
  return ans;
}

void DUPIsub3(DUPI ans, const DUPI x, const DUPI y)
{
  int d, dx, dy;

  dx = DUPIdeg(x);
  dy = DUPIdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (x->coeffs[d] - y->coeffs[d] == 0)) d--;
  if (ans->maxdeg < d)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }

  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] = - y->coeffs[d];
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = x->coeffs[d] - y->coeffs[d];
  return;
}

/***************************************************************************/

DUPI DUPImul(const DUPI x, const DUPI y)
{
  int dx, dy, dans;
  DUPI ans;
  dx = DUPIdeg(x);
  dy = DUPIdeg(y);
  if (dx < 0 || dy < 0) return DUPInew(-1); /* zero times anything is zero */
  dans = dx+dy;
  ans = DUPInew(dans);
  DUPImul3(ans, x, y);
  return ans;
}

void DUPImul3(DUPI ans, const DUPI x, const DUPI y)
{
  int i, dx, dy, dans;
  int32 *xc, *yc, *ac, *xci, *yci, *aci;

  if (ans == y)
  {
    if (ans == x) DUPIsquare(ans);
    else DUPImul3(ans, y, x);
    return;
  }
  /* we are now sure that ans and y are different */
  dx = DUPIdeg(x);
  dy = DUPIdeg(y);
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
  ac = ans->coeffs + dx;
  xc = x->coeffs;
  yc = y->coeffs;

  for (i=1; i<=dy; i++) ac[i] = 0;
  for (xci=xc+dx; xci>=xc; xci--, ac--)
  {
    if (!*xci) { *ac = 0; continue; }
    *ac = 0; /* perhaps xci == ac, so must test *xci before this line */
    for (yci=yc+dy, aci=ac+dy; yci>=yc; yci--, aci--)
    {
      *aci += *xci * *yci;
    }
  }

  ans->deg = dans;
  return;
}


/***************************************************************************/
/* This is an "in place" squaring function.                                */
/* The multiplication routine above cannot do "in place" squaring, and     */
/* while this routine can easily be modified to do general multiplication  */
/* it is noticeably slower than the routine above.                         */


void DUPIsquare(DUPI f)
{
  int df, da, *hi, *histart, *lo, *lostart, sum;
  int32 *fc;

  
  if (f->deg < 0) return;    /* trivial case */
  if (f->maxdeg < 2*f->deg)  /* error case */
  {
    JERROR(JERROR_DEG_TOO_LOW);
    return;
  }
  df = DUPIdeg(f);
  fc = f->coeffs;
  histart = fc + df;
  lostart = histart;
  for (da=2*df; da>=0; da--)
  {
    sum = 0;
    for (hi=histart, lo=lostart; hi > lo; hi--, lo++)
    {
      sum += *hi * *lo;
    }
    sum += sum;
    if (hi == lo) sum += *hi * *hi;
    fc[da] = sum;
    if (da > df) lostart--; else histart--;
  }

  f->deg = 2*df;
  return;
}


/****************************************************************************/

DUPI DUPIexpt(const DUPI base, const int power)
{
  DUPI ans;

  if (DUPIdeg(base) < 0 || power < 0)
  {
    if (power == 0) JERROR(JERROR_ZERO_TO_POWER_ZERO);
    return DUPInew(-1);
  }
  if (power == 0)
  {
    ans = DUPInew(0);
    ans->deg = 0;
    ans->coeffs[0] = 1;
    return ans;
  }
  if (power == 1) return DUPIcopy(base);
  ans = DUPInew(power*DUPIdeg(base));
  DUPIexpt3(ans, base, power);
  return ans;
}


void DUPIexpt3aux(DUPI ans, DUPI base, int n);

void DUPIexpt3(DUPI ans, const DUPI base, const int power)
{
  DUPI f;

  if (DUPIdeg(base) < 0 || power < 0)
  {
    if (power == 0) JERROR(JERROR_ZERO_TO_POWER_ZERO);
    ans->deg = -1;
    return;
  }
  if (ans->maxdeg < power*DUPIdeg(base))
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }
  if (power == 0)
  {
    ans->deg = 0;
    ans->coeffs[0] = 1;
    return;
  }
  if (power == 1)
  {
    DUPIcopy2(ans, base);
    return;
  }

  /* The next two lines are in case ans and base are identical */
  if (ans == base) f = DUPIcopy(base);
  else f = base;


  DUPIexpt3aux(ans, f, power);
  if (f != base) DUPIfree(f);
}

void DUPIexpt3aux(DUPI ans, DUPI base, int n)
{
  if (n == 1)
  {
    DUPIcopy2(ans, base);
    return;
  }
  DUPIexpt3aux(ans, base, n/2);
  DUPIsquare(ans);
  if ((n&1) == 0) return;
  DUPImul3(ans, ans, base);
}

/***************************************************************************/


void DUPIshift_add_raw(int32 *dest, const int32 *src, const int32 *srclast, const int32 scalar)
{
  if (scalar == 0) return;
  for (;src <= srclast; src++, dest++)
    *dest += scalar * *src;
}

void DUPIshift_add(DUPI f, const DUPI g, const uint32 deg, const int32 coeff)
{
  int i, dg;

  dg = g->deg;
  if (dg < 0) return;
  if (coeff == 0) return;
  if (dg+deg > f->maxdeg) return; /* ERROR */
  if (deg+dg > f->deg)
  {
    for (i=f->deg+1; i <= dg+deg; i++) f->coeffs[i] = 0;
    f->deg = dg+deg;
  }
  DUPIshift_add_raw(f->coeffs, g->coeffs, g->coeffs+dg, coeff);
  return;
}

/***************************************************************************/

void DUPIdiv2z(DUPI f, int n)
{
  int i;

  if (n == 0)
  {
      JERROR(JERROR_DIV_BY_ZERO);
      return;
  }
  for (i=0; i <= DUPIdeg(f); i++) f->coeffs[i] /= n;
}


DUPI DUPIdiv(const DUPI num, const DUPI den)
{
  int dnum, dden;
  DUPI quot, rem;
  dnum = num->deg;
  dden = den->deg;
  if (dden < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return DUPInew(-1);
  }
  if (dnum < 0) return DUPInew(-1);
  quot = DUPInew(dnum-dden);
  rem = DUPInew(dden-1);
  DUPIdiv4(quot, rem, num, den);
  DUPIfree(rem);
  return quot;
}


DUPI DUPIrem(const DUPI num, const DUPI den)
{
  int dnum, dden;
  DUPI quot, rem;

  dnum = num->deg;
  dden = den->deg;
  if (dden < 0)
    {
      JERROR(JERROR_DIV_BY_ZERO);
      return DUPInew(-1);
    }
  if (dnum < 0) return DUPInew(-1);
  quot = DUPInew(dnum-dden);
  rem = DUPInew(dden-1);
  DUPIdiv4(quot, rem, num, den);
  DUPIfree(quot);
  return rem;
}


void DUPIrem2(DUPI x, const DUPI m)
{
  int dx, dm;
  int32 lcm, q, *xc, *mc;
  dm = DUPIdeg(m);
  dx = DUPIdeg(x);
  if (dm < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return;
  }
  
  mc = m->coeffs;
  xc = x->coeffs;
  lcm = m->coeffs[dm];
  while (dx >= dm)
  {
    q = xc[dx] / lcm;
    DUPIshift_add_raw(xc+(dx-dm), mc, mc+dm, -q);
    do { --dx; } while ((dx >= 0) && (xc[dx] == 0));
  }
  x->deg = dx;
}


void DUPIdiv4(DUPI quot, DUPI rem, const DUPI num, const DUPI den)
{
  int dden, dnum, dquot, dtmp, i;
  int32 *q, *r, *n, *d, *tmp, lcd;

  if (quot == rem || quot == den || rem == den)
  {
    JERROR(JERROR_DIV4_ARGS);
    return;
  }
  dden = den->deg;
  if (dden < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return;
  }

  dnum = num->deg;
  if (dnum < 0)  /* quotient and remainder are 0 */
  {
    quot->deg = -1;
    rem->deg = -1;
    return; 
  }

  r = rem->coeffs;
  q = quot->coeffs;
  n = num->coeffs;
  d = den->coeffs;

  if (dnum < dden) /* quotient is zero, remainder = num */
  {
    if (rem->maxdeg < dnum)
    {
      JERROR(JERROR_DEG_TOO_LOW);
      return;
    }
    if (rem != num) for (i=0; i<=dnum; i++) r[i] = n[i];
    rem->deg = dnum;
    quot->deg = -1;
    return;
  }
  dquot = dnum-dden;
  if (quot->maxdeg < dquot)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    return;
  }
  quot->deg = dquot;

  /* create some temporary space */
  tmp = r;
  if (rem->maxdeg < dnum)
    tmp = (int32*)MALLOC((1+dnum)*sizeof(int32));
  /* copy numerator into the temporary space */
  dtmp = dnum;
  if (tmp != n) for (i=dtmp; i>=0; i--) tmp[i] = n[i];

  for (i=dquot; i>=0; i--) q[i] = 0;
  lcd = d[dden];
  while (dquot >= 0)
  {
    q[dquot] = tmp[dtmp] / lcd;
    DUPIshift_add_raw(tmp+dquot, d, d+dden, -q[dquot]);
    do { --dtmp; --dquot; } while ((dquot >= 0) && (tmp[dtmp] == 0));
  }

  /* get correct value in dtmp */
  while (dtmp >= 0 && tmp[dtmp] == 0) --dtmp;
  /* copy what's in tmp into rem */
  rem->deg = dtmp;
  if (dtmp >= 0 && tmp != r) 
    for (i=dtmp; i>= 0; i--) r[i] = tmp[i];
  if (tmp != r) FREE(tmp);
}
