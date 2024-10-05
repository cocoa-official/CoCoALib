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
#include "DUPFF.h"
#include "jaaerror.h"
#include "jalloc.h"

int DUPFFdeg(const DUPFF f)
{
  return f->deg;
}


DUPFF DUPFFnew(const int maxdeg)
{
  DUPFF ans = (DUPFF)MALLOC(sizeof(struct DUPFFstruct));
  ans->coeffs = NULL;
  if (maxdeg >= 0) ans->coeffs = (FFelem*)MALLOC((maxdeg+1)*sizeof(FFelem));
  ans->maxdeg = maxdeg;
  ans->deg = -1;
  return ans;
}

void DUPFFfree(DUPFF x)
{
  if (x->maxdeg >= 0) FREE(x->coeffs);
  FREE(x);
}

void DUPFFswap(DUPFF x, DUPFF y)
{
  FFelem *swap1;
  int swap2;

  swap1 = x->coeffs; x->coeffs = y->coeffs; y->coeffs = swap1;
  swap2 = x->deg;    x->deg    = y->deg;    y->deg    = swap2;
  swap2 = x->maxdeg; x->maxdeg = y->maxdeg; y->maxdeg = swap2;
}


DUPFF DUPFFcopy(const DUPFF x)
{
  int dx, i;
  DUPFF ans = (DUPFF)MALLOC(sizeof(struct DUPFFstruct));

  dx = x->deg;
  ans->deg = dx;
  ans->maxdeg = dx;
  if (dx < 0) { ans->coeffs = NULL; return ans; }
  ans->coeffs = (FFelem*)MALLOC((dx+1)*sizeof(FFelem));
  for (i=0; i<=dx; i++) ans->coeffs[i] = x->coeffs[i];
  return ans;
}


void DUPFFcopy2(DUPFF dest, const DUPFF src)
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


int DUPFFequal(const DUPFF x, const DUPFF y)
{
  int dx, dy;

  dx = x->deg;
  dy = y->deg;
  if (dx < dy) return DUPFFequal(y, x);
  for (; dx > dy; dx--) if (x->coeffs[dx] != 0) return 0;
  for (; dx >= 0; dx--) if (x->coeffs[dx] != y->coeffs[dx]) return 0;
  return 1;
}


/***************************************************************************/

FFelem DUPFFlc(const DUPFF f)
{
  if (f->deg < 0) return 0;
  return f->coeffs[f->deg];
}


void DUPFFmake_monic(DUPFF f)
{
  DUPFFdiv2ff(f, DUPFFlc(f));
}

void DUPFFdiv2ff(DUPFF f, const FFelem c)
{
  int i, df;
  FFelem crecip;

  if (c == 0) { JERROR(JERROR_DIV_BY_ZERO); return; }
  if (c == 1) return;
  df = DUPFFdeg(f);
  crecip = FFdiv(1, c);
  for (i=0; i <= df; i++) f->coeffs[i] = FFmul(f->coeffs[i], crecip);
}


void DUPFFmul2ff(DUPFF f, const FFelem c)
{
  int i, df;
 
  if (c == 0) { f->deg = -1; return; }
  if (c == 1) return;
  df = DUPFFdeg(f);
  for (i=0; i <= df; i++) f->coeffs[i] = FFmul(f->coeffs[i], c);
}

/***************************************************************************/


DUPFF DUPFFadd(const DUPFF x, const DUPFF y)
{
  DUPFF ans;
  int d, dx, dy;

  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (FFadd(x->coeffs[d], y->coeffs[d]) == 0)) d--;
  ans = DUPFFnew(d);
  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] = y->coeffs[d];
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = FFadd(x->coeffs[d], y->coeffs[d]);
  return ans;
}

void DUPFFadd3(DUPFF ans, const DUPFF x, const DUPFF y)
{
  int d, dx, dy;

  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (FFadd(x->coeffs[d], y->coeffs[d]) == 0)) d--;
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
    ans->coeffs[d] = FFadd(x->coeffs[d], y->coeffs[d]);
  return;
}

/***************************************************************************/


DUPFF DUPFFsub(const DUPFF x, const DUPFF y)
{
  DUPFF ans;
  int d, dx, dy;

  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (FFsub(x->coeffs[d], y->coeffs[d]) == 0)) d--;
  ans = DUPFFnew(d);
  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] = FFsub(0, y->coeffs[d]);
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = FFsub(x->coeffs[d], y->coeffs[d]);
  return ans;
}

void DUPFFsub3(DUPFF ans, const DUPFF x, const DUPFF y)
{
  int d, dx, dy;

  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
  d = (dx>dy)?dx:dy;
  if (dx == dy) while((d >= 0) &&
		      (FFsub(x->coeffs[d], y->coeffs[d]) == 0)) d--;
  if (ans->maxdeg < d)
  {
    JERROR(JERROR_DEG_TOO_LOW);
    ans->deg = -1;
    return;
  }

  ans->deg = d;
  for (; d>dx; d--) ans->coeffs[d] = FFsub(0, y->coeffs[d]);
  for (; d>dy; d--) ans->coeffs[d] = x->coeffs[d];
  for (; d>=0; d--)
    ans->coeffs[d] = FFsub(x->coeffs[d], y->coeffs[d]);
  return;
}

/***************************************************************************/

DUPFF DUPFFmul(const DUPFF x, const DUPFF y)
{
  int dx, dy, dans;
  DUPFF ans;

  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
  if (dx < 0 || dy < 0) return DUPFFnew(-1); /* zero times anything is zero */
  dans = dx+dy;
  ans = DUPFFnew(dans);
  DUPFFmul3(ans, x, y);
  return ans;
}


void DUPFFmul3(DUPFF ans, const DUPFF x, const DUPFF y)
{
  int d, dx, dy, dans, chunk;
  FFelem p = CurrentFF.prime, shift = CurrentFF.shift;
  FFelem sum, *xi, *yj, *xlast, *xstop;

  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
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

  chunk = 0;
  if (CurrentFF.k < dx && CurrentFF.k < dy) chunk = CurrentFF.k;
  for (d=dans; d >= 0; d--)
  {
    if (d >= dx) xlast = &x->coeffs[dx]; else xlast = &x->coeffs[d];
    if (d >= dy) xi = &x->coeffs[d-dy];  else xi = &x->coeffs[0];
    if (d >= dy) yj = &y->coeffs[dy];    else yj = &y->coeffs[d];
    
    sum = 0;
    if (chunk != 0)
      for (xstop=xi+chunk; xstop <= xlast; xstop += chunk)
      {
        for (; xi < xstop; xi++, yj--)
          sum += *xi * *yj;
        if (sum >= shift) sum -= shift;
      }
    for (; xi <= xlast; xi++, yj--)
      sum += *xi * *yj;
    ans->coeffs[d] = sum%p;
  }
  ans->deg = dans;
  return;
}

#if 0
void DUPFFmul3(DUPFF ans, const DUPFF x, const DUPFF y)
{
  int i, dx, dy, dans;
  FFelem logxci, p, tmp;
  FFelem *xc, *yc, *ac, *xci, *yci, *aci, *Log, *Exp;

  if (ans == y)
  {
    if (ans == x) DUPFFsquare(ans);
    else DUPFFmul3(ans, y, x);
    return;
  }
  /* we are now sure that ans and y are different */
  dx = DUPFFdeg(x);
  dy = DUPFFdeg(y);
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

  Log = CurrentFF.LogTable;
  Exp = CurrentFF.ExpTable;
  p = CurrentFF.prime;

  dans = dx+dy;
  ac = ans->coeffs + dx;
  xc = x->coeffs;
  yc = y->coeffs;

  for (i=1; i<=dy; i++) ac[i] = 0;
  for (xci=xc+dx; xci>=xc; xci--, ac--)
  {
    if (!*xci) { *ac = 0; continue; }
    logxci = Log[*xci];
    *ac = 0;
    for (yci=yc+dy, aci=ac+dy; yci>=yc; yci--, aci--)
    {
      if (!*yci) continue;
      tmp = *aci + Exp[logxci + Log[*yci]];
      *aci = (tmp>=p)?tmp-p:tmp;
    }
  }

  ans->deg = dans;
  return;
}
#endif

/***************************************************************************/
/* This is an "in place" squaring function.                                */
/* The multiplication routine above cannot do "in place" squaring, and     */
/* while this routine can easily be modified to do general multiplication  */
/* it is noticeably slower than the routine above.                         */


void DUPFFsquare(DUPFF f)
{
  int df, da;
  FFelem p;
  unsigned int *Log;
  FFelem *Exp, *fc, *hi, *histart, *lo, *lostart, sum;

  
  if (f->deg < 0) return;    /* trivial case */
  if (f->maxdeg < 2*f->deg)  /* error case */
  {
    JERROR(JERROR_DEG_TOO_LOW);
    return;
  }
  df = DUPFFdeg(f);
  fc = f->coeffs;
  histart = fc + df;
  lostart = histart;
  Log = CurrentFF.LogTable;
  Exp = CurrentFF.ExpTable;
  p = CurrentFF.prime;
  for (da=2*df; da>=0; da--)
  {
    sum = 0;
    for (hi=histart, lo=lostart; hi > lo; hi--, lo++)
    {
      if (*hi && *lo) sum += Exp[Log[*hi]+Log[*lo]];
      if (sum >= p) sum -=p;
    }
    sum += sum; if (sum >= p) sum -=p;
    if (hi == lo  && *hi) { sum += Exp[2*Log[*hi]]; if (sum >= p) sum -=p; }
    fc[da] = sum;
    if (da > df) lostart--; else histart--;
  }

  f->deg = 2*df;
  return;
}


/****************************************************************************/
/* This "binary" powering implementation is about twice as fast as naive    */
/* iteration but more complicated; for the moment I'll keep it.             */

DUPFF DUPFFexpt(const DUPFF base, const int power)
{
  DUPFF ans;

  if (DUPFFdeg(base) < 0 || power < 0)
  {
    if (power == 0) JERROR(JERROR_ZERO_TO_POWER_ZERO);
    return DUPFFnew(-1);
  }
  if (power == 0)
  {
    ans = DUPFFnew(0);
    ans->deg = 0;
    ans->coeffs[0] = 1;
    return ans;
  }
  if (power == 1) return DUPFFcopy(base);
  ans = DUPFFnew(power*DUPFFdeg(base));
  DUPFFexpt3(ans, base, power);
  return ans;
}


void DUPFFexpt3aux(DUPFF ans, DUPFF base, int n);

void DUPFFexpt3(DUPFF ans, const DUPFF base, const int power)
{
  DUPFF f;

  if (DUPFFdeg(base) < 0 || power < 0)
  {
    if (power == 0) JERROR(JERROR_ZERO_TO_POWER_ZERO);
    ans->deg = -1;
    return;
  }
  if (ans->maxdeg < power*DUPFFdeg(base))
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
    DUPFFcopy2(ans, base);
    return;
  }

  /* The next two lines are in case ans and base are identical */
  if (ans == base) f = DUPFFcopy(base);
  else f = base;


  DUPFFexpt3aux(ans, f, power);
  if (f != base) DUPFFfree(f);
}

void DUPFFexpt3aux(DUPFF ans, DUPFF base, int n)
{
  if (n == 1)
  {
    DUPFFcopy2(ans, base);
    return;
  }
  DUPFFexpt3aux(ans, base, n/2);
  DUPFFsquare(ans);
  if ((n&1) == 0) return;
  DUPFFmul3(ans, ans, base);
}

/***************************************************************************/


void DUPFFshift_add_raw(FFelem *dest, const FFelem *src, const FFelem *srclast, const FFelem scalar)
{
  FFelem tmp;
  FFelem p = CurrentFF.prime;
  unsigned int *Log = CurrentFF.LogTable, logfactor;
  FFelem *Exp = CurrentFF.ExpTable;

  if (scalar == 0) return;
  if (Log == NULL) { for (;src <= srclast; src++, dest++) *dest = (*dest + scalar* *src)%p; return; }
  logfactor = Log[scalar];
  for (;src <= srclast; src++, dest++)
  {
    if (*src == 0) continue;
    tmp = *dest + Exp[logfactor+Log[*src]];
    *dest = (tmp>=p)?tmp-p:tmp;
  }
}

void DUPFFshift_add(DUPFF f, const DUPFF g, int deg, const FFelem coeff)
{
  int i, dg;

  dg = g->deg;
  if (dg < 0) return;
  if (coeff == 0) return;
  if (dg+deg > f->maxdeg) { JERROR(JERROR_DEG_TOO_LOW); return; }
  for (i=DUPFFdeg(f)+1; i <= dg+deg; i++) f->coeffs[i] = 0;
  f->deg = i-1;
  DUPFFshift_add_raw(f->coeffs + deg, g->coeffs, g->coeffs+dg, coeff);
  for (i=f->deg; i >= 0 && f->coeffs[i] == 0; i--) {}
  f->deg = i;
  return;
}

/***************************************************************************/

DUPFF DUPFFdiv(const DUPFF num, const DUPFF den)
{
  int dnum, dden;
  DUPFF quot, rem;
  dnum = num->deg;
  dden = den->deg;
  if (dden < 0)
    {
      JERROR(JERROR_DIV_BY_ZERO);
      return DUPFFnew(-1);
    }
  if (dnum < 0) return DUPFFnew(-1);
  quot = DUPFFnew(dnum-dden);
  rem = DUPFFnew(dden-1);
  DUPFFdiv4(quot, rem, num, den);
  DUPFFfree(rem);
  return quot;
}


DUPFF DUPFFrem(const DUPFF num, const DUPFF den)
{
  int dnum, dden;
  DUPFF quot, rem;

  dnum = num->deg;
  dden = den->deg;
  if (dden < 0)
    {
      JERROR(JERROR_DIV_BY_ZERO);
      return DUPFFnew(-1);
    }
  if (dnum < 0) return DUPFFnew(-1);
  quot = DUPFFnew(dnum-dden);
  rem = DUPFFnew(dden-1);
  DUPFFdiv4(quot, rem, num, den);
  DUPFFfree(quot);
  return rem;
}


void DUPFFrem2(DUPFF x, const DUPFF m)
{
  int dx, dm;
  FFelem lcm, q, p, *xc, *mc;
  dm = DUPFFdeg(m);
  dx = DUPFFdeg(x);
  if (dm < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return;
  }
  
  mc = m->coeffs;
  xc = x->coeffs;
  p = CurrentFF.prime;
  lcm = m->coeffs[dm];
  while (dx >= dm)
  {
    q = FFdiv(xc[dx], lcm);
    DUPFFshift_add_raw(xc+(dx-dm), mc, mc+dm, p-q);
    do { --dx; } while ((dx >= 0) && (xc[dx] == 0));
  }
  x->deg = dx;
}


void DUPFFdiv4(DUPFF quot, DUPFF rem, const DUPFF num, const DUPFF den)
{
  int dden, dnum, dquot, dtmp, i, count;
  FFelem *q, *r, *n, *d, *tmp, lcdrecip, qq;
  FFelem p = CurrentFF.prime, k = CurrentFF.k, shift = CurrentFF.shift;

  if (quot == rem || quot == den || rem == den)
  {
    JERROR(JERROR_DIV4_ARGS);
    return;
  }
  dden = den->deg;
  dnum = num->deg;
  dquot = dnum-dden;
  if (dden < 0)
  {
    JERROR(JERROR_DIV_BY_ZERO);
    return;
  }
  if (dnum < dden) /* quotient is zero, remainder = num */
  {
    DUPFFcopy2(rem, num);
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
  r = rem->coeffs;
  q = quot->coeffs;
  n = num->coeffs;
  d = den->coeffs;


  /* if necessary, create some temporary space */
  tmp = r;
  if (rem->maxdeg < dnum)
    tmp = (FFelem*)MALLOC((1+dnum)*sizeof(FFelem));
  /* copy numerator into the temporary space */
  dtmp = dnum;
  if (tmp != n) for (i=dtmp; i>=0; i--) tmp[i] = n[i];

  lcdrecip = FFdiv(1, d[dden]);
  count = 0;
  if (dden < k) k = dquot+1; /* disable "shifting" if den has low degree */
  for(; dquot >= 0; --dquot, --dtmp)
  {
    tmp[dtmp] %= p;
    if (tmp[dtmp] == 0) { q[dquot] = 0; continue; }
    q[dquot] =  FFmul(tmp[dtmp], lcdrecip);
    qq = p - q[dquot];
    for (i=0; i < dden; i++)
      tmp[i+dquot] += qq*d[i];
    if (++count < k) continue;
    count = 0;
    for (i=k; i < dden; i++)
      if (tmp[i+dquot] > shift) tmp[i+dquot] -= shift;
  }
  for(i=0; i < dden; i++) tmp[i] %= p;

  /* get correct value in dtmp */
  while (dtmp >= 0 && tmp[dtmp] == 0) --dtmp;
  /* copy what's in tmp into rem */
  if (dtmp > rem->maxdeg) { JERROR(JERROR_DEG_TOO_LOW); return; }
  rem->deg = dtmp;
  if (tmp != r) 
    for (i=dtmp; i >= 0; i--) r[i] = tmp[i];
  if (tmp != r) FREE(tmp);
}


/***************************************************************************/

DUPFF DUPFFgcd(const DUPFF fin, const DUPFF gin)
{
  int df, dg;
  DUPFF f, g;

  df = DUPFFdeg(fin);
  dg = DUPFFdeg(gin);
  if (df < 0 || dg == 0) return DUPFFcopy(gin);  /* result is g if f=0 or dg=0 */
  if (dg < 0 || df == 0) return DUPFFcopy(fin);  /* result is f if g=0 or df=0 */

  /* make copies of the input ready for Euclid's algorithm */
  f = DUPFFcopy(fin);
  g = DUPFFcopy(gin);
  DUPFFgcd2(f, g);
  if (DUPFFdeg(f) < 0)  { DUPFFfree(f); return g; }
  else                  { DUPFFfree(g); return f; }
}

void DUPFFgcd2(DUPFF f, DUPFF g) /* DESTRUCTIVE */
{
  int df, dg;
  FFelem q, lcgrecip, *fc, *gc;
  FFelem p = CurrentFF.prime;

  if (DUPFFdeg(f) < DUPFFdeg(g)) DUPFFswap(f, g);

  dg = DUPFFdeg(g);
  while (dg > 0)
  {
    fc = f->coeffs;
    gc = g->coeffs;
    lcgrecip = FFdiv(1, gc[dg]);
    df = DUPFFdeg(f);
    while (df >= dg)
    {
      q = FFmul(fc[df], lcgrecip);
      DUPFFshift_add_raw(fc+(df-dg), gc, gc+dg, p-q);
      do { --df; } while ((df >= 0) && (fc[df] == 0));
    }
    f->deg = df;
    DUPFFswap(f, g);
    dg = df;
  }

  if (dg < 0) return;
  DUPFFswap(f, g);
}


DUPFF DUPFFexgcd(DUPFF *fcofac, DUPFF *gcofac, const DUPFF f, const DUPFF g)
{
  DUPFF u, v, uf, ug, vf, vg;
  FFelem q, lcu, lcvrecip, p;
  int du, dv;

  p = CurrentFF.prime;

  u = DUPFFcopy(f);
  v = DUPFFcopy(g);
  /* If f is "smaller" than g, swap roles of f and g. */
  if (DUPFFdeg(f) < DUPFFdeg(g)) DUPFFswap(u, v);
  du = DUPFFdeg(u);  if (du < 0) du = 0; /* both inputs are zero */
  dv = DUPFFdeg(v);  if (dv < 0) dv = 0; /* one input is zero */
  
  uf = DUPFFnew(dv); uf->coeffs[0] = 1; uf->deg = 0;
  ug = DUPFFnew(du);
  vf = DUPFFnew(dv);
  vg = DUPFFnew(du); vg->coeffs[0] = 1; vg->deg = 0;

  while (DUPFFdeg(v) > 0)
  {
    dv = DUPFFdeg(v);
    lcvrecip = FFdiv(1, v->coeffs[dv]);
    while (DUPFFdeg(u) >= dv)
    {
      du = DUPFFdeg(u);
      lcu = u->coeffs[du];
      q = FFmul(lcu, lcvrecip);
      DUPFFshift_add(u, v, du-dv, p-q);
      DUPFFshift_add(uf, vf, du-dv, p-q);
      DUPFFshift_add(ug, vg, du-dv, p-q);
    }
    DUPFFswap(u, v);
    DUPFFswap(uf, vf);
    DUPFFswap(ug, vg);
  }
  if (DUPFFdeg(v) == 0)
  {
    DUPFFswap(u, v);
    DUPFFswap(uf, vf);
    DUPFFswap(ug, vg);
  }
  DUPFFfree(vf);
  DUPFFfree(vg);
  DUPFFfree(v);

  /* If f was "smaller" than g, we swapped their roles, so swap uf and ug. */
  if (DUPFFdeg(f) < DUPFFdeg(g)) DUPFFswap(uf, ug);
  *fcofac = uf;
  *gcofac = ug;
  return u;
}


FFelem DUPFFeval(const DUPFF f, FFelem x)
{
  FFelem value;
  int i;

  value = 0;
  for (i=DUPFFdeg(f); i >= 0; i--)
    value = FFadd(FFmul(value, x), f->coeffs[i]);
  return value;
}
