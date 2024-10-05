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

#include "DUPZgcd.h"
#include "DUPZcontent.h"
#include "DUPFF.h"
#include "DUPZ_DUPFF.h"
#include "jaaerror.h"
#include "primes.h"
#include "DUPZcra.h"
#include "FF.h"


DUPZ DUPZgcd(const DUPZ fin, const DUPZ gin)
{
  int df, dg, dmin;
  DUPZ ans;

  df = DUPZdeg(fin);
  dg = DUPZdeg(gin);
  if (df < 0) return DUPZcopy(gin);
  if (dg < 0) return DUPZcopy(fin);

  dmin = (df<dg)?df:dg;
  ans = DUPZnew(dmin);
  DUPZgcd3(ans, fin, gin);
  return ans;
}

void DUPZgcd3(DUPZ ans, const DUPZ fin, const DUPZ gin)
{
  int df, dg, dmax, p, changed, xpowerf, xpowerg, xpower, reversed;
  mpz_t ac, fc, gc, tcans, lcans, modulus, excess_content;
  FFelem lch;
  DUPFF fp, gp;
  DUPZ f, g, copy;
  FF Fp;

  if (ans->maxdeg < 0) { JERROR(JERROR_DEG_TOO_LOW); return; }
  if (DUPZdeg(fin) < 0)  { DUPZcopy2(ans, gin); return; } /* 2 trivial cases */
  if (DUPZdeg(gin) < 0)  { DUPZcopy2(ans, fin); return; } /*                 */

  /* Remove powers of x from f and g */
  f = DUPZcopy(fin);
  g = DUPZcopy(gin);
  for (xpowerf=0; mpz_sgn(f->coeffs[xpowerf]) == 0; xpowerf++) {}
  for (xpowerg=0; mpz_sgn(g->coeffs[xpowerg]) == 0; xpowerg++) {}
  DUPZshift(f, -xpowerf); df = DUPZdeg(f);
  DUPZshift(g, -xpowerg); dg = DUPZdeg(g);
  xpower = (xpowerf < xpowerg)? xpowerf : xpowerg;

  /* Now compute integer content of the answer in ac */
  mpz_init(ac);
  mpz_init(fc);
  mpz_init(gc);
  DUPZcontent(fc, fin);  DUPZdiv2z(f, fc);      /* primitive parts in f, g */
  DUPZcontent(gc, gin);  DUPZdiv2z(g, gc);      /*                         */
  mpz_gcd(ac, fc, gc);
  mpz_clear(fc);
  mpz_clear(gc);

  if (df == 0 || dg == 0)                       /* another trivial case */
  {
    mpz_set(ans->coeffs[0], ac);
    ans->deg = 0;
    DUPZshift(ans, xpower);
    goto tidy_mem;
  }

  /* lcans is multiplicative upper bound for lc of ans */
  mpz_init_set_ui(lcans, 1);
  mpz_init_set_ui(tcans, 1);
  reversed = 0;
  if (mpz_cmp_si(DUPZlc(f), 1) == 0 ||
      mpz_cmp_si(DUPZlc(f), -1) == 0 ||
      mpz_cmp_si(DUPZlc(g), 1) == 0 ||
      mpz_cmp_si(DUPZlc(g), -1) == 0) goto no_reverse;
  if (mpz_cmp_si(f->coeffs[0], 1) == 0 ||
      mpz_cmp_si(f->coeffs[0], -1) == 0 ||
      mpz_cmp_si(g->coeffs[0], 1) == 0 ||
      mpz_cmp_si(g->coeffs[0], -1) == 0) goto reverse;
  mpz_gcd(lcans, DUPZlc(f), DUPZlc(g));
  mpz_gcd(tcans, f->coeffs[0], g->coeffs[0]);
  if (mpz_cmp(lcans, tcans) <= 0) goto no_reverse;
reverse:
  DUPZreverse(f);
  DUPZreverse(g);
  mpz_set(lcans, tcans);
  reversed = 1;
no_reverse:

  mpz_init_set_ui(modulus, 1);
  dmax = (df > dg)? df : dg;
  fp = DUPFFnew(dmax);   /* these are used as workspace */
  gp = DUPFFnew(dmax);   /*                             */
  copy = DUPZnew(dmax);  /*                             */
  p = df+dg; /* p not too big because FFctor gets costly */
  if (p > 999) p = 999;

  mpz_init(excess_content);
  ans->deg = -1; /* clear ans to zero */

  while (1)
  {
    /* If no small good prime exists, the "while" condition will div by 0 */
    do
    {
      p = nextprime(p);  /* becomes zero if we exhaust the small primes */
    }
    while (mpz_fdiv_ui(DUPZlc(f), p) == 0 || mpz_fdiv_ui(DUPZlc(g), p) == 0);

    Fp = FFctor(p);
    FFselect(Fp);
    DUPZ_to_DUPFF2(fp, f);
    DUPZ_to_DUPFF2(gp, g);
    DUPFFgcd2(fp, gp); /* destructive */
    if (DUPFFdeg(fp) == 0)
    {
      mpz_set_ui(ans->coeffs[0], 1);
      ans->deg = 0;
      FFdtor(Fp);
      break;
    }

    /* eliminate known bad modular images */
    if (DUPZdeg(ans) > 0 && DUPFFdeg(fp) > DUPZdeg(ans))
    {
      FFdtor(Fp);
      continue;
    }

    /* fix leading coeff */
    lch = mpz_fdiv_ui(lcans, p); 
    if (DUPFFlc(fp) != lch) DUPFFmul2ff(fp, FFdiv(lch, DUPFFlc(fp)));

    if (DUPFFdeg(fp) < DUPZdeg(ans)) /* discard all previous modular gcds */
    {
      mpz_set_ui(modulus, 1);
      ans->deg = -1;
    }

    changed = DUPZcra(ans, modulus, fp, p);
    FFdtor(Fp); /* finished with Fp now */
    mpz_mul_ui(modulus, modulus, p);
    if (changed) continue;
    if (mpz_cmp(lcans, modulus) > 0) continue;

    /* ans is almost certainly correct now, so check it */
    /* firstly remove any content (from forcing too big a leading coeff) */
    DUPZcontent(excess_content, ans);
    DUPZdiv2z(ans, excess_content);
    
    /* We could check divisibility of trailing coeffs here,               */
    /* but as we probably have the right answer anyway, it's not worth it */
    
    DUPZcopy2(copy, f); DUPZrem2(copy, ans);
    if (DUPZdeg(copy) >= 0) goto false_alarm;
    DUPZcopy2(copy, g); DUPZrem2(copy, ans);
    if (DUPZdeg(copy) < 0) break;

  false_alarm:
    /* it wasn't correct, so restore value and keep going */
    DUPZmul2z(ans, excess_content);
  }
  /* reverse if nec, put in the integer content, and restore power of x */
  if (reversed) DUPZreverse(ans);
  /* Fiddle content so result has positive leading coeff. */
  if (mpz_sgn(DUPZlc(ans)) < 0) mpz_neg(ac, ac);
  if (mpz_cmp_ui(ac, 1) != 0) DUPZmul2z(ans, ac);
  DUPZshift(ans, xpower);

  mpz_clear(lcans); mpz_clear(tcans);
  mpz_clear(excess_content);
  mpz_clear(modulus);
  DUPFFfree(fp); DUPFFfree(gp);
  DUPZfree(copy);

tidy_mem:
  DUPZfree(f);
  DUPZfree(g);
  mpz_clear(ac);
}
