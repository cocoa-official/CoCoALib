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

#include "DUPZexgcd.h"

#include "DUPFF.h"
#include "DUPZcra.h"
#include "WGD.h"
#include "jaaerror.h"
#include "jalloc.h"
#include "primes.h"
#include "DUPZ_DUPFF.h"

/* This function should be in its own file. */
static int DUPZmodular_to_rational(DUPZ ans, mpz_t D, DUPZ residue, mpz_t modulus)
{
  mpz_t Q;
  int deg, i, success_flag;
  mpz_t tmp, *denom;
  size_t n;
  const size_t guard_bits = 10;

  n = mpz_sizeinbase(modulus, 2);
  if (n < 2*(guard_bits+1)) return 0; /* heuristic below which we would reject anyway */
  deg = DUPZdeg(residue);
  if (ans->maxdeg < deg) { JERROR(JERROR_DEG_TOO_LOW); return 0; }
  denom = (mpz_t*)MALLOC((1+deg)*sizeof(mpz_t));
  for (i=0; i <= deg; i++) mpz_init(denom[i]);
  mpz_init_set_ui(Q, 1);
  mpz_mul_2exp(Q, Q, n/2);

  mpz_set_ui(D, 1);
  for (i=0; i <= deg; i++)
  {
    modular_to_rational(ans->coeffs[i], denom[i], Q, residue->coeffs[i], modulus);
    /* Bail out if any one coeff has a "large" denom, or if common denom > Q. */
    if (mpz_sizeinbase(denom[i], 2) > n/2 - guard_bits) { success_flag = 0; goto tidy_mem; }
    mpz_lcm(D, D, denom[i]);
    if (mpz_cmp(D, Q) > 0) { success_flag = 0; goto tidy_mem; }
  }

  /* Adjust ans so it uses the common denominator. */
  mpz_init(tmp);
  for (i=0; i <= deg; i++)
  {
    mpz_divexact(tmp, D, denom[i]);
    mpz_mul(ans->coeffs[i], ans->coeffs[i], tmp);
  }
  mpz_clear(tmp);
  ans->deg = deg;
  success_flag = 1;
  
 tidy_mem:
  mpz_clear(Q);
  for (i=0; i <= deg; i++) mpz_clear(denom[i]);
  FREE(denom);
  return success_flag;
}


DUPZ DUPZexgcd(DUPZ fcofac, DUPZ gcofac, const DUPZ f, const DUPZ g)
{
  int nprimes, gap; /* Used for the reconstruction strategy. */
  int df, dg, dmin, p, OK;
  mpz_t den, gcd_den, fcofac_den, gcofac_den, modulus, tmp;
  DUPFF fp, gp, gcdp, fcofacp, gcofacp;
  DUPZ gcd, gcd_num, fcofac_num, gcofac_num, tmp1, tmp2, tmp3;
  FF Fp;

  if (fcofac->maxdeg < DUPZdeg(g)-1) { JERROR(JERROR_DEG_TOO_LOW); return DUPZnew(-1); }
  if (gcofac->maxdeg < DUPZdeg(f)-1) { JERROR(JERROR_DEG_TOO_LOW); return DUPZnew(-1); }
  df = DUPZdeg(f);
  dg = DUPZdeg(g);
  dmin = (df < dg)? df : dg;

/* We do not bother to check for many special cases as this would complicate */
/* the code significantly probably without any benefit for our intended use. */
/* Possible special cases: f=0, g=0, x divides f or g, consider reversing... */
/* It might also be worth removing integer content, but we do not (yet) try. */

  mpz_init_set_ui(modulus, 1);
  fp = DUPFFnew(df);
  gp = DUPFFnew(dg);
  fcofac->deg = -1;
  gcofac->deg = -1;

  p = df+dg; /* p not too big because FFctor gets costly */
  if (p > 999) p = 999;

  gcd = DUPZnew(dmin);     /* Probably too large */
  gcd_num = DUPZnew(dmin); /*                    */
  fcofac_num = DUPZnew(dg-1);
  gcofac_num = DUPZnew(df-1);
  mpz_init(den);
  mpz_init(gcd_den);
  mpz_init(fcofac_den);
  mpz_init(gcofac_den);
  mpz_init(tmp);
  tmp1 = DUPZnew(df+dg-1);
  tmp2 = DUPZnew(df+dg-1);
  
  nprimes = 0;
  gap = 0;
  while (1)
  {
    do
    {
      p = nextprime(p);
    }
    while (mpz_fdiv_ui(DUPZlc(f), p) == 0 || mpz_fdiv_ui(DUPZlc(g), p) == 0);

    nprimes++;
    gap++;
    Fp = FFctor(p);
    FFselect(Fp);
    DUPZ_to_DUPFF2(fp, f);
    DUPZ_to_DUPFF2(gp, g);
    gcdp = DUPFFexgcd(&fcofacp, &gcofacp, fp, gp);
    /* eliminate known bad modular images */
    if (DUPZdeg(gcd) >= 0 && DUPFFdeg(gcdp) > DUPZdeg(gcd))
    {
DUPFFfree(gcdp);DUPFFfree(fcofacp);DUPFFfree(gcofacp);
      FFdtor(Fp);
      continue;
    }

    /* fix leading coeff */
    if (DUPFFlc(gcdp) != 1)
    {
      DUPFFmul2ff(fcofacp, FFdiv(1, DUPFFlc(gcdp)));
      DUPFFmul2ff(gcofacp, FFdiv(1, DUPFFlc(gcdp)));
      DUPFFmul2ff(gcdp, FFdiv(1, DUPFFlc(gcdp)));
    }
    
    if (DUPFFdeg(gcdp) < DUPZdeg(gcd)) /* discard all previous modular gcds */
    {
      nprimes=0; gap=0;
      mpz_set_ui(modulus, 1);
      gcd->deg = -1;
    }

    DUPZcra(gcd, modulus, gcdp, p);
    DUPZcra(fcofac, modulus, fcofacp, p);
    DUPZcra(gcofac, modulus, gcofacp, p);
    FFdtor(Fp); /* finished with Fp now */
DUPFFfree(gcdp);DUPFFfree(fcofacp);DUPFFfree(gcofacp);
    mpz_mul_ui(modulus, modulus, p);

    mpz_set_ui(den, 0);
    if (gap*gap < nprimes) continue;
    gap = 0;
    OK = DUPZmodular_to_rational(gcd_num, gcd_den, gcd, modulus);
    if (!OK) continue;

    /* It is imperative to check the gcd really divides f and g.          */
    /* I've just wasted two days tracking this bug; previously I didn't   */
    /* perform this check and assumed that the checks below sufficed.     */
    /* I should use an intelligent way of checking for divisibility...    */
    DUPZcopy2(tmp1, f); DUPZrem2(tmp1, gcd_num);
    if (DUPZdeg(tmp1) >= 0) continue;
    DUPZcopy2(tmp1, g); DUPZrem2(tmp1, gcd_num);
    if (DUPZdeg(tmp1) >= 0) continue;

    OK = DUPZmodular_to_rational(fcofac_num, fcofac_den, fcofac, modulus);
    if (!OK) continue;

    OK = DUPZmodular_to_rational(gcofac_num, gcofac_den, gcofac, modulus);
    if (!OK) continue;

    mpz_lcm(den, gcd_den, fcofac_den);
    mpz_lcm(den, den, gcofac_den);
    mpz_divexact(tmp, den, gcd_den);
    if (mpz_cmp_ui(tmp, 1) != 0) DUPZmul2z(gcd_num, tmp);
    mpz_divexact(tmp, den, fcofac_den);
    if (mpz_cmp_ui(tmp, 1) != 0) DUPZmul2z(fcofac_num, tmp);
    mpz_divexact(tmp, den, gcofac_den);
    if (mpz_cmp_ui(tmp, 1) != 0) DUPZmul2z(gcofac_num, tmp);

    /* Check if gcd = f*fcofac + g*gcofac */
    DUPZmul3(tmp1, f, fcofac_num);
    DUPZmul3(tmp2, g, gcofac_num);
    tmp3 = DUPZadd(tmp1, tmp2); /* better to use DUPZadd3? */
    OK = DUPZequal(tmp3, gcd_num);
    DUPZfree(tmp3);
    if (OK) break; /* We have the solution. */
  }

  DUPZcopy2(fcofac, fcofac_num);
  DUPZcopy2(gcofac, gcofac_num);

  DUPFFfree(fp); DUPFFfree(gp);
  DUPZfree(gcd); DUPZfree(fcofac_num); DUPZfree(gcofac_num); DUPZfree(tmp1); DUPZfree(tmp2);
  mpz_clear(den); mpz_clear(gcd_den); mpz_clear(fcofac_den); mpz_clear(gcofac_den); mpz_clear(modulus);
  mpz_clear(tmp);

  return gcd_num;
}
