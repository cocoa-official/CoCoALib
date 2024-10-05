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

#include <stdio.h>

#include "WGD.h"

typedef unsigned long UL;

/* This implementation is based on the paper by Collins & Encarnacion (JSC) */
/* which explains a faster way (akin to Lehmer's integer GCD method).       */
/* Reference: JSC 20 pp.287-297 (1995).                                     */

/* This function replaces A and B by linear combinations specified in M.    */
/* new_A = M[0]*A-M[1]*B, new_B = M[2]*A-M[3]*B but if M[4] is not 1 then   */
/* new_A and new_B are negated.                                             */

static void do_lincomb(mpz_t A, mpz_t B, UL M[5])
{
  mpz_t tmp1, tmp2, tmp3;
  mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);
  mpz_mul_ui(tmp1, A, M[0]);
  mpz_mul_ui(tmp2, B, M[1]);
  mpz_mul_ui(tmp3, A, M[2]);
  mpz_sub(A, tmp1, tmp2);
  mpz_mul_ui(tmp2, B, M[3]);
  mpz_sub(B, tmp2, tmp3);
  mpz_clear(tmp1);  mpz_clear(tmp2);  mpz_clear(tmp3);
  if (M[4] != 1) return;
  mpz_neg(A, A);
  mpz_neg(B, B);
}

/* This routine performs some iterations of Euclid's algorithm on the       */
/* "leading" digits of the big number arguments.  It performs almost as many*/
/* iterations as possible before returning.  It may perform no iterations.  */

static void DPCC(UL M[5], UL a, UL b)
{
  UL u1=1, u2=0, v1=0, v2=1;
  UL q, u3, v3, c;
  int sign = 1;
  
  while(1)
  {
    q = a/b;
    c = a - q*b;
    u3 = u1+q*u2;
    v3 = v1+q*v2;
    sign = -sign;
    if (c < v3 || b-c < v2+v3) break;
    u1=u2; u2=u3;
    v1=v2; v2=v3;
    a=b;   b=c;
  }
  M[0]=u1; M[1]=v1;
  M[2]=u2; M[3]=v2;
  M[4]=sign;
  return;
}




int modular_to_rational(mpz_t num, mpz_t den, const mpz_t Q, const mpz_t res, const mpz_t mod)
{
  mpz_t u_cofac, u, v_cofac, v, tmp, quot, rem, absv_cofac;
  int error;
  UL U, V;
  UL M[5];
  const size_t W = 8*sizeof(UL); /* number of bits in a UL */
  size_t v_bits, u_bits, Q_bits, v_cofac_bits, M_bits;

  /* Error return if input values are stupid */
  if (mpz_sgn(mod) == 0 || mpz_sgn(Q) <= 0) return 1;
  mpz_init(v);
  mpz_mul_2exp(v, Q, 1);
  if (mpz_cmpabs(v, mod) >= 0) { mpz_clear(v); return 1; }
  mpz_mod(v, res, mod); /* make sure v is reduced modulo mod */

  /* Dispose of the trivial case where v is small (and positive). */
  v_bits = mpz_sizeinbase(v, 2);
  M_bits = mpz_sizeinbase(mod, 2);
  Q_bits = mpz_sizeinbase(Q, 2);
  if (v_bits + Q_bits < M_bits - 1)
  {
    mpz_set_ui(den, 1); mpz_set(num, v); mpz_clear(v); return 0;
  }
  mpz_init(u);
  mpz_abs(u, mod);
  /* Dispose of the trivial case where v-M is small (and negative). */
  if (v_bits >= M_bits-1)
  {
    mpz_sub(v, v, u); /* v = v-mod */
    v_bits = mpz_sizeinbase(v, 2);
    if (v_bits + Q_bits < M_bits - 1)
    {
      mpz_set_ui(den, 1); mpz_set(num, v); mpz_clear(v); mpz_clear(u); return 0;
    }
    mpz_add(v, v, u); /* restore value of v */
  }

  /* Not a trivial case, so perform main algorithm. */
  mpz_init_set_ui(u_cofac, 0);
  mpz_init_set_ui(v_cofac, 1);

  mpz_init(tmp);
  mpz_init(quot);
  mpz_init(rem);
  mpz_init_set_ui(absv_cofac, 1);
  while (mpz_sgn(v) != 0 && mpz_cmp(absv_cofac, Q) <= 0)
  {
    /* Here we have u > v > 0 always */
    u_bits = mpz_sizeinbase(u, 2);
    v_bits = mpz_sizeinbase(v, 2);
    v_cofac_bits = mpz_sizeinbase(v_cofac, 2);

    M[1]=0;
    /* Decide whether to do some single digit iterations: only if           */
    /* u and v are both large, of a similar size, and v_cofac is far from Q */
    /* NB: The difference u_bits - v_bits is guaranteed non-negative.       */
    /* NB: W/2-2 is guaranteed positive, so no problems with unsigned values*/
    if (u_bits > W && u_bits - v_bits < W/2-2 && v_cofac_bits + W/2 < Q_bits)
    {
      size_t shift = u_bits-W;
      mpz_tdiv_q_2exp(tmp, u, shift); U = mpz_get_ui(tmp);
      mpz_tdiv_q_2exp(tmp, v, shift); V = mpz_get_ui(tmp);
      DPCC(M, U, V);
    }
    if (M[1]!=0) /* true only if DPCC did at least 1 iteration */
    {
      do_lincomb(u, v, M);
      do_lincomb(u_cofac, v_cofac, M);
      mpz_abs(absv_cofac, v_cofac);
      continue;
    }

    /* We were not able to do a single digit step so do a full step. */
    mpz_fdiv_qr(quot, rem, u, v);
    mpz_mul(tmp, quot, v_cofac);
    mpz_sub(tmp, u_cofac, tmp);
    mpz_set(u_cofac, v_cofac);
    mpz_set(v_cofac, tmp);
    mpz_set(u, v);
    mpz_set(v, rem);
    mpz_abs(absv_cofac, v_cofac);
  }

  mpz_gcd(tmp, u, u_cofac); /* this GCD is guaranteed to be positive */
  error = (mpz_cmp_ui(tmp, 1) != 0); /* non-zero if our result will be wrong */
  if (mpz_sgn(u_cofac) > 0) { mpz_set(num, u); mpz_set(den, u_cofac); }
  else                      { mpz_neg(num, u); mpz_neg(den, u_cofac); }

  /* Check that the numerator is small enough to guarantee uniqueness.*/
  if (!error)
  {
    mpz_mul(tmp, u, Q);
    mpz_mul_2exp(tmp, tmp, 1);
    error = (mpz_cmpabs(tmp, mod) >= 0);
  }

  mpz_clear(u_cofac);  mpz_clear(u);
  mpz_clear(v_cofac);  mpz_clear(v);
  mpz_clear(quot); mpz_clear(rem);
  mpz_clear(tmp); mpz_clear(absv_cofac);

  return error;
}


int modular_to_rational2(mpz_t num, mpz_t den, int log2_den, const mpz_t res, const mpz_t mod)
{
  mpz_t Q;
  int ans;

  mpz_init_set_ui(Q, 1);
  mpz_mul_2exp(Q, Q, log2_den);
  ans = modular_to_rational(num, den, Q, res, mod);
  mpz_clear(Q);
  return ans;
}
