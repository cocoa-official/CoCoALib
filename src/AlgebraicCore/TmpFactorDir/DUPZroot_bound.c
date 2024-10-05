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


#include "DUPZroot_bound.h"
#include "mpz_log.h"
#include <math.h>
/* Must define our own log for integers -- some compilers get it wrong. */
#include "logi.h"

/* This finds an upper bound on the magnitude of the roots of f using    */
/* Knuth's bound -- see "Seminumerical Algorithms" ex. 4.6.2-20 (1969 ed) */

static double DUPZroot_bound_knuth(const DUPZ f)
{
  int df, i;
  double ans, loglcf, logfi;

  df = DUPZdeg(f);
  loglcf = mpz_log(f->coeffs[df]);
  ans = -loglcf;
  for (i = 0; i < df; i++)
  {
    if (mpz_sgn(f->coeffs[i]) == 0) continue;
    logfi = mpz_log(f->coeffs[i]) - loglcf;
    if (ans < logfi/(df-i)) ans = logfi/(df-i);
  }
  return ans + logi(2);
}


/* This corresponds to exercise 4.6.2-19 in "Seminumerical Algorithms" (1969 ed) */

static double DUPZroot_bound_zassenhaus(const DUPZ f)
{
  int df, i;
  double ans, loglcf, term, log_binomial;
 
  df = DUPZdeg(f);
  loglcf = mpz_log(f->coeffs[df]);
  ans = - loglcf - df*logi(2);
  log_binomial = 0.0;
  for (i = df-1; i >= 0; i--)
  {
    log_binomial += logi(i+1) - logi(df-i);
    if (mpz_sgn(f->coeffs[i]) == 0) continue;
    term = mpz_log(f->coeffs[i]) - loglcf - log_binomial;
    if (ans < term/(df-i)) ans = term/(df-i);
  }
  return ans - log(exp(logi(2)/df)-1);

}


/* Compute both Knuth's root bound, and Zassenhaus's; and return the smaller */
double DUPZroot_bound(const DUPZ f)
{
  double K, Z;
  K = DUPZroot_bound_knuth(f);
  Z = DUPZroot_bound_zassenhaus(f);
  if (K < Z) return K;
  return Z;
}
