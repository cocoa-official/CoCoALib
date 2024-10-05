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


/* This is essentially copied from "A Course in Computational Algebraic
   Number Theory" by H Cohen, Springer GTM 138, page 125 (second corrected
   printing).  Probably the "Chinese Remainder" method given in Wang&Trager
   (SIAM J Comp vol 8, no 3,pp.300-305) is better. */

/* need stddef for NULL */
#include <stddef.h>
#include "DUPZsqfrd.h"
#include "DUPZformal_deriv.h"
#include "DUPZgcd.h"
#include "DUPZcontent.h"
#include "jaaerror.h"


/***************************************************************************/
/* Function to test whether a DUPZ is squarefree.                          */
/* Result is 1 if it is square-free, and 0 otherwise.                      */

int DUPZsqfr(const DUPZ f)
{
  DUPZ f1, tmp;
  int df, ans;

  df = DUPZdeg(f);
  f1 = DUPZnew(df);
  DUPZformal_deriv2(f1, f);
  tmp = DUPZgcd(f, f1);
  ans = DUPZdeg(tmp) == 0;
  DUPZfree(tmp);
  DUPZfree(f1);
  return ans;
}

/***************************************************************************/
/* Compute square-free decomposition (after discarding any content).       */
/* W A R N I N G: if leading coeff of f is negative then the product       */
/* W A R N I N G: of the sqfr components will be -f (and not f).           */
/*                (otherwise how to handle cases like -x^2?)               */

DUPZlist DUPZsqfrd(const DUPZ f)
{
  DUPZ PowerProd, shadow, tmp1, tmp2;
  DUPZ deriv, gcd, newshadow, junk, cmpt;
  int df, power;
  DUPZlist ans;
  mpz_t content;

  ans = NULL;
  if (DUPZdeg(f) <= 0) return ans; /* trivial case */

  df = DUPZdeg(f);

  /* Make a working copy in PowerProd, and remove any content */
  PowerProd = DUPZcopy(f);
  mpz_init(content);
  DUPZcontent(content, PowerProd);
  DUPZdiv2z(PowerProd, content);
  mpz_clear(content);

  /* some aliases for tmp1 and tmp2 to make the following more readable */
  tmp1 = DUPZnew(df);  gcd = tmp1;   newshadow = tmp1;
  tmp2 = DUPZnew(df);  junk = tmp2;  cmpt = tmp2;  deriv = tmp2;
  shadow = DUPZnew(df);
  DUPZformal_deriv2(deriv, PowerProd);     /* deriv = d/dx(PowerProd)    */
  DUPZgcd3(gcd, deriv, PowerProd);         /* gcd = gcd(deriv, original) */
  if (DUPZdeg(gcd) == 0)                   /* SPECIAL CASE: square-free  */
  {
    ans = DUPZlist_ctor(PowerProd, 1);
    goto tidy_mem;
  }
  DUPZdiv4(shadow, junk, PowerProd, gcd);  /* shadow = PowerProd/gcd     */
  DUPZswap(PowerProd, gcd);                /* PowerProd = gcd            */
  power = 0;

  while (DUPZdeg(shadow) > 0)
  {
    power++;
    DUPZgcd3(newshadow, shadow, PowerProd);
    if (DUPZdeg(shadow) != DUPZdeg(newshadow))
    {
      DUPZdiv4(cmpt, shadow, shadow, newshadow); /* shadow is trashed here */
      ans = DUPZlist_append(ans, DUPZlist_ctor(DUPZcopy(cmpt), power));
    }
    DUPZswap(shadow, newshadow);            /* shadow = newshadow */
    DUPZdiv4(PowerProd, junk, PowerProd, shadow);
  }

  DUPZfree(PowerProd);
tidy_mem:
  DUPZfree(tmp1);
  DUPZfree(tmp2);
  DUPZfree(shadow);

  return ans;
}
