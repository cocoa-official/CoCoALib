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
   printing) */

/* need stddef for NULL */
#include <stddef.h>
#include "DUPFFsqfrd.h"
#include "DUPFFderiv.h"

/***************************************************************************/

/* Function to test whether a DUPFF is squarefree.                         */
/* Result is 1 if it is square-free, and 0 otherwise.                      */

int DUPFFsqfr(const DUPFF f)
{
  DUPFF f1, fcopy;
  int df, ans;

  df = DUPFFdeg(f);
  f1 = DUPFFnew(df);
  DUPFFformal_deriv2(f1, f);
  fcopy = DUPFFcopy(f);
  DUPFFgcd2(fcopy, f1);
  ans = DUPFFdeg(fcopy) == 0;
  DUPFFfree(fcopy);
  DUPFFfree(f1);
  return ans;
}

/***************************************************************************/


DUPFFlist DUPFFsqfrd(const DUPFF f)
{
  DUPFF PowerProd, shadow, tmp1, tmp2;
  DUPFF deriv, gcd, oldshadow, newshadow, original, junk, cmpt;
  FFelem p = CurrentFF.prime;
  long df, i, ppower, power;
  DUPFFlist ans;

  ans = NULL;
  if (DUPFFdeg(f) < 0) return ans; /* trivial case */

  df = DUPFFdeg(f);

  ppower = 1;
  PowerProd = DUPFFcopy(f);
  /* some aliases for tmp1 and tmp2 to make the following more readable */
  tmp1 = DUPFFnew(df);  deriv = tmp1; gcd = tmp1; oldshadow = tmp1; newshadow = tmp1;
  tmp2 = DUPFFnew(df);  original = tmp2;  junk = tmp2; cmpt = tmp2;
  shadow = DUPFFnew(df);
  while (DUPFFdeg(PowerProd) > 0)
  {
    DUPFFformal_deriv2(deriv, PowerProd);     /* deriv = d/dx(PowerProd)    */
    DUPFFcopy2(original, PowerProd);          /* original = PowerProd       */
    DUPFFgcd2(deriv, original);               /* gcd = gcd(deriv, original) */
    DUPFFdiv4(shadow, junk, PowerProd, gcd);  /* shadow = PowerProd/gcd     */
    DUPFFswap(PowerProd, gcd);                /* PowerProd = gcd            */
    power = 0;

    while (DUPFFdeg(shadow) > 0)
    {
      ++power;
      if (power%p == 0)
      {
	DUPFFdiv4(PowerProd, junk, PowerProd, shadow); /* PowerProd /= shadow */
	continue;
      }
      DUPFFcopy2(original, PowerProd);
      DUPFFcopy2(oldshadow, shadow);
      DUPFFgcd2(oldshadow, original);       /*  newshadow = gcd(PowerProd, shadow) */
      if (DUPFFdeg(shadow) != DUPFFdeg(newshadow))
      {
	DUPFFdiv4(cmpt, shadow, shadow, newshadow); /* shadow is trashed here */
	/* In the line below it would be slightly faster to swap the args to append */
	ans = DUPFFlist_append(ans, DUPFFlist_ctor(DUPFFcopy(cmpt), power*ppower));
      }
      DUPFFswap(shadow, newshadow);
      DUPFFdiv4(PowerProd, junk, PowerProd, shadow);
    }
    ppower *= p;
    for (i=p; i <= DUPFFdeg(PowerProd); i+=p)
      PowerProd->coeffs[i/p] = PowerProd->coeffs[i];
    PowerProd->deg /= p;
  }

  DUPFFfree(tmp1);
  DUPFFfree(tmp2);
  DUPFFfree(shadow);
  DUPFFfree(PowerProd);

  return ans;
}
