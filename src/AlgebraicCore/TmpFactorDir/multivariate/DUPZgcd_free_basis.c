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

/* Function to compute a square-free GCD-free basis of a list of univariate */
/* polynomials over the integers.                                           */

/* Currently we use a naive implementation which returns only the GCD-free  */
/* basis.  A more sophisticated version should return also the factorizations*/
/* of the input polynomials in terms of the basis elements -- a list of     */
/* integer exponents suffices to describe each factorization.               */

#include "DUPZgcd_free_basis.h"
#include "DUPZgcd.h"
#include "DUPZsqfrd.h"

DUPZlist DUPZgcd_free_basis(const DUPZlist polys)
{
  DUPZlist basis, sqfrs;
  DUPZlist poly_iter, sqfr_iter, basis_iter; /* list iterators */
  DUPZ poly, gcd, zero;
  int maxdeg;

  /* Need deg bound for gcd, use maxdeg = max of degs of input polys */
  maxdeg = 0;
  for (poly_iter = polys; poly_iter; poly_iter = poly_iter->next)
    if (DUPZdeg(poly_iter->poly) > maxdeg)
      maxdeg = DUPZdeg(poly_iter->poly);
  gcd = DUPZnew(maxdeg);
  zero = DUPZnew(0); /* only used in calls to DUPZdiv4 */
  basis = NULL;

  for (poly_iter = polys; poly_iter; poly_iter = poly_iter->next)
  {
    /* We are guaranteed that each sqfr component is content-free */
    sqfrs = DUPZsqfrd(poly_iter->poly);
    for (sqfr_iter = sqfrs; sqfr_iter; sqfr_iter = sqfr_iter->next)
    {
      poly = sqfr_iter->poly;
      for (basis_iter = basis; DUPZdeg(poly) > 0 && basis_iter; basis_iter = basis_iter->next)
      {
	DUPZgcd3(gcd, poly, basis_iter->poly);
	if (DUPZdeg(gcd) == 0) continue;
	if (DUPZdeg(gcd) == DUPZdeg(basis_iter->poly))
        {
          DUPZcopy2(poly, zero); /* should really set it to 1 */
          break;
        }
	DUPZdiv4(poly, zero, poly, gcd);
//	if (DUPZdeg(gcd) == DUPZdeg(basis_iter->poly))
//	{
//	  if (mpz_cmp(DUPZlc(gcd), DUPZlc(basis_iter->poly)) != 0)
//	    DUPZcopy2(basis_iter->poly, gcd);
//	  continue;
//	}
	DUPZdiv4(basis_iter->poly, zero, basis_iter->poly, gcd);
	basis = DUPZlist_append(DUPZlist_ctor(DUPZcopy(gcd), 1), basis);
      }
      if (DUPZdeg(poly) > 0)
	basis = DUPZlist_append(DUPZlist_ctor(DUPZcopy(poly), 1), basis);
    }
    DUPZlist_dtor(sqfrs);
  }
  DUPZfree(gcd);
  return basis;
}
