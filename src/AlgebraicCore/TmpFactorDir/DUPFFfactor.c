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

#include "DUPFFfactor.h"
#include <stddef.h>
#include <stdlib.h>
#include "DUPFFsqfrd.h"
#include "berlekamp.h"

DUPFFlist DUPFFfactor(const DUPFF f)
{
  DUPFFlist sqfr_cmpts, ans, iter, iter2, cmpt_factors;

  sqfr_cmpts = DUPFFsqfrd(f);
  ans = NULL;

  /* Apply Berlekamp to each component... */
  for (iter=sqfr_cmpts; iter; iter = iter->next)
  {
    cmpt_factors = Berlekamp(iter->poly);
    for (iter2 = cmpt_factors; iter2; iter2 = iter2->next)
      iter2->deg = iter->deg;
    ans = DUPFFlist_append(ans, cmpt_factors);
  }

  DUPFFlist_dtor(sqfr_cmpts);

  return ans;
}
