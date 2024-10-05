#ifndef POINCARECPP_H
#define POINCARECPP_H

//   Copyright (c)  2007  Anna Bigatti,  John Abbott

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

/*-----------------------------------------------------------------------------
	Part of project CoCoA
-----------------------------------------------------------------------------*/

#include "AnnaUtils.h"
#include "poincare.h"  // "TermList.h"
#include "CoCoA/DenseUPolyRing.H"

namespace CoCoA
{

  class PPMonoid;       // fwd decl -- defined in PPMonoid.H
  class SparsePolyRing; // fwd decl -- defined in SparsePolyRing.H

/**********************************************************************/
  DenseUPolyRing StartPoincareQQt(int IndetsNo);
  void StartPoincare(int IndetsNo);
  void MakeHPPowerList(std::vector<RingElem>& PL, const DenseUPolyRing& HSRing, size_t MaxDeg);
  RingElem TLPoincareNumeratorCPP(const DenseUPolyRing& HSRing, TermList theTList);
  RingElem TLPoincareNumeratorCPP(const SparsePolyRing& HSRing, const PPMonoid& PPM, TermList theTList); // HSRing for output, PPM for input
/**********************************************************************/
}

#endif /* POINCARECPP_H */
