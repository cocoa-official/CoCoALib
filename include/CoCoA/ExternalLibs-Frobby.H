#ifndef CoCoA_frobby_H
#define CoCoA_frobby_H

//   Copyright (c)  2008  John Abbott,  Anna M. Bigatti
//   Original author: 2008 Bjarke Hammersholt Roune

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

#include "CoCoA/PREPROCESSOR_DEFNS.H"

#ifdef CoCoA_WITH_FROBBY

#include <vector>

namespace CoCoA
{
  class ConstRefPPMonoidElem; // forward declaration -- defined in PPMonoid.H 
  class ideal; // forward declaration -- defined in ideal.H
  class RingElem; // forward declaration -- defined in ring.H

  const char* FrbVersion();

  long FrbDimension(const ideal& I);

  ideal FrbAlexanderDual(const ideal& I, ConstRefPPMonoidElem pp);
  ideal FrbAlexanderDual(const ideal& I);
  ideal FrbMaximalStandardMonomials(const ideal& I);

  void FrbIrreducibleDecomposition(std::vector<ideal>& components, const ideal& I);
  void FrbPrimaryDecomposition(std::vector<ideal>& components, const ideal& I);
  void FrbAssociatedPrimes(std::vector<ideal>& primes, const ideal& I);

  RingElem FrbMultigradedHilbertPoincareNumerator(const ideal& I);
  RingElem FrbTotalDegreeHilbertPoincareNumerator(const ideal& I);
  RingElem FrbTotalDegreeHilbertPoincareNumerator(const ideal& I, const RingElem& base);

} // end of namespace CoCoA

#endif // CoCoA_WITH_FROBBY

#endif
