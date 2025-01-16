#ifndef CoCoA_SparsePolyOps_cyclotomic_H
#define CoCoA_SparsePolyOps_cyclotomic_H

//   Copyright (c)  2022  John Abbott,  Anna M. Bigatti
//   Original author: Nico Mexis (2022)

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


#include "CoCoA/factorization.H"
#include "CoCoA/ring.H"

#include <vector>
// using std::vector;

namespace CoCoA
{

  RingElem cyclotomic(long n, ConstRefRingElem x);  ///< need n >= 1 and x must be an indet

  //----------------------
  // Below are fns to detect if a poly is cyclo

  // If f is (monic) cyclo, returns its index; o/w returns 0
  unsigned long CyclotomicTest(ConstRefRingElem f);  // 

  // If f cyclo, return its index; o/w return some integer
  unsigned long CyclotomicIndex(ConstRefRingElem f);


  // Return list of indexes of cyclo factors; may contain some false positives (usu. only 3, 4 and/or 6)
  std::vector<long> CyclotomicFactorIndexes(ConstRefRingElem f);
  std::vector<long> CyclotomicFactorIndexes2(ConstRefRingElem f);

  // Returns a factorization containing all cyclotomic factors of f
  // (non-cyclo factors are discarded; so RemainingFactor is always 1)
  // Factorization is empty iff f has no cyclotomic factors.
  factorization<RingElem> CyclotomicFactors(ConstRefRingElem f);

  // Not sure if I really want to export this fn
  factorization<RingElem> CyclotomicFactors_BeukersSmyth(ConstRefRingElem f);

} // end of namespace CoCoA

#endif
