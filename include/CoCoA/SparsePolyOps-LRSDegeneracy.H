#ifndef CoCoA_SparsePolyOps_LRSDegeneracy_H
#define CoCoA_SparsePolyOps_LRSDegeneracy_H

//   Copyright (c)  2022-2023  John Abbott,  Anna M. Bigatti
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


#include "CoCoA/ring.H"
#include "CoCoA/VerificationLevel.H"

#include <vector>
// using std::vector;


namespace CoCoA
{

  // Checks for LRS-degeneracy of the given order
  bool IsLRSDegenerateOrder(ConstRefRingElem f, const unsigned long n);
  bool IsLRSDegenerateOrder(RingElem f, const unsigned long n, VerificationLevel VerLev);

  // Test whether LRS-degenerate; may be faster than LRSDegeneracyOrder.
  bool IsLRSDegenerate(RingElem f);

  // Find order of LRS-degeneracy: either least order, or all orders
  unsigned long LRSDegeneracyOrder(ConstRefRingElem f);
  std::vector<unsigned long> LRSDegeneracyOrders(ConstRefRingElem f);
  std::vector<unsigned long> LRSDegeneracyOrders(ConstRefRingElem f, VerificationLevel VerLev); // recommend lev = 3


} // end of namespace CoCoA

#endif
