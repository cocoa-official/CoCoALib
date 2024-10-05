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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-LRSDegeneracy.H,v 1.9 2024/01/25 08:40:15 abbott Exp $
// $Log: SparsePolyOps-LRSDegeneracy.H,v $
// Revision 1.9  2024/01/25 08:40:15  abbott
// Summary: Improved comments; minor cleaning
//
// Revision 1.8  2023/12/09 20:29:49  abbott
// Summary: Split "eval" fns out of SparsePolyOps-cyclotomic, into SparsePolyOps-eval; added example & doc
//
// Revision 1.7  2023/10/05 18:06:30  abbott
// Summary: Added LRSDegeneracyOrders with VerifLev
//
// Revision 1.6  2023/03/06 22:06:56  abbott
// Summary: Added new fn LRSDegeneracyOrders (returns a vector)
//
// Revision 1.5  2023/02/23 20:49:10  abbott
// Summary: Added new fn LRSDegeneracyOrder (was IsLRSDegenerate); major code revision
//
// Revision 1.4  2023/01/31 12:33:34  abbott
// Summary: Added IsLRSDegenerateOrderMod; some tidying
//
// Revision 1.3  2022/11/30 15:19:43  abbott
// Summary: Removed unused flag arg
//
// Revision 1.2  2022/11/24 21:11:29  abbott
// Summary: Added copyright notice
//


#endif
