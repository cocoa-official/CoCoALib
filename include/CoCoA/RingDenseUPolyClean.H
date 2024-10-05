#ifndef CoCoA_RingDenseUPolyClean_H
#define CoCoA_RingDenseUPolyClean_H

//   Copyright (c)  2007  Anna Bigatti

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


namespace CoCoA
{

  class symbol;         // forward declaration -- defined in symbol.H
  class ring;           // forward declaration -- defined in ring.H
  class DenseUPolyRing;  // forward declaration -- defined in DensePolyRing.H

  // Pseudo-ctors for (sparse) polynomial rings.

  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing); ///< default indet name is x
  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing, const symbol& IndetName);
  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing, const symbol& IndetName, long MinCapacity);

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingDenseUPolyClean.H,v 1.6 2022/02/18 14:11:48 abbott Exp $
// $Log: RingDenseUPolyClean.H,v $
// Revision 1.6  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2011/03/14 10:30:38  abbott
// Removed superfluous #include directive.
//
// Revision 1.4  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.3  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/19 10:04:23  bigatti
// -- RingDenseUPolyClean now allow to specify the MinCapacity for all
//    coeff vectors (to avoid too many reallocations)
//
// Revision 1.1  2007/10/05 15:28:56  bigatti
// -- added abstract class DenseUPolyRing for representing dense
//    univariate polynomials
// -- added concrete class RingDenseUPolyClean, the cleanest
//    implementation
// -- just for testing, still horribly inefficient and incomplete
//

#endif
