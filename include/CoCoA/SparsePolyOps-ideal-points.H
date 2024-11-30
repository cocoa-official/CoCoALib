#ifndef CoCoA_IdealOfPoints_H
#define CoCoA_IdealOfPoints_H

//   Copyright (c)  2012-2013  John Abbott, and Anna M. Bigatti

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

#include <vector>
//using std::vector;

namespace CoCoA
{

  class ConstMatrixView; // fwd decl -- defined in MatrixView.H
  class SparsePolyRing;  // fwd decl -- defined in PolyRing.H
  class RingElem;        // fwd decl -- defined in ring.H
  class ideal;           // fwd decl -- defined in ideal.H

  // A set of points is represented as a matrix whose i-th row contains the coords of the i-th point.
///???  std::vector<RingElem> BM(const SparsePolyRing& P, const ConstMatrixView& pts);

  ideal IdealOfPoints(const SparsePolyRing& P, const ConstMatrixView& pts);
  ideal IdealOfProjectivePoints(const SparsePolyRing& P, const ConstMatrixView& pts);

} // end of namespace CoCoA

#endif
