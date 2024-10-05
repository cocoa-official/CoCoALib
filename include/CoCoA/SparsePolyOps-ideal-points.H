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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-ideal-points.H,v 1.3 2022/02/18 14:11:49 abbott Exp $
// $Log: SparsePolyOps-ideal-points.H,v $
// Revision 1.3  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2019/10/15 12:57:55  bigatti
// -- renamed files for ideals
//
// Revision 1.1  2018/04/06 14:59:57  bigatti
// -- renamed IdealOfPoints.H
//
// Revision 1.4  2017/03/01 18:35:00  bigatti
// -- added IdealOfProjectivePoints
//
// Revision 1.3  2014/04/17 13:37:19  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.2  2013/05/20 16:00:50  abbott
// Commented out fn BM -- gave a name clash in IdealOfPoints.C
//
// Revision 1.1  2013/01/21 15:30:36  abbott
// Renamed files called BuchbergerMoeller* to IdealOfPoints*.
//
//

#endif
