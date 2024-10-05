#ifndef CoCoA_SparsePolyOps_MinPoly_H
#define CoCoA_SparsePolyOps_MinPoly_H

//   Copyright (c)  2016-2018  John Abbott and Anna M. Bigatti

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


// Header file for functions MinPoly etc.
#include "CoCoA/ring.H"
#include "CoCoA/VerificationLevel.H"

#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{

  class ConstMatrixView;  // fwd decl -- defined in matrix.H
//  class ConstRefRingElem; // fwd decl -- defined in ring.H
  class RingElem;         // fwd decl -- defined in ring.H
  class RingElemAlias;    // fwd decl -- defined in ring.H
  class matrix;           // fwd decl -- defined in matrix.H
  class ideal;            // fwd decl -- defined in ideal.H
  class PPMonoidElem;     // fwd decl -- defined in PPMonoid.H
  class PolyRing;         // fwd decl -- defined in PolyRing.H
  class VerificationLevel;// fwd decl -- defined in VerificationLevel.H

  matrix FrobeniusMat(const ideal& I);
  matrix FrobeniusMat(const ideal& I, const std::vector<PPMonoidElem>& QB2);

  RingElem MinPoly(ConstMatrixView M, ConstRefRingElem x);
  RingElem MinPoly(ConstRefRingElem f, ConstRefRingElem x);

  RingElem MinPoly(ConstRefRingElem f, ConstRefRingElem x);
  RingElem MinPolyMat(ConstRefRingElem f, ConstRefRingElem x);
  RingElem MinPolyDef(ConstRefRingElem f, ConstRefRingElem x);
  RingElem MinPolyElim(ConstRefRingElem f, ConstRefRingElem x);

  //  RingElem CharPolyQuot(ConstRefRingElem f_orig, const ideal& I, ConstRefRingElem x);
  RingElem MinPolyQuot(ConstRefRingElem f, const ideal& I, ConstRefRingElem x, VerificationLevel NumChecks=guaranteed());
//  RingElem MinPolyQuot(ConstRefRingElem f, const ideal& I, ConstRefRingElem x);
  RingElem MinPolyQuotMat(ConstRefRingElem f, const ideal& I, ConstRefRingElem x);
  RingElem MinPolyQuotDef(ConstRefRingElem f, const ideal& I, ConstRefRingElem x);
  RingElem MinPolyQuotDefLin(ConstRefRingElem f, const ideal& I, ConstRefRingElem x);
  RingElem MinPolyQuotElim(ConstRefRingElem f, const ideal& I, ConstRefRingElem x);

  std::vector<RingElem> ShapeLemma(const ideal& I);

  /////////// ??????? 
  RingElem LiftPolyFromFpxToQQx(const PolyRing& QQx, ConstRefRingElem f);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-MinPoly.H,v 1.4 2024/03/15 16:09:21 bigatti Exp $
// $Log: SparsePolyOps-MinPoly.H,v $
// Revision 1.4  2024/03/15 16:09:21  bigatti
// Summary: just a comment
//
// Revision 1.3  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2019/12/21 16:46:49  abbott
// Summary: Changed MinPolyQuot to single fn with optional param
//
// Revision 1.1  2018/04/06 15:09:39  bigatti
// -- renamed MinPoly.H
//
// Revision 1.13  2018/03/15 14:20:42  bigatti
// -- removed empty space
//
// Revision 1.12  2018/03/13 18:03:52  bigatti
// -- MinPolyModular: now using VerificationLevel class
//
// Revision 1.11  2018/03/13 17:43:38  bigatti
// -- now MinPolyQuot takes a verification level
//    instead of being called MinPolyQuotHeuristic
//
// Revision 1.10  2018/02/22 16:51:09  bigatti
// -- added MinPolyQuotHeuristic
//
// Revision 1.9  2018/01/17 10:52:18  abbott
// Summary: Corrected includes and fwd decls
//
// Revision 1.8  2017/06/26 13:17:24  bigatti
// -- exported utility function LiftPolyFromFpxToQQx
//
// Revision 1.7  2016/10/27 13:07:24  bigatti
// -- added MinPolyQuotDefLin
// -- cleaned up code, unified notation
//
// Revision 1.6  2016/10/24 12:10:35  bigatti
// -- added ShapeLemma (first draft)
//
// Revision 1.5  2016/06/27 14:50:28  bigatti
// -- now FrobeniusMat may take two args
//
// Revision 1.4  2016/06/20 15:24:42  bigatti
// -- renamed MinPolyXX --> MinPolyQuotXX
// -- added MinPolyDef(f, x)
// -- added default to MinPolyDef in both cases
//
// Revision 1.3  2016/06/10 15:55:00  bigatti
// -- only one arg for FrobeniusMat
//
// Revision 1.2  2016/04/14 11:34:58  bigatti
// -- added FrobeniusMat
//
// Revision 1.1  2016/03/18 11:54:00  bigatti
// -- first import
//

#endif
