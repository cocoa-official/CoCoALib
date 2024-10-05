#ifndef CoCoA_ApproxPts2_H
#define CoCoA_ApproxPts2_H

//   Copyright (c)  2006,2008  John Abbott,  Anna M. Bigatti
//   Main authors: Laura Torrente (assisted by John Abbott)

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
// using std::vector;
#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class PPMonoidElem;     // fwd decl -- defined in PPMonoid.H
  class RingElem;         // fwd decl -- defined in ring.H
  class RingElemAlias; // fwd decl -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;
  class SparsePolyRing;   // fwd decl -- defined in SparsePolyRing.H


  namespace ApproxPts
  {

    typedef std::vector<RingElem> PointR; // all coords must be in same ring (which must be an ordered field)
    typedef std::vector<double> PointDbl;

    void SOI(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const std::vector<PointR>& OrigPts, const std::vector<RingElem>& OrigTolerance, ConstRefRingElem gamma);
    void SOI(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const SparsePolyRing& P, const std::vector<PointR>& OrigPts, const std::vector<RingElem>& OrigTolerance, ConstRefRingElem gamma);
    void SOITwinFloat(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const SparsePolyRing& P, const std::vector<PointR>& OrigPts, const std::vector<RingElem>& OrigTolerance, ConstRefRingElem OrigGamma);
    void SOI(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const std::vector<PointDbl>& OrigPts, const std::vector<double>& OrigTolerance, double gamma);

    void NBM(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const std::vector<PointR>& OrigPts, const std::vector<RingElem>& OrigTolerance);
    void NBM(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const SparsePolyRing& P, const std::vector<PointR>& OrigPts, const std::vector<RingElem>& OrigTolerance);
    void NBMTwinFloat(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const std::vector<PointR>& OrigPts, const std::vector<RingElem>& OrigTolerance);
    void NBM(std::vector<PPMonoidElem>& QB, std::vector<RingElem>& BBasis, std::vector<RingElem>& AlmostVanishing, const std::vector<PointDbl>& OrigPts, const std::vector<double>& OrigTolerance);

// Next fn is TEMPORARY!!!!! Added for CoCoA School 2013-06
    void VanishPoly(std::vector<PointR>& newpts,
                    std::vector<RingElem>& LowDegreePoly,
                    //std::vector<PPMonoidElem>& QB,
                    const SparsePolyRing& P,
                    const std::vector<PointR>& pts,
                    const std::vector<RingElem>& tolVector,
                    ConstRefRingElem tol);

  } // end of namespace preprocess

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ApproxPts2.H,v 1.9 2022/02/18 14:11:44 abbott Exp $
// $Log: ApproxPts2.H,v $
// Revision 1.9  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2021/01/07 14:09:11  abbott
// Summary: Corrected copyright
//
// Revision 1.7  2017/04/07 20:04:24  abbott
// Summary: Corrected fn prototype
//
// Revision 1.6  2017/03/28 12:55:10  bigatti
// -- added SOI
//
// Revision 1.5  2013/06/06 07:00:01  bigatti
// -- added SparsePolyRing in VanishPoly
//
// Revision 1.4  2013/06/05 17:44:40  abbott
// Added temporary fn VanishPoly.
//
// Revision 1.3  2012/10/24 10:33:16  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.2  2012/07/31 11:58:21  abbott
// Added full interface to NBM -- it now takes a SparsePolyRing (in which the result will reside).
//
// Revision 1.1  2008/11/23 18:58:32  abbott
// Major overhaul to preprocessing and SOI/NBM code.
// Split SOI/NBM off into a separate file.
// Preprocessing is now "rational" (but internally guided by C++ doubles).
// SOI/NBM now each have 3 similar interfaces: one purely rational, one for
// input which is represented as doubles, and one which converts the input
// to RingTwinFloat values and produces a result which is over some RingTwinFloat
// (the precision is increased automatically until an answer is obtained).
//
// Revision 1.6  2008/09/12 13:28:43  bigatti
// -- new: NBM implementation
//
// Revision 1.5  2008/06/04 18:27:37  abbott
// Modified the server interface for "SOI": it now accepts a 3rd arg (gamma).
//
// Revision 1.4  2008/05/30 14:20:43  abbott
// SOI now returns also the "almost vanishing" polynomials.
//
// Revision 1.3  2008/05/29 15:46:29  bigatti
// -- added Approximate Border Basis (by Abbott,Torrente)
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2006/12/06 17:09:37  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/10/06 10:10:56  cocoa
// Corrected position of #ifndef etc.
//
// Revision 1.2  2006/06/21 17:05:47  cocoa
// Major overhaul of approx point preprocessing algms.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.2  2006/05/22 15:52:16  cocoa
// Added preprocess-disg algorithm to ApproxPts.
// Sundry minor improvements.
//
// Revision 1.1  2006/05/12 13:16:30  cocoa
// Added functions for preprocessing approximate points.
//
//

#endif
