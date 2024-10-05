#ifndef CoCoA_ApproxPts_H
#define CoCoA_ApproxPts_H

//   Copyright (c)  2006,2008,2013  John Abbott,  Anna M. Bigatti
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

  class RingElem;     // fwd decl -- defined in ring.H

  namespace ApproxPts
  {

    typedef std::vector<RingElem> PointR; // all coords must be in same ring (which must be an ordered field)

  } // end of namespace ApproxPts;

  void PreprocessPts(std::vector<ApproxPts::PointR>& NewPts,
                     std::vector<long>& weights,
                     const std::vector<ApproxPts::PointR>& OrigPts,
                     std::vector<RingElem> tolerance);

  void PreprocessPtsGrid(std::vector<ApproxPts::PointR>& NewPts,
                         std::vector<long>& weights,
                         const std::vector<ApproxPts::PointR>& OrigPts,
                         std::vector<RingElem> tolerance);

  void PreprocessPtsAggr(std::vector<ApproxPts::PointR>& NewPts,
                         std::vector<long>& weights,
                         const std::vector<ApproxPts::PointR>& OrigPts,
                         std::vector<RingElem> tolerance);

  void PreprocessPtsSubdiv(std::vector<ApproxPts::PointR>& NewPts,
                           std::vector<long>& weights,
                           const std::vector<ApproxPts::PointR>& OrigPts,
                           std::vector<RingElem> tolerance);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ApproxPts.H,v 1.15 2022/02/18 14:11:44 abbott Exp $
// $Log: ApproxPts.H,v $
// Revision 1.15  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.14  2021/01/07 14:08:49  abbott
// Summary: Corrected copyright
//
// Revision 1.13  2013/03/27 18:24:33  abbott
// Added approx point preprocessing to C5; also changed names of the fns, and updated doc.
//
// Revision 1.12  2011/03/11 11:06:48  bigatti
// -- changed size_t --> long
//
// Revision 1.11  2008/11/24 17:14:34  abbott
// Removed unused typedef.
//
// Revision 1.10  2008/11/23 18:58:32  abbott
// Major overhaul to preprocessing and SOI/NBM code.
// Split SOI/NBM off into a separate file.
// Preprocessing is now "rational" (but internally guided by C++ doubles).
// SOI/NBM now each have 3 similar interfaces: one purely rational, one for
// input which is represented as doubles, and one which converts the input
// to RingTwinFloat values and produces a result which is over some RingTwinFloat
// (the precision is increased automatically until an answer is obtained).
//
// Revision 1.9  2008/11/20 10:47:30  abbott
// Now the code actually computes the list of weights.
//
// Revision 1.8  2008/11/20 10:03:12  abbott
// Preprocessing fns now return the list of weights associated to each representative.
//
// Revision 1.7  2008/11/06 12:50:44  abbott
// Moved definitions of square and round to utils.H from ApproxPts.H
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

//  LocalWords:  PointDbl
