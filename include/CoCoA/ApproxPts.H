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

#endif

//  LocalWords:  PointDbl
