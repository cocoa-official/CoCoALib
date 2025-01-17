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

#endif
