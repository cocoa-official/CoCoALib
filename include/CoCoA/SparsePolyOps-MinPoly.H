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

#endif
