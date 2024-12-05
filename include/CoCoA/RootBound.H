#ifndef CoCoA_RootBound_H
#define CoCoA_RootBound_H

//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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

namespace CoCoA
{

  class BigRat; // fwd decl -- defined in BigRat.H

  BigRat RootBound(ConstRefRingElem f, long NumIters = -1);
  BigRat RootBound2(ConstRefRingElem f, long NumIters = -1);

  BigRat RootBound_Cauchy(ConstRefRingElem f);
  BigRat RootBound_Lagrange(ConstRefRingElem f);
  BigRat RootBound_Birkhoff(ConstRefRingElem f);
  BigRat RootBound_LMS(ConstRefRingElem f);
  BigRat RootBound_simple(ConstRefRingElem f);

  RingElem RootBoundTransform(ConstRefRingElem f);
//  RingElem RootBound_preprocess(ConstRefRingElem f);  // Should be public???

} // end of namespace CoCoA

#endif
