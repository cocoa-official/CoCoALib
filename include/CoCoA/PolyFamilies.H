#ifndef CoCoA_PolyFamilies_H
#define CoCoA_PolyFamilies_H

//   Copyright (c)  2017-2018  John Abbott,  Anna M. Bigatti

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

#include <vector>
// using std::vector;

namespace CoCoA
{

  // Orthogonal polynomial families
  RingElem ChebyshevPoly(long n, ConstRefRingElem x);  // 1st type
  RingElem ChebyshevPoly2(long n, ConstRefRingElem x); // 2nd type

  RingElem HermitePoly(long n, ConstRefRingElem x);    // physics
  RingElem HermitePoly2(long n, ConstRefRingElem x);   // probability

  RingElem LaguerrePoly(long n, ConstRefRingElem x);

  RingElem BesselPoly(long n, ConstRefRingElem x);

  // See Wikipedia for info about Dickson polys
  RingElem DicksonPoly(ConstRefRingElem x, long n, ConstRefRingElem alpha);  // 1st type
  RingElem DicksonPoly2(ConstRefRingElem x, long n, ConstRefRingElem alpha); // 2nd type

  // (generalized) Swinnerton-Dyer polynomial
  RingElem SwinnertonDyerPoly(ConstRefRingElem x, const std::vector<long>& L);
  RingElem SwinnertonDyerPoly(ConstRefRingElem x, const std::vector<BigRat>& L);

} // end of namespace CoCoA

#endif
