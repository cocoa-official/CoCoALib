#ifndef CoCoA_SmallFqUtils_H
#define CoCoA_SmallFqUtils_H

//   Copyright (c)  2010  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/SmallFpImpl.H"

#include <vector>
//using std::vector;

namespace CoCoA
{
class SmallFpImpl; // fwd decl -- defd in SmallFpImpl.H

  std::vector<SmallFpImpl::value> PolyToVec(const RingElem& m, const SmallFpImpl& FFp);
  RingElem VectorToPoly(const ring& Fpx, const std::vector<SmallFpImpl::value>& v);

  std::vector<SmallFpImpl::value> FindGroupGenerator(long p, int d);
  
} // end of namespace CoCoA

#endif
