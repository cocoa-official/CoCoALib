#ifndef CoCoA_RingExtAlg_H
#define CoCoA_RingExtAlg_H

//   Copyright (c)  2018  John Abbott,  Anna Bigatti

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


//#include "CoCoA/SparsePolyRing.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

  class ring;            // fwd decl -- defined in ring.H
  class SparsePolyRing;  // fwd decl -- defined in SparsePolyRing.H
  class symbol;          // fwd decl -- defined in symbol.H

  // Pseudo-ctors for (sparse) polynomial rings.

  SparsePolyRing NewExtAlgebra(const ring& CoeffRing, long NumIndets);
  SparsePolyRing NewExtAlgebra(const ring& CoeffRing, const std::vector<symbol>& names);


} // end of namespace CoCoA

#endif
