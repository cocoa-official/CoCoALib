#ifndef CoCoA_RingDistrMPolyClean_H
#define CoCoA_RingDistrMPolyClean_H

//   Copyright (c)  2005-2007,2010-2011  John Abbott, and Anna M. Bigatti

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
//using std::vector;

namespace CoCoA
{

  class symbol;         // forward declaration -- defined in symbol.H
  class PPOrdering;     // forward declaration -- defined in PPOrdering.H
  class PPOrderingCtor; // forward declaration -- defined in PPOrdering.H
  class PPMonoid;       // forward declaration -- defined in PPMonoid.H
  class SparsePolyRing; // forward declaration -- defined in SparsePolyRing.H
  class ring;           // forward declaration -- defined in ring.H

  // Pseudo-ctors for (sparse) polynomial rings.

  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const PPMonoid& PPM);
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames); ///< Default ordering is StdDegRevLex
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, long NumIndets); ///< Default ordering is StdDegRevLex
  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& ord);

} // end of namespace CoCoA

#endif
