#ifndef CoCoA_RingWeyl_H
#define CoCoA_RingWeyl_H

//   Copyright (c)  2003-2007, 2011  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/SparsePolyRing.H"

#include <vector>
//using std::vector;

namespace CoCoA
{

  /*-----------------------------------------------------------------*/
  /** \include RingWeyl.txt  */
  /*-----------------------------------------------------------------*/

  // Pseudo-ctors for (sparse) polynomial rings.

  std::vector<symbol> WANAMES(const std::vector<symbol>& names); ///< {x,y,z,..} --> {x,y,z,..,dx,dy,dz,..}

  SparsePolyRing NewWeylAlgebra(const ring& CoeffRing, long NumTrueIndets, const std::vector<long>& ElimIndets);
  SparsePolyRing NewWeylAlgebra(const ring& CoeffRing, const std::vector<symbol>& names, const std::vector<long>& ElimIndets);
//   const RingElem& indet(const RingWeyl& RW, long var);
//   const RingElem& derivation(const RingWeyl& RW, long var);


} // end of namespace CoCoA

#endif
