#ifndef CoCoA_RingDenseUPolyClean_H
#define CoCoA_RingDenseUPolyClean_H

//   Copyright (c)  2007  Anna Bigatti

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


namespace CoCoA
{

  class symbol;         // forward declaration -- defined in symbol.H
  class ring;           // forward declaration -- defined in ring.H
  class DenseUPolyRing;  // forward declaration -- defined in DensePolyRing.H

  // Pseudo-ctors for (sparse) polynomial rings.

  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing); ///< default indet name is x
  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing, const symbol& IndetName);
  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing, const symbol& IndetName, long MinCapacity);

} // end of namespace CoCoA

#endif
