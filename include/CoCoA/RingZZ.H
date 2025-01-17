#ifndef CoCoA_RingZZ_H
#define CoCoA_RingZZ_H

//   Copyright (c)  2005-2007,2010  John Abbott, and Anna M. Bigatti

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


//---------------------------------------------------------------------------
// RingZZ() is the CoCoA ring representing the integers.


namespace CoCoA
{

  class ring; // fwd decl, defined in ring.H.
  class RingHom; // fwd decl, defined in RingHom.H.

  const ring& RingZZ();    ///< The unique ring of integers.
  bool IsZZ(const ring& R);
  RingHom ZZEmbeddingHom(const ring& codomain);
  RingHom ZZEmbeddingHom(const ring& Z, const ring& codomain);

} // end of namespace CoCoA

#endif
