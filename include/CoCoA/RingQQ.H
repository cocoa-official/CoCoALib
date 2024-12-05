#ifndef CoCoA_RingQQ_H
#define CoCoA_RingQQ_H

//   Copyright (c)  2003--2008  John Abbott, and Anna M. Bigatti

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


#include "CoCoA/FractionField.H"

// RingQQ() is the CoCoA ring representing the rationals.
namespace CoCoA
{

  class RingHom; // fwd decl, defined in RingHom.H

  const FractionField& RingQQ(); ///< The unique field of rationals.
  bool IsQQ(const ring& R);
  RingHom QQEmbeddingHom(const ring& codomain); ///< NB this can produce a *partial* homomorphism.

} // end of namespace CoCoA

#endif
