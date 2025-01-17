#ifndef CoCoA_HomomorphismOps_H
#define CoCoA_HomomorphismOps_H

//   Copyright (c)  2017  John Abbott, Anna Bigatti

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

  class RingHom;  // fwd decl -- defined in RingHom.H
  class ideal;    // fwd decl -- defined in ideal.H

  ideal ker(const RingHom& phi);
  ideal ker_H(const RingHom& phi);
  bool IsInjective(const RingHom& phi);
  bool IsSurjective(const RingHom& phi);
  bool IsInImage(const RingHom& phi, ConstRefRingElem y);

  RingElem preimage(const RingHom& phi, ConstRefRingElem y); // a preimage, throws exception if none exists
  RingElem preimage0(const RingHom& phi, ConstRefRingElem y); // a preimage, or 0 if none exists

  
} // end of namespace CoCoA

#endif
