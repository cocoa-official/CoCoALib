#ifndef CoCoA_ideal_RadicalMembership_H
#define CoCoA_ideal_RadicalMembership_H

//   Copyright (c)  2017  John Abbott,  Anna Bigatti

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

  class ideal; // fwd decl -- defined in ideal.H

  bool IsInRadical(ConstRefRingElem f, const ideal& I);
  bool IsInRadical(const ideal& I, const ideal& J);
  long MinPowerInIdeal(RingElem f, const ideal& J);

} // end of namespace CoCoA

#endif
