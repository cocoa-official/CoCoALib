#ifndef CoCoA_RingFpDouble_H
#define CoCoA_RingFpDouble_H

//   Copyright (c)  2005-2009,2011  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/GlobalManager.H"

namespace CoCoA
{

  class BigInt;       // fwd decl -- defined in BigInt.H
  class MachineInt;   // fwd decl -- defined in MachineInt.H
  class QuotientRing; // fwd decl -- defined in QuotientRing.H
  class ideal;        // fwd decl -- defined in ring.H
  class ring;         // fwd decl -- defined in ring.H

  QuotientRing NewRingFpDouble(const MachineInt& p, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
  QuotientRing NewRingFpDouble(const BigInt& P);
  QuotientRing NewRingFpDouble(const ideal& I);

  bool IsGoodForRingFpDouble(const MachineInt& p);
  bool IsGoodForRingFpDouble(const BigInt& P);
  bool IsGoodForRingFpDouble(const ideal& I);

} // end of namespace CoCoA

#endif
