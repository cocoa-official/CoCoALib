#ifndef CoCoA_RingFp_H
#define CoCoA_RingFp_H

//   Copyright (c)  2001-2009,2011  John Abbott, and Anna M. Bigatti

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
  class SmallFpImpl;  // fwd decl -- defined in SmallFpImpl.H

  QuotientRing NewRingFp(const MachineInt& p, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
  QuotientRing NewRingFp(const BigInt& P);
  QuotientRing NewRingFp(const ideal& I);

  bool IsGoodForRingFp(const MachineInt& p);
  bool IsGoodForRingFp(const BigInt& P);
  bool IsGoodForRingFp(const ideal& I);

  bool IsRingFp(const ring& R); ///< DO NOT USE this function unless you really know what you're doing!

  const SmallFpImpl& ModularArith(const ring& R);

} // end of namespace CoCoA

#endif
