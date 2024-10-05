#ifndef CoCoA_RingFqVec_H
#define CoCoA_RingFqVec_H

//   Copyright (c)  2015  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/config.H"
#include "CoCoA/BigInt.H"

namespace CoCoA
{

  class ring;  // fwd decl -- defined in ring.H
  class MachineInt; // fwd decl -- defined in MachineInt.H
  
  ring NewRingFqVec(const MachineInt& p, const MachineInt& d);
  ring NewRingFqVec(const BigInt& p, const MachineInt& d);
  ring NewRingFqVec(const MachineInt& p, const BigInt& d);
  ring NewRingFqVec(const BigInt& p, const BigInt& d);

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingFqVec.H,v 1.4 2022/02/18 14:11:48 abbott Exp $
// $Log: RingFqVec.H,v $
// Revision 1.4  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.2  2016/01/27 14:00:26  abbott
// Summary: Added new pseudo-ctors which accept BigInts
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//

#endif
