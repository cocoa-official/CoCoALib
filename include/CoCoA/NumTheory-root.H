#ifndef CoCoA_NumTheory_root_H
#define CoCoA_NumTheory_root_H

//   Copyright (c)  2022  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/BigInt.H"

namespace CoCoA
{

  class BigInt; // fwd decl -- defined in BigInt.H

  BigInt StarRoot(BigInt N, long UPBexp);
  long CertifyNotPower(const BigInt& N, long k);


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/NumTheory-root.H,v 1.2 2022/02/18 14:11:47 abbott Exp $
// $Log: NumTheory-root.H,v $
// Revision 1.2  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2022/01/20 19:15:57  abbott
// Summary: Added new fns StarRoot and CertifyNotPower
//
//
//

#endif
