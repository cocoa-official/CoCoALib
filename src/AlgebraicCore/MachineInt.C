//   Copyright (c)  2007-2009  John Abbott and Anna M. Bigatti

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


#include "CoCoA/MachineInt.H"
#include "CoCoA/convert.H"

#include <iostream>

namespace CoCoA
{

  #ifdef CoCoA_OLD_MACHINEINT


#ifdef CoCoA_32BIT_LONG
  // These two are not inline to avoid a circular include problem.

  // If long is 32 bit but size_t is 64 bit, activate these definitions;
  // **BUT** note that the value stored is only 32 bits!
  MachineInt::MachineInt(long long n):
      myValue(IntegerCast<long>(n)),
      IamNegative(n<0)
  {}

  MachineInt::MachineInt(unsigned long long n):
      myValue(IntegerCast<unsigned long>(n)),
      IamNegative(false)
  {}
#endif


  std::ostream& operator<<(std::ostream& out, const MachineInt& n)
  {
    if (!out) return out;  // short-cut for bad ostreams
    using std::operator<<; // need this to avoid infinite recursion as current fn hides std::operator<<
    if (IsNegative(n))
      return out << AsSignedLong(n);
/////      return out << AsSignedLongLong(n);
    return out << AsUnsignedLong(n);
/////    return out << AsUnsignedLongLong(n);
  }


  // Checks that  lwb <= val <= upb
  bool IsInRange(const MachineInt& lwb, const MachineInt& val, const MachineInt& upb) noexcept
  {
    if (IsNegative(val))
    {
      if (!IsNegative(lwb))  return false;
      const signed long VAL = AsSignedLong(val);
/////      const signed long long VAL = AsSignedLongLong(val);
      if (AsSignedLong(lwb) > VAL ) return false;
/////      if (AsSignedLongLong(lwb) > VAL) return false;
      if (!IsNegative(upb))  return true;
      return VAL <= AsSignedLong(upb);
/////      return VAL <= AsSignedLongLong(upb);
    }
    // Here we know that val >= 0.
    if (IsNegative(upb))  return false;
    const unsigned long VAL = AsUnsignedLong(val);
/////    const unsigned long long VAL = AsUnsignedLongLong(val);
    if (AsUnsignedLong(upb) < VAL)  return false;
/////    if (AsUnsignedLongLong(upb) < VAL)  return false;
    if (IsNegative(lwb))  return true;
    return AsUnsignedLong(lwb) <= VAL;
/////    return AsUnsignedLongLong(lwb) <= VAL;
  }


#else
  
  // Alternative (faster?) defn of MachineInt
  std::ostream& operator<<(std::ostream& out, const MachineInt& n)
  {
    if (!out) return out;  // short-cut for bad ostreams
    using std::operator<<; // need this to avoid infinite recursion as current fn hides std::operator<<
    return out << AsSignedLongLong(n);
  }

  // Checks that  lwb <= val <= upb
  /*???inline???*/
  bool IsInRange(const MachineInt& lwb, const MachineInt& val, const MachineInt& upb) noexcept
  {
    return ((AsSignedLongLong(lwb) <= AsSignedLongLong(val)) &&
            (AsSignedLongLong(val) <= AsSignedLongLong(upb)));
  }

#endif

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MachineInt.C,v 1.9 2024/05/03 12:09:37 abbott Exp $
// $Log: MachineInt.C,v $
// Revision 1.9  2024/05/03 12:09:37  abbott
// Summary: Commented out all code to do with "long long (redmine 1804)
//
// Revision 1.8  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.7  2021/02/10 19:41:17  abbott
// Summary: Updated "alternative impl" of MachineInt; added noexcept (redmine 1572); some tidying
//
// Revision 1.6  2021/01/07 15:07:03  abbott
// Summary: Corrected copyright
//
// Revision 1.5  2017/04/05 14:31:46  abbott
// Summary: Added new machineint impl
//
// Revision 1.4  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2013/02/15 16:31:19  abbott
// Moved IsInRange here from "convert".
//
// Revision 1.2  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.1  2011/11/09 14:06:12  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.9  2011/08/27 21:49:48  abbott
// Added two file local fns called "cmp" (just for 2 longs or 2 unsigned longs).
// Slightly simplified defn of "cmp" for MachineInt.
//
// Revision 1.8  2011/08/23 16:17:37  abbott
// Corrected & simplified defn of RoundDiv; added comment about rounding halves.
//
// Revision 1.7  2010/03/05 18:39:49  abbott
// Added SmallPower function -- currently undefined behaviour if overflow occurs!!
//
// Revision 1.6  2009/12/23 18:53:52  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.5  2009/10/08 13:39:47  abbott
// Renamed "round" into "RoundDiv".
// Added some new versions of "RoundDiv".
// Added a test for "RoundDiv".
//
// Revision 1.4  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInt.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
// Revision 1.3  2008/12/11 10:47:36  abbott
// Fixed bug in IsZero (it appeared only when CoCoA_DEBUG was set).
// Some cleaning.
//
// Revision 1.2  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/05/21 12:57:28  abbott
// New class for passing machine integers as args; includes some simple
// operations on machine integers (cmp, gcd, IsNegative,...).  Operations
// between ZZ and machine integers modified to use the new class.  Inexact
// integer division (of a ZZ) by a negative value now triggers an error;
// new error for reporting inexact integer division by a negative value.
//
//
