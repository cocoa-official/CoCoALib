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
