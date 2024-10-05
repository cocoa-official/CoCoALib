#ifndef CoCoA_MachineInt_H
#define CoCoA_MachineInt_H

//   Copyright (c)  2007,2008,2024  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/long32or64.H"
#include "CoCoA/NumTheory-prime.H"

#include <limits>
//using std::numeric_limits;
#include <iosfwd>
//using std::ostream;

#define CoCoA_OLD_MACHINEINT
#ifdef CoCoA_OLD_MACHINEINT

#include "CoCoA/ULong2Long.H"
#include "CoCoA/ULongLong2LongLong.H"

namespace CoCoA
{


  // This class is deliberately inline.
  class MachineInt
  {
  public:
    MachineInt(signed char n) noexcept: myValue(n), IamNegative(n<0) {}
    MachineInt(unsigned char n) noexcept: myValue(n), IamNegative(false) {}
    MachineInt(short n) noexcept: myValue(n), IamNegative(n<0) {}
    MachineInt(unsigned short n) noexcept: myValue(n), IamNegative(false) {}
    MachineInt(int n) noexcept: myValue(n), IamNegative(n<0) {}
    MachineInt(unsigned int n) noexcept: myValue(n), IamNegative(false) {}
    MachineInt(long n) noexcept: myValue(n), IamNegative(n<0) {}
    MachineInt(unsigned long n) noexcept: myValue(n), IamNegative(false) {}
#ifdef CoCoA_32BIT_LONG
    // If long is 32 bit but size_t is 64 bit, activate these definitions;
    // **BUT** note that the value stored is only 32 bits!
    MachineInt(long long n);           // not inline to avoid include problem
    MachineInt(unsigned long long n);  // (ditto)
#endif
    MachineInt(SmallPrime n) noexcept: myValue(static_cast<long>(n)), IamNegative(false) {}
    friend bool IsZero(const MachineInt& n) noexcept;
    friend bool IsOne(const MachineInt& n) noexcept;
    friend bool IsMinusOne(const MachineInt& n) noexcept;
    friend bool IsNegative(const MachineInt& n) noexcept;
    friend bool IsSignedLong(const MachineInt& n) noexcept;
    friend bool IsUnsignedLong(const MachineInt& n) noexcept;
    friend unsigned long AsUnsignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE); /// noexcept(sizeof(unsigned long) == sizeof(unsigned long long));
    friend long AsSignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE); /// noexcept(sizeof(unsigned long) == sizeof(unsigned long long));
/////    friend bool IsSignedLongLong(const MachineInt& n) noexcept;
/////    friend bool IsUnsignedLongLong(const MachineInt& n) noexcept;
/////    friend unsigned long long AsUnsignedLongLong(const MachineInt& n) noexcept;
/////    friend long long AsSignedLongLong(const MachineInt& n) noexcept;
////    friend unsigned long negate(const MachineInt& n);
/////    friend unsigned long long ull_abs(const MachineInt& n) noexcept;
    friend unsigned long uabs(const MachineInt& n) noexcept;
  private: // data members
/////    unsigned long long myValue;
    unsigned long myValue;
    bool IamNegative;
  };

  bool IsInRange(const MachineInt& lwb, const MachineInt& val, const MachineInt& upb) noexcept;   // Checks that  lwb <= val <= upb

  inline bool IsZero(const MachineInt& n) noexcept
  { return (n.myValue == 0); }

  inline bool IsOne(const MachineInt& n) noexcept
  { return (n.myValue == 1); }

  inline bool IsMinusOne(const MachineInt& n) noexcept
  { return (n.IamNegative && (n.myValue == static_cast<unsigned long>(-1L))); }
/////  { return (n.IamNegative && (n.myValue == static_cast<unsigned long long>(-1LL))); }

  inline  bool IsNegative(const MachineInt& n) noexcept
  { return n.IamNegative; }

  inline  bool IsSignedLong(const MachineInt& n) noexcept
  { return n.IamNegative || (n.myValue <= static_cast<unsigned long>(std::numeric_limits<long>::max())); }

  inline long AsSignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE) //noexcept(sizeof(unsigned long) == sizeof(unsigned long long))
  {
    CoCoA_ASSERT(IsSignedLong(n));
    return static_cast<long>(ULong2Long(n.myValue));
/////    return static_cast<long>(ULongLong2LongLong(n.myValue));  // OVERFLOW BUG????
  }

  inline  bool IsUnsignedLong(const MachineInt& n) noexcept
  { return !(n.IamNegative) && (n.myValue <= std::numeric_limits<unsigned long>::max()); }

  inline unsigned long AsUnsignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE) //noexcept(sizeof(unsigned long) == sizeof(unsigned long long))
  {
    CoCoA_ASSERT(IsUnsignedLong(n));
    return static_cast<unsigned long>(n.myValue);
  }

/////  inline  bool IsSignedLongLong(const MachineInt& n) noexcept
/////  { return n.IamNegative || (n.myValue <= static_cast<unsigned long long>(std::numeric_limits<long long>::max())); }

/////  inline long long AsSignedLongLong(const MachineInt& n) noexcept
/////  {
/////    CoCoA_ASSERT(IsSignedLongLong(n));
/////    return ULongLong2LongLong(n.myValue);
/////  }

/////  inline  bool IsUnsignedLongLong(const MachineInt& n) noexcept
/////  { return !(n.IamNegative); }

/////  inline unsigned long long AsUnsignedLongLong(const MachineInt& n) noexcept
/////  {
/////    CoCoA_ASSERT(IsUnsignedLongLong(n));
/////    return n.myValue;
/////  }

  // inline unsigned long negate(const MachineInt& n)
  // { 
  //   CoCoA_ASSERT(n.IamNegative);
  //   return -ULong2Long(n.myValue);
  // }

/////  inline unsigned long long ull_abs(const MachineInt& n) noexcept
/////  { if (n.IamNegative) return -ULongLong2LongLong(n.myValue); else return n.myValue; }

/////  inline unsigned long ul_abs(const MachineInt& n) noexcept(sizeof(unsigned long) == sizeof(unsigned long long))
/////  {
/////    if (sizeof(unsigned long) == sizeof(unsigned long long)) return ull_abs(n);
/////    const unsigned long long a = ull_abs(n);
/////    if (a > std::numeric_limits<unsigned long>::max())
/////      CoCoA_THROW_ERROR(ERR::ArgTooBig, "uabs");
/////    return static_cast<unsigned long>(a);
/////  }

  /*TEMPORARY*/  inline unsigned long uabs(const MachineInt&n) noexcept
                 /////{return ul_abs(n);}
  {
    if (n.IamNegative)
      return -ULong2Long(n.myValue);
    else
      return n.myValue;
  }


  std::ostream& operator<<(std::ostream& out, const MachineInt& n);

} // end of namespace CoCoA

#else
  
#include "CoCoA/error.H"

namespace CoCoA
{

  // This class is deliberately inline.
  class MachineInt
  {
  public:
    MachineInt(signed char n) noexcept: myValue(n) {}
    MachineInt(unsigned char n) noexcept: myValue(n) {}
    MachineInt(short n) noexcept: myValue(n) {}
    MachineInt(unsigned short n) noexcept: myValue(n) {}
    MachineInt(int n) noexcept: myValue(n) {}
    MachineInt(unsigned int n) noexcept: myValue(n) {}
    MachineInt(long n) noexcept: myValue(n) {}
    MachineInt(unsigned long n): myValue(ULong2Long(n))
      {
        if (n > static_cast<unsigned long>(std::numeric_limits<long>::max()))
          CoCoA_THROW_ERROR1(ERR::ArgTooBig);
      }
/////    MachineInt(long long n) noexcept: myValue(n) {}
/////    MachineInt(unsigned long long n): myValue(n) { if (n > static_cast<unsigned long long>(std::numeric_limits<long long>::max())) CoCoA_THROW_ERROR1(ERR::ArgTooBig); }
    MachineInt(SmallPrime n) noexcept: myValue(static_cast<long>(n)) {}
    friend bool IsZero(const MachineInt& n) noexcept;
    friend bool IsOne(const MachineInt& n) noexcept;
    friend bool IsMinusOne(const MachineInt& n) noexcept;
    friend bool IsNegative(const MachineInt& n) noexcept;
    friend bool IsSignedLong(const MachineInt& n) noexcept;
    friend bool IsUnsignedLong(const MachineInt& n) noexcept;
    friend unsigned long AsUnsignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE); //noexcept(sizeof(unsigned long) == sizeof(unsigned long long));
    friend long AsSignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE); //noexcept(sizeof(unsigned long) == sizeof(unsigned long long));
/////    friend bool IsSignedLongLong(const MachineInt& n) noexcept;
/////    friend bool IsUnsignedLongLong(const MachineInt& n) noexcept;
/////    friend unsigned long long AsUnsignedLongLong(const MachineInt& n) noexcept;
/////    friend long long AsSignedLongLong(const MachineInt& n) noexcept;
///    friend unsigned long negate(const MachineInt& n);
    friend unsigned long uabs(const MachineInt& n) noexcept;
  private: // data members
    long myValue;
  };

  bool IsInRange(const MachineInt& lwb, const MachineInt& val, const MachineInt& upb) noexcept;   // Checks that  lwb <= val <= upb

  inline bool IsZero(const MachineInt& n) noexcept
  { return (n.myValue == 0); }

  inline bool IsOne(const MachineInt& n) noexcept
  { return (n.myValue == 1); }

  inline bool IsMinusOne(const MachineInt& n) noexcept
  { return (n.myValue == -1); }

  inline  bool IsNegative(const MachineInt& n) noexcept
  { return n.myValue < 0; }

  inline bool IsSignedLong(const MachineInt& n) noexcept
  { return (n.myValue >= numeric_limits<long>::min() && n.myValue <= numeric_limits<long>::max()); }

  inline bool IsSignedLong(const MachineInt& /*n*/) noexcept
  { return (n.myValue >= 0 && n.myValue <= numeric_limits<unsigned long>::max()); }

  inline unsigned long AsUnsignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE) //noexcept(sizeof(unsigned long) == sizeof(unsigned long long))
  {
    CoCoA_ASSERT(IsUnsignedLong(n));
    return n.myValue; // auto conversion to ulong is safe
  }

  inline long AsSignedLong(const MachineInt& n) noexcept(!CoCoA_DEBUG_MODE) //noexcept(sizeof(unsigned long) == sizeof(unsigned long long))
  {
    CoCoA_ASSERT(IsSignedLong(n));    
    return n.myValue;
  }

/////  inline bool IsSignedLongLong(const MachineInt& /*n*/) noexcept
/////  { return true; }

/////  inline bool IsSignedLongLong(const MachineInt& n) noexcept
/////  { return (n.myValue >= 0); }

/////  inline unsigned long long AsUnsignedLongLong(const MachineInt& n) noexcept(sizeof(unsigned long) == sizeof(unsigned long long))
/////  {
/////    CoCoA_ASSERT(!IsUnsignedLongLong(n));
/////    return n.myValue; // auto conversion to ulong is safe
/////  }

/////  inline long long AsSignedLongLong(const MachineInt& n) noexcept//(sizeof(unsigned long) == sizeof(unsigned long long))
/////  {
/////    return n.myValue;
/////  }

  // inline bool IsInRange(const MachineInt& lwb, const MachineInt& val, const MachineInt& upb)  noexcept  // Checks that  lwb <= val <= upb
  // {
  //   return (AsSignedLong(lwb) <= AsSignedLong(val)) &&
  //          (AsSignedLong(val) <= AsSignedLong(upb));
  // }
//???  long CheckRange(const MachineInt& val, const MachineInt& lwb, const MachineInt& upb);

  // inline unsigned long negate(const MachineInt& n)
  // { 
  //   CoCoA_ASSERT(n.IamNegative);
  //   return -ULong2Long(n.myValue);
  // }

  // BUG??? correct when value is most negative long???
  inline unsigned long uabs(const MachineInt& n) noexcept
  { if (n.myValue >= 0) return n.myValue; else return -n.myValue; }

  std::ostream& operator<<(std::ostream& out, const MachineInt& n);

} // end of namespace CoCoA
#endif


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/MachineInt.H,v 1.15 2024/07/31 07:52:46 abbott Exp $
// $Log: MachineInt.H,v $
// Revision 1.15  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.14  2024/05/03 12:06:31  abbott
// Summary: Comented out all parts to do with "long long" (redmine 1804)
//
// Revision 1.13  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.12  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.11  2021/01/07 14:15:20  abbott
// Summary: Corrected copyright
//
// Revision 1.10  2020/10/23 07:51:22  abbott
// Summary: Added IsOne and IsMinusOne
//
// Revision 1.9  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.8  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.7  2018/02/27 11:03:46  abbott
// Summary: Added conversion from SmallPrime
//
// Revision 1.6  2017/04/05 14:31:14  abbott
// Summary: Added new machineint impl
//
// Revision 1.5  2015/10/09 18:28:44  abbott
// Summary: Corrected redmine reference
//
// Revision 1.4  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.3  2013/02/15 16:31:19  abbott
// Moved IsInRange here from "convert".
//
// Revision 1.2  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.1  2011/11/09 13:38:58  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.10  2011/09/06 13:39:59  abbott
// Changed include directive: now includes ULong2Long.H instead of utils.H
//
// Revision 1.9  2011/08/27 20:48:02  abbott
// Inserted forgotten calls to ULong2Long.
//
// Revision 1.8  2010/03/05 18:39:49  abbott
// Added SmallPower function -- currently undefined behaviour if overflow occurs!!
//
// Revision 1.7  2009/10/08 13:45:58  abbott
// Added new "RoundDiv" function.
//
// Revision 1.6  2009/06/05 12:14:55  abbott
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
// Revision 1.5  2008/12/12 16:32:49  abbott
// Improved comments.
//
// Revision 1.4  2008/12/11 20:06:47  abbott
// Added IsSignedLong function; needed elsewhere, and should be present anyway.
//
// Revision 1.3  2008/12/11 10:47:36  abbott
// Fixed bug in IsZero (it appeared only when CoCoA_DEBUG was set).
// Some cleaning.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
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

#endif
