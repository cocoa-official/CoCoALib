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

#endif
