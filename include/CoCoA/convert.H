#ifndef CoCoA_convert_H
#define CoCoA_convert_H

//   Copyright (c)  2007,2009,2013  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/MachineInt.H" // for IntegerCast (IsInRange)
#include "CoCoA/error.H"

#include <limits>
// using std::numeric_limits;  <--- needed only for IntegerCast
#include <type_traits>
// using std::make_unsigned_t;   <--- needed only in IntegerCast

namespace CoCoA
{

  class BigInt;  // fwd decl, defined in BigInt.H
  class BigRat;  // fwd decl, defined in BigRat.H

  class RingElemAlias; // fwd decl, defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;

  // from BigInt
  bool IsConvertible(long& lhs, const BigInt& src) noexcept;
  bool IsConvertible(int& lhs, const BigInt& src) noexcept;
  bool IsConvertible(unsigned long& lhs, const BigInt& src) noexcept;
  bool IsConvertible(unsigned int& lhs, const BigInt& src) noexcept;

  // from BigRat
  bool IsConvertible(long& lhs, const BigRat& src);          // NOT noexcept with curr impl
  bool IsConvertible(int& lhs, const BigRat& src);           // NOT noexcept with curr impl
  bool IsConvertible(unsigned long& lhs, const BigRat& src); // NOT noexcept with curr impl
  bool IsConvertible(unsigned int& lhs, const BigRat& src);  // NOT noexcept with curr impl

  // double with long, BigInt, BigRat
  bool IsConvertible(long& n, double z) noexcept;       ///< convert double z into a long, fails if the the value is not an integer or is out of range
  bool IsConvertible(BigInt& N, double z);              ///< conversion fails if the the value is not an integer
  bool IsConvertible(BigRat& Q, double z);              ///< convert double z into an exact rational ==> den will be a power of 2.
  bool IsConvertible(double& z, const BigInt& N);       ///< conversion may discard information!
  bool IsConvertible(double& z, const BigRat& Q);       ///< conversion may discard information!

  // RingElem to long, BigInt, BigRat
  bool IsConvertible(long& n, ConstRefRingElem x);
  bool IsConvertible(BigInt& N, ConstRefRingElem x);
  bool IsConvertible(BigRat& N, ConstRefRingElem x);


  template<typename ToType, typename FromType>
  inline ToType ConvertTo(const FromType& val) { ToType copy; if (!IsConvertible(copy, val))  CoCoA_THROW_ERROR1(ERR::BadConvert);  return copy; }

  template<typename ToType, typename FromType>
  inline ToType ConvertTo(const FromType& val, const CoCoA::ErrorInfo& ErrMesg) { ToType copy; if (!IsConvertible(copy, val)) ThrowException(ErrMesg); return copy; }


//------------------------------------------------------------------
// Replacement for rather poor NumericCast (now deleted).
  
// Next template adapted from https://stackoverflow.com/questions/998571/
  
/// Cast integer of type "FromType" to integer of type "ToType" checking for overflow.
/// In case of overflow, throw an error ERR::BadConvert
  template<typename FromType, typename ToType>
  constexpr ToType IntegerCastFn(FromType val)
  {
    static_assert(std::is_integral<FromType>::value,  "IntegerCast valid only for integers");
    static_assert(std::is_integral<ToType>::value,    "IntegerCast valid only for integers");
    static_assert(!std::is_same<ToType, bool>::value,    "IntegerCast valid only for integers");
    static_assert(!std::is_same<FromType, bool>::value,  "IntegerCast valid only for integers");

    if /*constexpr*/ (std::is_same<FromType, ToType>::value)  return val;
    constexpr bool FromSigned = std::is_signed<FromType>::value;
    constexpr bool ToSigned =   std::is_signed<ToType>::value;
    constexpr bool BothSigned = FromSigned && ToSigned;
    constexpr bool BothUnsigned = !FromSigned && !ToSigned;

    constexpr ToType ToMax = std::numeric_limits<ToType>::max();
    constexpr ToType ToMin = std::numeric_limits<ToType>::min();

    if /*constexpr*/ (BothUnsigned)  // need C++17
    {
       using widen = std::conditional_t<(sizeof(FromType) > sizeof(ToType)), FromType, ToType>;
       if (val > widen(ToMax))  CoCoA_THROW_ERROR1(ERR::BadConvert);
       return ToType(val);

    }
    if /*constexpr*/ (BothSigned)  // need C++17
    {
        using widen = std::conditional_t<(sizeof(FromType) > sizeof(ToType)), FromType, ToType>;
        if (val > widen(ToMax) || val < widen(ToMin))
          CoCoA_THROW_ERROR1(ERR::BadConvert);
        return ToType(val);
    }
    if /*constexpr*/ (FromSigned && !ToSigned) // need C++17
    {
//      using widen = std::make_unsigned_t<std::conditional_t<(sizeof(FromType) > sizeof(ToType)), FromType, ToType>>;
      if (val < 0 || std::make_unsigned_t<FromType>(val) > ToMax)
        CoCoA_THROW_ERROR1(ERR::BadConvert);
      return ToType(val);
    }
    // [[remaining case]] if constexpr (!FromSigned && ToSigned)
    {
//      using widen = std::make_unsigned_t<std::conditional_t<(sizeof(FromType) > sizeof(ToType)), FromType, ToType>>;
      if (val > std::make_unsigned_t<ToType>(ToMax))
        CoCoA_THROW_ERROR1(ERR::BadConvert);
      return ToType(val);
    }
  }

  template<typename ToType>
  class IntegerCast
  {
  public:
    template<typename FromType>
    IntegerCast(const FromType& val): myVal(IntegerCastFn<FromType, ToType>(val)) {}
    operator ToType() const noexcept { return myVal; }
  public: // disable copy ctor and assignment
    IntegerCast(const IntegerCast<ToType>&) = delete;
    IntegerCast<ToType>& operator=(const IntegerCast<ToType>&) = delete;
  private: // data member
    ToType myVal;
  };

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/convert.H,v 1.28 2024/07/31 07:52:46 abbott Exp $
// $Log: convert.H,v $
// Revision 1.28  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.27  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.26  2021/10/30 17:26:49  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.25  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.24  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.23  2020/06/19 19:39:21  abbott
// Summary: Now all throws go through new template fn ThrowException; seems to work
//
// Revision 1.22  2020/06/17 15:49:21  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.21  2014/04/02 15:55:45  abbott
// Summary: Reorganized file contents (more sensible order now, I hope)
// Author: JAA
//
// Revision 1.20  2014/01/28 09:51:03  abbott
// Added new template fn ConvertTo which accepts 2nd arg containing error info
// to throw in case of failure.
//
// Revision 1.19  2013/03/26 14:56:06  abbott
// Updated the conversion fns (in ptic removed procedure "convert");
// numerous consequential changes.
//
// Revision 1.18  2013/03/25 17:21:49  abbott
// Replaced ConvertTo by a simpler more obvious impl which works better
// (why didn't I think of that before?)
//
// Revision 1.17  2013/02/15 16:32:44  abbott
// Removed IsInRange (moved to MachineInt).
// Added IsConvertible(BigInt, RingElem) and IsConvertible(BigRat, RingElem).
//
// Revision 1.16  2011/11/09 15:09:08  bigatti
// -- removed fwd decl MachineInt (included)
//
// Revision 1.15  2011/11/09 15:02:29  bigatti
// -- added include MachineInt
//
// Revision 1.14  2011/11/09 13:50:01  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.13  2011/08/27 20:48:39  abbott
// Added fwd decl for MachineInt, & removed the include directive.
//
// Revision 1.12  2011/08/24 10:22:45  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.11  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.10  2011/08/12 16:07:48  abbott
// Hastily added conversion fns for BigInt (same as those for ZZ).
//
// Revision 1.9  2011/03/22 20:26:12  abbott
// Added  IsConvertible(long&, double).
// Activated commented out ConvertTo template.
// Checking in because some other files need the new code.
//
// Revision 1.8  2010/04/23 13:21:45  abbott
// Conversion from string to ZZ now defaults to decimal (and not the quirky C standard).
//
// Revision 1.7  2009/12/23 22:27:28  abbott
// Added conversions from QQ to machine integers.
//
// Revision 1.6  2009/12/23 18:53:52  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.5  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.4  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.3  2008/07/21 07:59:11  abbott
// Improved two comments.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1  2007/03/23 18:38:43  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
//

#endif
