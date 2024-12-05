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

  template<typename ToType, typename FromType>
  inline ToType ConvertTo(const FromType& val, const CoCoA::ErrorContext& ErrCtx) { ToType copy; if (!IsConvertible(copy, val))  CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::BadConvert, ErrCtx);  return copy; }


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

#endif
