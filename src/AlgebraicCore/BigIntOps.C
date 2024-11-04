//   Copyright (c)  2012-2018  John Abbott and Anna Maria Bigatti

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

#include "CoCoA/BigIntOps.H"

#include "CoCoA/config.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/long32or64.H"
#include "CoCoA/utils.H"
#include "CoCoA/utils-gmp.H"

#include <cmath>
//using std::log;
//using std::atan;
#include <limits>
using std::numeric_limits;


namespace CoCoA
{

  namespace //anonymous
  {
    // Compute  round(n/d)  rounding halves up.
    // NB cannot simply do  (n+d/2)/d  because of possible overflow!
    inline unsigned long uround_half_up(unsigned long n, unsigned long d) noexcept
    {
      CoCoA_ASSERT(d != 0);
      const unsigned long q = n/d;
      if (n-q*d <= (d-1)/2) return q;
      return q+1; // cannot overflow
    }

    // 2022-02-07  Code is good but commented out to avoid compiler warning.
    // // Compute  round(n/d)  rounding halves down.
    // // NB cannot simply do  (n+(d-1)/2)/d  because of possible overflow!
    // inline unsigned long uround_half_down(unsigned long n, unsigned long d) noexcept
    // {
    //   CoCoA_ASSERT(d != 0);
    //   const unsigned long q = n/d;
    //   if (n-q*d <= d/2) return q;
    //   return q+1; // cannot overflow
    // }

  } // end of anonymous namespace


  BigInt abs(const BigInt& N)
  {
    BigInt ans;
    mpz_abs(mpzref(ans), mpzref(N));
    return ans;
  }


  BigInt operator-(const BigInt& N)
  {
    BigInt ans(N);
    negate(ans);
    return ans;
  }


  BigInt operator+(const BigInt& N1, const BigInt& N2)
  {
    BigInt ans;
    mpz_add(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  BigInt operator+(const BigInt& N1, const MachineInt& n2)
  {
    BigInt ans;
    if (IsNegative(n2))
      mpz_sub_ui(mpzref(ans), mpzref(N1), uabs(n2));
    else
      mpz_add_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  BigInt operator+(const MachineInt& n1, const BigInt& N2)
  {
    BigInt ans;
    if (IsNegative(n1))
      mpz_sub_ui(mpzref(ans), mpzref(N2), uabs(n1));
    else
      mpz_add_ui(mpzref(ans), mpzref(N2), AsUnsignedLong(n1));
    return ans;
  }


  BigInt operator-(const BigInt& N1, const BigInt& N2)
  {
    BigInt ans;
    mpz_sub(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  BigInt operator-(const BigInt& N1, const MachineInt& n2)
  {
    BigInt ans;
    if (IsNegative(n2))
      mpz_add_ui(mpzref(ans), mpzref(N1), uabs(n2));
    else
      mpz_sub_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  BigInt operator-(const MachineInt& n1, const BigInt& N2)
  {
    BigInt ans;
    if (IsNegative(n1))
      mpz_add_ui(mpzref(ans), mpzref(N2), uabs(n1));
    else
      mpz_sub_ui(mpzref(ans), mpzref(N2), AsUnsignedLong(n1));
    negate(ans);
    return ans;
  }


  BigInt operator*(const BigInt& N1, const BigInt& N2)
  {
    BigInt ans;
    // if (!mpz_mul_check(mpzref(ans), mpzref(N1), mpzref(N2)))
    // {
    //   CoCoA_THROW_ERROR(ERR::ArgTooBig, "overflow in big integer multiplication");
    // }
    mpz_mul(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  BigInt operator*(const BigInt& N1, const MachineInt& n2)
  {
    BigInt ans;
    if (IsNegative(n2))
      mpz_mul_si(mpzref(ans), mpzref(N1), AsSignedLong(n2));
    else
      mpz_mul_ui(mpzref(ans), mpzref(N1), AsUnsignedLong(n2));
    return ans;
  }

  BigInt operator*(const MachineInt& n1, const BigInt& N2)
  {
    BigInt ans;
    if (IsNegative(n1))
      mpz_mul_si(mpzref(ans), mpzref(N2), AsSignedLong(n1));
    else
      mpz_mul_ui(mpzref(ans), mpzref(N2), AsUnsignedLong(n1));
    return ans;
  }


  BigInt operator/(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigInt / BigInt");
    BigInt ans;
    mpz_tdiv_q(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  BigInt operator/(const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigInt / MachineInt");
    BigInt ans;
    mpz_tdiv_q_ui(mpzref(ans), mpzref(N1), uabs(n2));
    if (IsNegative(n2))
      negate(ans);

    return ans;
  }

  // Impl is simple rather than fast.
  BigInt operator/(const MachineInt& n1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "MachineInt / BigInt");
    BigInt ans(n1);
    ans /= N2;
    return ans;
  }


  BigInt operator%(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroModulus, "BigInt % BigInt");
    BigInt ans;
    mpz_tdiv_r(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  // This impl is valid for any second arg, but would then
  // have to return an unsigned long.  I have decided to forbid
  // 2nd args which are too big to fit into signed long; this
  // guarantees that the result will fit into signed long.
  // An alternative would be simply to IntegerCast the result to
  // long (which will throw if the result doesn't fit).
  long operator%(const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "BigInt % MachineInt");
    if (!IsSignedLong(n2))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "BigInt % MachineInt");
    const long ans = mpz_tdiv_ui(mpzref(N1), uabs(n2)); // abs value of remainder
    if (N1 < 0)
      return -ans; // CAREFUL: read GMP doc for mpz_tdiv_ui
    return ans;
  }

  // Impl is simple rather than fast.
  BigInt operator%(const MachineInt& n1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "MachineInt % BigInt");
    return BigInt(n1)%N2;
  }


  BigInt LeastNNegRemainder(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroModulus, "LeastNNegRemainder(BigInt, BigInt)");
    BigInt ans;
    mpz_mod(mpzref(ans), mpzref(N1), mpzref(N2));
    return ans;
  }

  long LeastNNegRemainder(const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2) || !IsSignedLong(n2)) CoCoA_THROW_ERROR(ERR::BadModulus, "LeastNNegRemainder");
    return mpz_fdiv_ui(mpzref(N1), uabs(n2)); // harmless silent conversion to long
  }

  BigInt LeastNNegRemainder(const MachineInt& n1, const BigInt& N2)
  {
    // Simple rather than fast;
    return LeastNNegRemainder(BigInt(n1), N2);
  }

  long LeastNNegRemainder(const MachineInt& n1, const MachineInt& n2)
  {
    if (IsZero(n2) || !IsSignedLong(n2)) CoCoA_THROW_ERROR(ERR::BadModulus, "LeastNNegRemainder");
    const unsigned long N2 = uabs(n2);
    const unsigned long ans = uabs(n1)%N2;
    if (ans != 0 && IsNegative(n1)) return N2-ans; // harmless silent conversion to long
    return ans; // harmless silent conversion to long
  }


  BigInt SymmRemainder(const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroModulus, "SymmRemainder(BigInt, BigInt)");
    BigInt ans;
    mpz_tdiv_r(mpzref(ans), mpzref(N1), mpzref(N2));
    BigInt HalfN2; mpz_tdiv_q_2exp(mpzref(HalfN2), mpzref(N2), 1); mpz_abs(mpzref(HalfN2), mpzref(HalfN2));
    if (CmpAbs(ans, HalfN2) < 0) return ans;
    if (ans == HalfN2) return ans;
    if (sign(ans) == sign(N2)) ans -= N2;
    else ans += N2;
    return ans;
  }

  long SymmRemainder(const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2) || !IsSignedLong(n2)) CoCoA_THROW_ERROR(ERR::BadModulus, "SymmRemainder");
    const unsigned long N2 = uabs(n2);
    unsigned long ans = mpz_fdiv_ui(mpzref(N1), N2);
    if (ans > N2/2) return -IntegerCast<long>(N2-ans); // cannot overflow (on a binary computer)
    return ans;
  }

  BigInt SymmRemainder(const MachineInt& n1, const BigInt& N2)
  {
    // Simple rather than fast.
    return SymmRemainder(BigInt(n1), N2);
  }


  long SymmRemainder(const MachineInt& n1, const MachineInt& n2)
  {
    if (IsZero(n2) || !IsSignedLong(n2)) CoCoA_THROW_ERROR(ERR::BadModulus, "SymmRemainder");
    const unsigned long N2 = uabs(n2);
    unsigned long ans = uabs(n1)%N2;
    if (IsNegative(n1)) ans = N2-ans;
    if (ans > N2/2) return -IntegerCast<long>(N2-ans); // cannot overflow (on a binary computer)
    return ans;
  }


  namespace // anonymous
  {

    void power_OverflowCheck(const BigInt& base, unsigned long exponent)
    {
      // Note: OVERFLOW_BITS defd in config.H
      if (IsZero(base) || IsOne(base) || IsMinusOne(base)) return;
      if (IsZero(exponent) || IsOne(exponent)) return;
      if (exponent > OVERFLOW_BITS)
        CoCoA_THROW_ERROR(ERR::ExpTooBig, "power");
      const long BitsBase = FloorLog2(base);
      if (BitsBase < (OVERFLOW_BITS/static_cast<long>(exponent))/2) return; // since exponent <= OVERFLOW_BITS, static_cast is safe
      const double LogBase = LogAbs(base);
      if (LogBase*exponent > OVERFLOW_BITS)      // values are "double"
        CoCoA_THROW_ERROR(ERR::ArgTooBig, "power");
    }
    
  }

  BigInt power(const BigInt& base, const BigInt& exponent, BigIntPowerOverflowCheck OverFlowCheck /*default: ENABLED*/)
  {
    if (exponent < 0)
      CoCoA_THROW_ERROR(ERR::NegExp, "power(BigInt, BigInt)");
    // trivial cases
    if (IsZero(exponent)) return BigInt(1);
    if (IsZero(base) || IsOne(base) || IsOne(exponent)) return base;
    if (IsMinusOne(base)) { if (IsEven(exponent)) return BigInt(1); else return BigInt(-1); }
    unsigned long exp; if (!IsConvertible(exp, exponent))
      CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(BigInt, BigInt)");
    if (OverFlowCheck == BigIntPowerOverflowCheck::ENABLED)  power_OverflowCheck(base, exp);
    BigInt ans;
    mpz_pow_ui(mpzref(ans), mpzref(base), exp);
    return ans;
  }

  BigInt power(const BigInt& base, const MachineInt& exponent, BigIntPowerOverflowCheck OverFlowCheck /*default: ENABLED*/)
  {
    if (IsNegative(exponent))
      CoCoA_THROW_ERROR(ERR::NegExp, "power(BigInt, MachineInt)");
    const unsigned long exp = AsUnsignedLong(exponent);
    // trivial cases
    if (exp == 0) return BigInt(1); // Recall that 0^0 is 1
    if (IsZero(base) || IsOne(base) || exp == 1) return base;
    if (IsMinusOne(base)) { if (IsEven(exp)) return BigInt(1); else return BigInt(-1); }
    if (OverFlowCheck == BigIntPowerOverflowCheck::ENABLED)  power_OverflowCheck(base, exp);
    BigInt ans;
    mpz_pow_ui(mpzref(ans), mpzref(base), exp);
    return ans;
  }

  BigInt power(const MachineInt& base, const BigInt& exponent, BigIntPowerOverflowCheck OverFlowCheck /*default: ENABLED*/)
  {
    const char* const FnName = "power(MachineInt, BigInt)";
    if (exponent < 0)
      CoCoA_THROW_ERROR(ERR::NegExp, FnName);
    // trivial cases
    if (IsZero(exponent)) return BigInt(1);
    if (IsZero(base) || IsOne(base) || IsOne(exponent)) return BigInt(base);
    if (IsMinusOne(base)) { if (IsEven(exponent)) return BigInt(1); else return BigInt(-1); }
    unsigned long exp; if (!IsConvertible(exp, exponent))
      CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(BigInt, BigInt)");
    if (OverFlowCheck == BigIntPowerOverflowCheck::ENABLED)  power_OverflowCheck(BigInt(base), exp);
    const bool NegateResult = (IsNegative(base) && IsOdd(exp));
    const unsigned long B = uabs(base);
    BigInt ans;
    mpz_ui_pow_ui(mpzref(ans), B, exp);
    if (NegateResult)
      negate(ans);
    return ans;
  }

  BigInt power(const MachineInt& base, const MachineInt& exponent, BigIntPowerOverflowCheck OverFlowCheck /*default: ENABLED*/)
  {
    const char* const FnName = "power(MachineInt, MachineInt)";
    if (IsNegative(exponent))
      CoCoA_THROW_ERROR(ERR::NegExp, FnName);
    const unsigned long exp = AsUnsignedLong(exponent);
    // trivial cases
    if (exp == 0) return BigInt(1);
    if (IsZero(base) || IsOne(base) || exp == 1) return BigInt(base);
    if (IsMinusOne(base)) { if (IsEven(exp)) return BigInt(1); else return BigInt(-1); }
    if (OverFlowCheck == BigIntPowerOverflowCheck::ENABLED)  power_OverflowCheck(BigInt(base), exp);

    const bool NegateResult = (IsNegative(base) && IsOdd(exp));
    const unsigned long B = uabs(base);
    BigInt ans;
    mpz_ui_pow_ui(mpzref(ans), B, exp);
    if (NegateResult)
      negate(ans);
    return ans;
  }


  long SmallPower(const MachineInt& base, const MachineInt& exponent)
  {
    if (IsNegative(exponent)) CoCoA_THROW_ERROR(ERR::NegExp, "SmallPower");
    if (!IsSignedLong(base)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "SmallPower: 1st arg");
    if (IsZero(exponent)) return 1; // Order is important for  0^0
    if (IsZero(base)) return 0;     // ...
    long b = AsSignedLong(base);
    if (b == 1) return 1;
    unsigned long e = AsUnsignedLong(exponent);
    if (b == -1) return (IsEven(e)?1:-1);
    const bool negative = (IsOdd(e) && IsNegative(base));
    long ans = 1;
    while (e != 1)
    {
      if (e&1) ans *= b;
      e >>= 1;
      b *= b;
    }
    ans *= b;
    if (negative) return -ans;
    return ans;
  }


  // --------------------------------------------
  // SumBigInt members

//   void SumBigInt::myAdd(const BigInt& N)
//   {
// //    std::scoped_lock lock(myMutex); // needs C++17
//     std::lock_guard<std::mutex> lock(myMutex); // cannot use auto :-(
//     this->operator+=(N);
//   }


  const BigInt& SumBigInt::myTotal() const
  {
//???   std::lock_guard<std::mutex> lock(myMutex);
    if (IsZero(myNegativeSum))  return myPositiveSum;
    if (IsZero(myPositiveSum))  return myNegativeSum;
    myPositiveSum += myNegativeSum;  myNegativeSum = 0;
    if (sign(myPositiveSum) >= 0)  return myPositiveSum;
    swap(myPositiveSum, myNegativeSum);
    return myNegativeSum;
  }


  std::ostream& operator<<(std::ostream& out, const SumBigInt& S)
  {
    if (!out)  return out;
    out << "SumBigInt(myNegativeSum=" << FloorLog2(S.myNegativeSum) << " bits, myPositiveSum=" << FloorLog2(S.myPositiveSum) << " bits)";
    return out;
  }


  ///////////////////////////////////////////////////////

  void add(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    mpz_add(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void add(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      mpz_sub_ui(mpzref(lhs), mpzref(N1), uabs(n2));
    else
      mpz_add_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  void add(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    add(lhs, N2, n1);
  }


  void sub(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    mpz_sub(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void sub(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      mpz_add_ui(mpzref(lhs), mpzref(N1), uabs(n2));
    else
      mpz_sub_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  void sub(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    sub(lhs, N2, n1);
    negate(lhs);
  }


  void mul(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    mpz_mul(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void mul(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsNegative(n2))
      mpz_mul_si(mpzref(lhs), mpzref(N1), AsSignedLong(n2));
    else
      mpz_mul_ui(mpzref(lhs), mpzref(N1), AsUnsignedLong(n2));
  }

  void mul(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    mul(lhs, N2, n1);
  }


  void div(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "div(BigInt, BigInt, BigInt)");
    mpz_tdiv_q(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void div(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_THROW_ERROR(ERR::DivByZero, "div(BigInt, BigInt, MachineInt)");
    mpz_tdiv_q_ui(mpzref(lhs), mpzref(N1), uabs(n2));
    if (IsNegative(n2))
      negate(lhs);
  }

  // Impl is simple rather than fast.
  void div(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    if (N2 == 0)
      CoCoA_THROW_ERROR(ERR::DivByZero, "div(BigInt, MachineInt, BigInt)");
    div(lhs, BigInt(n1), N2);
  }


  void mod(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroModulus, "mod(BigInt, BigInt, BigInt)");
    mpz_tdiv_r(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

  void mod(BigInt& lhs, const BigInt& N1, const MachineInt& n2)
  {
    if (IsZero(n2))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroModulus, "mod(BigInt, BigInt, MachineInt)");
    mpz_tdiv_r_ui(mpzref(lhs), mpzref(N1), uabs(n2));
  }

  // Impl is simple rather than fast.
  void mod(BigInt& lhs, const MachineInt& n1, const BigInt& N2)
  {
    if (IsZero(N2))
      CoCoA_THROW_ERROR(ERR::ReqNonZeroModulus, "mod(BigInt, MachineInt, BigInt)");
    mod(lhs, BigInt(n1), N2);
  }


  void quorem(BigInt& quo, BigInt& rem, const BigInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_THROW_ERROR(ERR::DivByZero, "quorem(BigInt, BigInt, BigInt, BigInt)");
    mpz_tdiv_qr(mpzref(quo), mpzref(rem), mpzref(num), mpzref(den));
  }

  void quorem(BigInt& quo, BigInt& rem, const BigInt& num, const MachineInt& den)
  {
    if (IsZero(den))
      CoCoA_THROW_ERROR(ERR::DivByZero, "quorem(BigInt, BigInt, BigInt, MachineInt)");
    if (IsNegative(den))
    {
      mpz_tdiv_qr_ui(mpzref(quo), mpzref(rem), mpzref(num), uabs(den));
      negate(quo);
    }
    else
      mpz_tdiv_qr_ui(mpzref(quo), mpzref(rem), mpzref(num), AsUnsignedLong(den));
  }

  // Impl is simple rather than fast.
  void quorem(BigInt& quo, BigInt& rem, const MachineInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_THROW_ERROR(ERR::DivByZero, "quorem(BigInt, BigInt, MachineInt, BigInt)");
    quorem(quo, rem, BigInt(num), den);
  }


  void divexact(BigInt& lhs, const BigInt& N1, const BigInt& N2)
  {
    CoCoA_ASSERT(N1%N2 == 0);
    mpz_divexact(mpzref(lhs), mpzref(N1), mpzref(N2));
  }

//???  void divexact(BigInt& lhs, const MachineInt& N1, const BigInt& N2);
//???  void divexact(BigInt& lhs, const BigInt& N1, const MachineInt& N2);


  bool IsZero(const BigInt& N) noexcept
  {
    return mpz_sgn(mpzref(N))==0;
  }

  bool IsOne(const BigInt& N) noexcept
  {
    return mpz_cmp_ui(mpzref(N), 1)==0;
  }

  bool IsMinusOne(const BigInt& N) noexcept
  {
    return mpz_cmp_si(mpzref(N),-1)==0;
  }


  bool IsPowerOf2(const BigInt& N) noexcept
  {
    if (N <= 0) return false;
    const std::size_t LogN = mpz_sizeinbase(mpzref(N), 2)-1;
    const std::size_t Last1Bit = mpz_scan1(mpzref(N), 0);
    return (Last1Bit == LogN);
  }

  bool IsPowerOf2(const MachineInt& n) noexcept
  {
    if (IsZero(n) || IsNegative(n)) return false;
    const unsigned long N = AsUnsignedLong(n);
    return !(N & (N-1));
  }


  bool IsDivisible(const MachineInt& N, const MachineInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_THROW_ERROR(ERR::DivByZero, "IsDivisible");
    return uabs(N)%uabs(D) == 0;
  }

  bool IsDivisible(const MachineInt& N, const BigInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_THROW_ERROR(ERR::DivByZero, "IsDivisible");
    // Simple rather than fast...
    return IsDivisible(BigInt(N), D);
  }

  bool IsDivisible(const BigInt& N, const MachineInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_THROW_ERROR(ERR::DivByZero, "IsDivisible");
    return mpz_divisible_ui_p(mpzref(N), uabs(D));
  }

  bool IsDivisible(const BigInt& N, const BigInt& D) // is N divisibile by D?
  {
    if (IsZero(D)) CoCoA_THROW_ERROR(ERR::DivByZero, "IsDivisible");
    return mpz_divisible_p(mpzref(N), mpzref(D));
  }


  int cmp(const BigInt& N1, const BigInt& N2) noexcept
  {
    return sign(mpz_cmp(mpzref(N1), mpzref(N2)));
  }

  int cmp(const BigInt& N1, const MachineInt& n2) noexcept
  {
    if (IsNegative(n2))
      return sign(mpz_cmp_si(mpzref(N1), AsSignedLong(n2)));
    else
      return sign(mpz_cmp_ui(mpzref(N1), AsUnsignedLong(n2)));
  }

  int cmp(const MachineInt& n1, const BigInt& N2) noexcept
  {
    return -cmp(N2, n1);
  }

  int cmp(const MachineInt& n1, const MachineInt& n2) noexcept
  {
    if (IsNegative(n1))
    {
      if (IsNegative(n2))
        return cmp(AsSignedLong(n1), AsSignedLong(n2));
      else
        return -1;
    }
    // n1 is non-negative
    if (IsNegative(n2)) return 1;
    // Both n1 and n2 are non-negative
    return cmp_ul(AsUnsignedLong(n1), AsUnsignedLong(n2));
  }


  int CmpAbs(const BigInt& N1, const BigInt& N2)  noexcept
  {
    return sign(mpz_cmpabs(mpzref(N1), mpzref(N2)));
  }

  int CmpAbs(const BigInt& N1, const MachineInt& n2) noexcept
  {
    return sign(mpz_cmpabs_ui(mpzref(N1), uabs(n2)));
  }

  int CmpAbs(const MachineInt& n1, const BigInt& N2) noexcept
  {
    return -sign(mpz_cmpabs_ui(mpzref(N2), uabs(n1)));
  }

  int CmpAbs(const MachineInt& n1, const MachineInt& n2) noexcept
  {
    const unsigned long a1 = uabs(n1);
    const unsigned long a2 = uabs(n2);
    if (a1 > a2) return 1;
    if (a1 < a2) return -1;
    return 0;
  }



  bool operator==(const BigInt& N1, const BigInt& N2) noexcept
  {
    return cmp(N1, N2) == 0;
  }

  bool operator==(const BigInt& N1, const MachineInt& n2) noexcept
  {
    return cmp(N1, n2) == 0;
  }

  bool operator==(const MachineInt& n1, const BigInt& N2) noexcept
  {
    return cmp(n1, N2) == 0;
  }


  bool operator!=(const BigInt& N1, const BigInt& N2) noexcept
  {
    return cmp(N1, N2) != 0;
  }

  bool operator!=(const BigInt& N1, const MachineInt& n2) noexcept
  {
    return cmp(N1, n2) != 0;
  }

  bool operator!=(const MachineInt& n1, const BigInt& N2) noexcept
  {
    return cmp(n1, N2) != 0;
  }


  bool operator> (const BigInt& N1, const BigInt& N2) noexcept
  {
    return cmp(N1, N2) > 0;
  }

  bool operator> (const BigInt& N1, const MachineInt& n2) noexcept
  {
    return cmp(N1, n2) > 0;
  }

  bool operator> (const MachineInt& n1, const BigInt& N2) noexcept
  {
    return cmp(n1, N2) > 0;
  }


  bool operator>=(const BigInt& N1, const BigInt& N2) noexcept
  {
    return cmp(N1, N2) >= 0;
  }

  bool operator>=(const BigInt& N1, const MachineInt& n2) noexcept
  {
    return cmp(N1, n2) >= 0;
  }

  bool operator>=(const MachineInt& n1, const BigInt& N2) noexcept
  {
    return cmp(n1, N2) >= 0;
  }


  bool operator< (const BigInt& N1, const BigInt& N2) noexcept
  {
    return cmp(N1, N2) < 0;
  }

  bool operator< (const BigInt& N1, const MachineInt& n2) noexcept
  {
    return cmp(N1, n2) < 0;
  }

  bool operator< (const MachineInt& n1, const BigInt& N2) noexcept
  {
    return cmp(n1, N2) < 0;
  }


  bool operator<=(const BigInt& N1, const BigInt& N2) noexcept
  {
    return cmp(N1, N2) <= 0;
  }

  bool operator<=(const BigInt& N1, const MachineInt& n2) noexcept
  {
    return cmp(N1, n2) <= 0;
  }
			
  bool operator<=(const MachineInt& n1, const BigInt& N2) noexcept
  {
    return cmp(n1, N2) <= 0;
  }
			

  double mantissa(long& exp, const BigInt& N) noexcept
  {
    return mpz_get_d_2exp(&exp, mpzref(N));
  }

  // double mantissa(const BigInt& N) noexcept
  // {
  //   long DiscardExponent;
  //   return mpz_get_d_2exp(&DiscardExponent, mpzref(N));
  // }


  long exponent(const BigInt& N) noexcept
  {
    return 1+FloorLog2(N); // BUG?? overflow???
  }


  long MultiplicityOf2(const BigInt& N)
  {
    if (IsZero(N)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "MultiplicityOf2");
    return mpz_scan1(mpzref(N), 0);
  }

  long MultiplicityOf2(const MachineInt& n)
  {
    if (IsZero(n)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "MultiplicityOf2");
    unsigned long val = uabs(n);
/////    unsigned long long val = ull_abs(n);
    int i=0;
    while ((val&1) == 0)  { val /= 2; ++i; }
    return i;
  }


  double log(const BigInt& N)
  {
    if (sign(N) <= 0)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "log(BigInt)");
// Next 2 lines useful???  Extra code for no benefit?
///    if (exponent(N) < numeric_limits<double>::max_exponent)
///      return std::log(ConvertTo<double>(N));
    long exp;
    const double mantissa = mpz_get_d_2exp(&exp, mpzref(N));
    const double log2 = std::log(2.0);
    return std::log(std::abs(mantissa)) + exp*log2;
  }


  double LogAbs(const BigInt& N)
  {
    const int s = sign(N);
    if (s == 0) CoCoA_THROW_ERROR(ERR::ReqNonZero, "LogAbs(BigInt)");
    if (s > 0) return log(N);
    return log(-N);
  }

  namespace // anonymous
  {
    // Ancient code from 2006 or earlier (found in a file called topbit.C)
    
    // Goes about twice as fast if inlined.  Goes about 25% faster with table of 256 entries.
    inline size_t topbit(unsigned long n)
    {
      static const unsigned char topbittab[16] = {0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3};
//   static const unsigned char topbittab[256] = {0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
//                                                5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
//                                                6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
//                                                6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
//                                                7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
//                                                7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
//                                                7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
//                                                7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7};
//FASTER???      if (n <= 0xf) return topbittab[n];
      unsigned char pos = 0;
#if !defined(CoCoA_32BIT_LONG)
      if (n > 0xffffffff) { n >>= 32; pos += 32; }
#endif
      if (n > 0x0000ffff) { n >>= 16; pos += 16; }
      if (n > 0x000000ff) { n >>= 8;  pos += 8; }
      if (n > 0x0000000f) { n >>= 4;  pos += 4; }  // remove this line if you use the 256 size table
      return pos + topbittab[n];
    }

  } // end of namespace anonymous

  long FloorLog2(const MachineInt& n)
  {
    if (IsZero(n))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "FloorLog2");
    return topbit(uabs(n));
  }

  long FloorLog2(const BigInt& N) /*noexcept*/
  {
    if (IsZero(N))
      CoCoA_THROW_ERROR(ERR::ReqNonZero, "FloorLog2");
    const size_t sz = mpz_sizeinbase(mpzref(N),2)-1;
    // CHECK FOR OVERFLOW  - only necessary if long is 32-bit (M$ users deserve what they get?)
    // if (sz > static_cast<size_t>(numeric_limits<long>::max()))
    //   CoCoA_THROW_ERROR(ERR::ArgTooBig, "FloorLog2");
    return static_cast<long>(sz);
  }

  long FloorLog10(const MachineInt& n) { return FloorLogBase(n,10); }
  long FloorLog10(const BigInt& N) { return FloorLogBase(N,10); }


  long FloorLogBase(const MachineInt& n, const MachineInt& base)
  {
    return FloorLogBase(BigInt(n), BigInt(base));
  }

  long FloorLogBase(const MachineInt& n, const BigInt& base)
  {
    return FloorLogBase(BigInt(n), base);
  }

  long FloorLogBase(const BigInt& N, const MachineInt& base)
  {
    return FloorLogBase(N, BigInt(base));
  }

  long FloorLogBase(const BigInt& N, const BigInt& base)
  {
    if (base < 2) CoCoA_THROW_ERROR(ERR::BadArg, "FloorLogBase: base must be at least 2");
    if (IsZero(N)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "FloorLogBase");
    if (IsOne(abs(N))) return 0;
    if (base == 2) return FloorLog2(N);
    const double ApproxLog = LogAbs(N)/log(base);
    // Check whether ApproxLog is far from integer.
    const double delta = 5 * ApproxLog * numeric_limits<double>::epsilon(); // estimate max poss numerical error in ApproxLog
    const long candidate = static_cast<long>(floor(ApproxLog+delta)); // probably right but could be too big by 1
    if (std::abs(ApproxLog - candidate) > delta)
    {
      // ApproxLog was not too close to an integer, so we are sure that candidate is correct.
#ifdef CoCoA_DEBUG
      {
        // LINE BELOW IS BigInt pwr = power(base, candidate);  BUT AVOIDS OVERFLOW CHECK
        BigInt pwr; mpz_pow_ui(mpzref(pwr), mpzref(base), candidate); // silent conversion of candidate to ulong is safe
        CoCoA_ASSERT(CmpAbs(N,pwr) >= 0);
      }
#endif
      return candidate;
    }
    // ApproxLog was close to integer, so we make a slow-but-sure check.
    // LINE BELOW IS const BigInt pwr = power(base, candidate);  BUT AVOIDS OVERFLOW CHECK
    BigInt pwr; mpz_pow_ui(mpzref(pwr), mpzref(base), candidate); // silent conversion of candidate to ulong is safe
    const int test = CmpAbs(N, pwr);
    if (test >= 0)  return candidate;
    return candidate-1;
  }


  // This fn comes out a long horrible mess because I cannot return a
  // machine int as a BigInt.  There must be a design error somewhere.
  BigInt RangeFactorial(const MachineInt& lo, const MachineInt& hi)
  {
    BigInt ans;
    // If zero is in the range, result is zero
    if ((IsZero(lo) || IsNegative(lo)) && !IsNegative(hi)) return ans;
    ans = 1;
    // First some checks for empty products.
    if (IsNegative(hi) && !IsNegative(lo)) return ans;
    // We now know that hi and lo have the same sign.
    // Two more checks for empty products.
    if (IsNegative(lo) && AsSignedLong(lo) > AsSignedLong(hi)) return ans;
    if (!IsNegative(lo) && AsUnsignedLong(lo) > AsUnsignedLong(hi)) return ans;
    // We now know that the product is not empty & does not span 0.
    unsigned long l = uabs(lo);
    unsigned long h = uabs(hi);
    if (IsNegative(hi))
    {
      std::swap(l, h);
      if (((h^l)&1) == 0) negate(ans);
    }

    if (h-l > 15)
    {
      const unsigned long mid = l + (h-l)/2; // equiv to (l+h)/2 but avoids overflow
      // const BigInt FirstHalf = RangeFactorial(l,mid);
      // CheckForInterrupt("RangeFactorial");
      // const BigInt SecondHalf = RangeFactorial(1+mid,h);
      // CheckForInterrupt("RangeFactorial");
      // return FirstHalf*SecondHalf;
      return RangeFactorial(l,mid) * RangeFactorial(1+mid,h);
    }
    // From here on 0 < l < h.
    for (; l <= h; ++l)
    {
      ans *= l;
//      mpz_mul_ui(mpzref(ans), mpzref(ans), l);
    }
    return ans;
  }


  namespace // anonymous
  {
    void factorial_OverflowCheck(long n)
    {
      CoCoA_ASSERT(n >= 0);
      if (n < 10000000) return; // assume no overflow up to n=1000000
      /// SLUG BUG SLUG should just check n directly!!!
      const double log_est = LogFactorial(n)/std::log(2.0);
      if (log_est < OVERFLOW_BITS) return; // OVERFLOW_BITS defn in config.H
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "factorial");
    }
    
  }

  BigInt factorial(const BigInt& N)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "factorial(BigInt)");
    long n;
    if (!IsConvertible(n, N))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "factorial(BigInt)");
    return factorial(n);
    // factorial_OverflowCheck(n);
    // BigInt ans;
    // mpz_fac_ui(mpzref(ans), n);
    // return ans;
  }


  BigInt factorial(const MachineInt& n)
  {
    if (IsNegative(n))
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "factorial(MachineInt)");

    factorial_OverflowCheck(AsSignedLong(n));
    BigInt ans;
    mpz_fac_ui(mpzref(ans), AsUnsignedLong(n));
    return ans;
  }


// From Wikipedia; apparently due to Ramanujan
// Experimentally max error is 4.6*10^(-8) at n=24.
double LogFactorial(const MachineInt& N)
{
  using std::log;
  using std::atan;
  if (IsNegative(N)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "LogFactorial");
  const long n = AsSignedLong(N);
  const static double LogSqrtPi = log(4*atan(1.0))/2;
  if (n > 23)
    return n*log(double(n))-n + log(n*(1+n*(4.0+8*n)))/6 + LogSqrtPi;
  double fact = 1.0;
  for (long k=2; k <= n; ++k)
    fact *= k;
  return log(fact);
}

double LogFactorial(const BigInt& N)
{
  if (N < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "LogFactorial");
  double n;
  if (!IsConvertible(n, N)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "LogFactorial");
  {
    long SmallN;
    if (IsConvertible(SmallN, N)) return LogFactorial(SmallN);
  }
  using std::log;
  using std::atan;
  const static double LogSqrtPi = log(4*atan(1.0))/2;
  return n*log(n)-n + log(n*(1+n*(4.0+8*n)))/6 + LogSqrtPi;
}



  BigInt primorial(const BigInt& N)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "primorial(BigInt)");
    unsigned long n;
    if (!IsConvertible(n, N))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "primorial(BigInt)");
//    return factorial(n);
    BigInt ans;
    mpz_primorial_ui(mpzref(ans), n);
    return ans;
  }


  BigInt primorial(const MachineInt& n)
  {
    if (IsNegative(n))
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "factorial(MachineInt)");
    BigInt ans;
    mpz_primorial_ui(mpzref(ans), AsUnsignedLong(n));
    return ans;
  }


namespace // anonymous
{
  void binomial_OverflowCheck(const BigInt& N, long r)
  {
    constexpr long OVERFLOW_BITS = 1000000000;
    CoCoA_ASSERT(r >= 0 && r <= N);
    if (N < OVERFLOW_BITS) return; // N small enough that overflow cannot occur
    if (r == 0 || r == 1) return;
    if (r*log(N) < OVERFLOW_BITS) return;
    double log_num = LogFactorial(N);
    double log_den1 = LogFactorial(r);
    double log_den2 = LogFactorial(N-r);
    if (log_num - (log_den1+log_den2) < OVERFLOW_BITS) return;
    CoCoA_THROW_ERROR(ERR::ArgTooBig, "binomial");
  }
  
} // end of namespace anonymous

  BigInt binomial(const BigInt& N, const BigInt& R)
  {
    if (R < 0 || CmpAbs(R,N) > 0) return BigInt(0);
//    if (R < 0)
//      CoCoA_THROW_ERROR(ERR::BadArg, "binomial(BigInt,BigInt): 2nd arg must be non-neg");
    if (N < 0)
    {
      if (IsEven(R)) return binomial(R-N-1,R);
      else return -binomial(R-N-1,R);
    }
    if (R > N)
      return BigInt(0);

    // General case:  N >= R >= 0
    unsigned long r;
    const bool RIsSmall = (2*R < N);
    if ((RIsSmall && !IsConvertible(r, R)) ||
        (!RIsSmall && !IsConvertible(r, N-R)))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "binomial(BigInt, BigInt)");

    binomial_OverflowCheck(N,r);
    BigInt ans;
    mpz_bin_ui(mpzref(ans), mpzref(N), r); // we know that r >= 0.
    return ans;
  }


  BigInt binomial(const MachineInt& n, const BigInt& R)
  {
    if (R < 0 || CmpAbs(R,n) > 0) return BigInt(0);
//    if (R < 0)
//      CoCoA_THROW_ERROR(ERR::BadArg, "binomial(MachineInt,BigInt): 2nd arg must be non-neg");
    return binomial(BigInt(n), R);
  }


  BigInt binomial(const BigInt& N, const MachineInt& r)
  {
    if (IsNegative(r) || CmpAbs(r,N) > 0) return BigInt(0);
//    if (IsNegative(r))
//      CoCoA_THROW_ERROR(ERR::BadArg, "binomial(BigInt,MachineInt): 2nd arg must be non-neg");
    binomial_OverflowCheck(abs(N),AsSignedLong(r));
    BigInt ans;
    mpz_bin_ui(mpzref(ans), mpzref(N), AsUnsignedLong(r));
    return ans;
  }


  BigInt binomial(const MachineInt& n, const MachineInt& r)
  {
    if (IsNegative(r) || uabs(r) > uabs(n)) return BigInt(0);
//    if (IsNegative(r))
//      CoCoA_THROW_ERROR(ERR::BadArg, "binomial(MachineInt,MachineInt): 2nd arg must be non-neg");

    if (uabs(n) > 1000000) binomial_OverflowCheck(BigInt(uabs(n)),AsSignedLong(r));
    BigInt ans;
    if (IsNegative(n))
      mpz_bin_ui(mpzref(ans), mpzref(BigInt(n)), AsUnsignedLong(r));
    else
      mpz_bin_uiui(mpzref(ans), AsUnsignedLong(n), AsUnsignedLong(r));
    return ans;
  }


  BigInt fibonacci(const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "fibonacci(BigInt)");
    return fibonacci(n);
  }


  BigInt fibonacci(const MachineInt& n)
  {
    const bool EvenNegative = IsNegative(n) && IsEven(n);
    BigInt ans;
    unsigned long arg = uabs(n);
    mpz_fib_ui(mpzref(ans), arg);
    if (EvenNegative)
      negate(ans);
    return ans;
  }


  // NB halves are rounded away from zero.
  BigInt RoundDiv(const BigInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_THROW_ERROR(ERR::DivByZero, "RoundDiv(BigInt,BigInt)");
    BigInt q,r;
    quorem(q, r, num, den);
    if (CmpAbs(2*r,den) >= 0) // ROUND AWAY FROM ZERO
//    if (CmpAbs(2*r,den) > 0) // ROUND TOWARDS ZERO
      q += sign(r)*sign(den);
    return q;
  }

  BigInt RoundDiv(const BigInt& num, const MachineInt& den)
  {
    if (IsZero(den))
      CoCoA_THROW_ERROR(ERR::DivByZero, "RoundDiv(BigInt,MachInt)");
    return RoundDiv(num, BigInt(den));
  }

  BigInt RoundDiv(const MachineInt& num, const BigInt& den)
  {
    if (IsZero(den))
      CoCoA_THROW_ERROR(ERR::DivByZero, "RoundDiv(MachInt,BigInt)");
    return RoundDiv(BigInt(num), den);
  }

  // NB Halves are rounded away from zero.
  long RoundDiv(const MachineInt& n, const MachineInt& d)
  {
    if (IsZero(d))
      CoCoA_THROW_ERROR(ERR::DivByZero, "RoundDiv(n,d)");
    if (IsZero(n)) return 0;
    // Henceforth n and d are both non-zero.
    const long s = sign(n)*sign(d);
    const unsigned long q = uround_half_up(uabs(n), uabs(d)); // ROUND AWAY FROM ZERO
//    const unsigned long q = uround_half_down(uabs(n), uabs(d)); // ROUND TOWARDS ZERO
    if (s > 0)  return IntegerCast<long>(q);
    static const long MinLong = numeric_limits<long>::min();
    if (q + MinLong == 0) return MinLong; // to avoid overflow if result is MinLong.
    return -IntegerCast<long>(q);

  }


  BigInt FloorRoot(const MachineInt& n, const MachineInt& r)
  {
    return FloorRoot(BigInt(n), r);
  }

  BigInt FloorRoot(const MachineInt& n, const BigInt& R)
  {
    if (IsNegative(n))
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "FloorRoot: 1st arg");
    if (R < 1)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "FloorRoot: 2nd arg");
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "FloorRoot: 2nd arg"); // could just return 1???
    return FloorRoot(BigInt(n), r);
  }

  BigInt FloorRoot(const BigInt& N, const MachineInt& r)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "FloorRoot: 1st arg");
    if (IsNegative(r) || IsZero(r))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "FloorRoot: 2nd arg");
    BigInt ans;
    mpz_root(mpzref(ans), mpzref(N), AsUnsignedLong(r));
    return ans;
  }

  BigInt FloorRoot(const BigInt& N, const BigInt& R)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "FloorRoot: 1st arg");
    if (R < 1)
      CoCoA_THROW_ERROR(ERR::ReqPositive, "FloorRoot: 2nd arg");
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "FloorRoot: 2nd arg"); // could just return 1???
    return FloorRoot(N, r);
  }


  bool IsExactFloorRoot(long& ans, const MachineInt& n, const MachineInt& r)
  {
    BigInt root;
    const bool IsExact = IsExactFloorRoot(root, BigInt(n), r);
    if (!IsConvertible(ans, root))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsExactFloorRoot");
    return IsExact;
  }

  bool IsExactFloorRoot(BigInt& ans, const MachineInt& n, const MachineInt& r)
  {
    return IsExactFloorRoot(ans, BigInt(n), r);
  }

  bool IsExactFloorRoot(long& ans, const MachineInt& n, const BigInt& R)
  {
//    if (IsNegative(n))
//      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "IsExactFloorRoot: 1st arg must be non-negative");
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsExactFloorRoot: root index must be >=1 & fit into a ulong");
    BigInt root;
    const bool IsExact = IsExactFloorRoot(root, BigInt(n), r);
    if (!IsConvertible(ans, root))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsExactFloorRoot");
    return IsExact;
  }

  bool IsExactFloorRoot(BigInt& ans, const MachineInt& n, const BigInt& R)
  {
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsExactFloorRoot: root index must be >=1 & fit into a ulong");
    return IsExactFloorRoot(ans, BigInt(n), r);
  }

  bool IsExactFloorRoot(BigInt& ans, const BigInt& N, const MachineInt& r)
  {
    if (N < 0)
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "IsExactFloorRoot: 2nd arg");
    if (IsNegative(r) || IsZero(r))
      CoCoA_THROW_ERROR(ERR::ReqPositive, "IsExactFloorRoot: 3rd arg");
//    if (N < 0 && IsEven(AsUnsignedLong(r)))
//      CoCoA_THROW_ERROR(ERR::BadArg, "IsExactFloorRoot: even root of negative number");
    const bool IsExact = mpz_root(mpzref(ans), mpzref(N), AsUnsignedLong(r));
    return IsExact;
  }

  bool IsExactFloorRoot(BigInt& ans, const BigInt& N, const BigInt& R)
  {
    unsigned long r;
    if (!IsConvertible(r, R))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsExactFloorRoot: root index must be >=1 & fit into a ulong");
    return IsExactFloorRoot(ans, N, r);
  }


  long FloorSqrt(const MachineInt& n)
  {
    if (IsNegative(n)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "FloorSqrt(n)");
    return ConvertTo<long>(FloorSqrt(BigInt(n))); // conversion will succeed.
    // long ans;
    // IsConvertible(ans, FloorSqrt(BigInt(n))); // conversion will succeed.
    // return ans;
  }

  BigInt FloorSqrt(const BigInt& N)
  {
    if (N < 0) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "FloorSqrt(N)");
    BigInt ans;
    mpz_sqrt(mpzref(ans), mpzref(N));
    return ans;
  }


  bool IsSquare(const MachineInt& n)
  {
    return IsSquare(BigInt(n));
  }

  bool IsSquare(const BigInt& N) noexcept
  {
    return mpz_perfect_square_p(mpzref(N));
  }


  bool IsPower(const MachineInt& n)
  {
    return IsPower(BigInt(n));
  }

  bool IsPower(const BigInt& N) noexcept
  {
    return mpz_perfect_power_p(mpzref(N));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/BigIntOps.C,v 1.29 2024/07/02 15:29:22 bigatti Exp $
// $Log: BigIntOps.C,v $
// Revision 1.29  2024/07/02 15:29:22  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
// Specialized some BadArg with explicit codes
// Tidied up some CoCoA_THROW_ERROR second argument (string)
//
// Revision 1.28  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.27  2024/05/03 12:12:19  abbott
// Summary: Removed "long long" (redmine 1804)
//
// Revision 1.26  2024/04/12 16:55:52  abbott
// Summary: Changed name of CPP symbol which signifies that longs are 32 bits (redmine 1661/1804)
//
// Revision 1.25  2024/03/25 08:04:26  abbott
// Summary: Modifications to work  on MinGW (& 32 bitters) (redmine 1661)
//
// Revision 1.24  2024/03/15 18:53:59  abbott
// Summary: Corrected operator% for BigInt
//
// Revision 1.23  2024/01/25 21:06:00  abbott
// Summary: Removed thread-safety from SumBigInt/SumBigRat (costly/complicated to do it right)
//
// Revision 1.22  2023/12/09 21:37:32  abbott
// Summary: Tidied SumBigInt/SumBigRat; added doc
//
// Revision 1.21  2023/12/01 20:25:28  abbott
// Summary: Added flag to disable overflow check in power (for integers)
//
// Revision 1.20  2023/11/11 20:29:53  abbott
// Summary: Cleaned FloorLogBase; now avoids overflow check when debugging is on
//
// Revision 1.19  2023/03/09 22:34:18  abbott
// Summary: Commented out unused param -- to keep compiler quiet
//
// Revision 1.18  2022/11/16 20:03:58  abbott
// Summary: Two minor changes (hints from clang compiler)
//
// Revision 1.17  2022/10/25 17:23:49  abbott
// Summary: Added myAdd; threadsafe version of op+=; need mutex now
//
// Revision 1.16  2022/10/10 09:52:44  abbott
// Summary: Added SumBigInt class
//
// Revision 1.15  2022/08/08 17:53:45  abbott
// Summary: SizeInBase is obsolete (redmine 1693)
//
// Revision 1.14  2022/05/14 15:57:44  abbott
// Summary: FloorLog2 now calls mpz_sizeinbase (redmine 1680)
//
// Revision 1.13  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.12  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.11  2021/09/22 13:56:27  abbott
// Summary: Overflow check (power & factorial) (redmin 1585)
//
// Revision 1.10  2021/08/03 12:27:33  abbott
// Summary: Removed "useless" const qualifier on rtn values (redmine 1606)
//
// Revision 1.9  2021/08/03 08:34:49  abbott
// Summary: Revised impl of binomial (redmine 1605)
//
// Revision 1.8  2021/04/01 20:02:48  abbott
// Summary: Added overflow check to power fn (redmine 1585)
//
// Revision 1.7  2021/03/04 16:59:37  abbott
// Summary: log is not constexpr
//
// Revision 1.6  2021/02/10 19:39:59  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.5  2021/01/31 10:04:51  abbott
// Summary: Added noexcept; added LogAbs; variation to mantissa (not yet complete)
//
// Revision 1.4  2020/10/06 19:14:30  abbott
// Summary: SmallPower now handles some simple case specially (incl power of -1)
//
// Revision 1.3  2020/06/17 15:49:21  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.2  2019/09/16 14:30:34  abbott
// Summary: Changed iroot -> FloorRoot, IsEactIroot -> IsExactFloorRoot; added primorial
//
// Revision 1.1  2018/05/18 12:09:25  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.20  2017/11/29 20:30:17  abbott
// Summary: Added MultiplicityOf2C
//
// Revision 1.19  2017/10/16 14:48:05  abbott
// Summary: Corrected internal conversion to long (was unsigned long) in fibonacci
//
// Revision 1.18  2017/02/08 17:01:26  abbott
// Summary: Made uround_half_up and uround-half-down inline
//
// Revision 1.17  2016/08/02 12:49:34  abbott
// Summary: Renamed NumDigits to SizeInBase; updated doc.
//
// Revision 1.16  2015/11/24 12:47:35  abbott
// Summary: Renamed "isqrt" --> "FloorSqrt"
//
// Revision 1.15  2015/11/23 18:21:52  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.14  2015/11/05 14:20:06  abbott
// Summary: Changed rtn type of LeastNonNegRemainder to long when modulus is MachineInt
//
// Revision 1.13  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.12  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.11  2014/05/16 14:15:29  abbott
// Summary: Rewrote DoundDiv for MachineInts to handle rounding away from zero correctly when result is MinLong
// Author: JAA
//
// Revision 1.10  2014/05/16 12:26:11  abbott
// Summary: Changed rounding (in RoundDiv): it is now "away from zero"
// Author: JAA
//
// Revision 1.9  2014/05/05 18:21:00  abbott
// Removed another M_PI (replaced by 4*atan(1.0))
//
// Revision 1.8  2014/04/30 16:08:42  abbott
// Summary: Removed dependency on nonstd M_PI; now uses atan; cast arg of log to double
// Author: JAA
//
// Revision 1.7  2014/04/10 15:30:33  abbott
// Summary: Corrected impl of IsPowerOf2
// Author: JAA
//
// Revision 1.6  2014/04/08 13:04:13  abbott
// Summary: Added new fn IsPowerOf2
// Author: JAA
//
// Revision 1.5  2014/03/06 15:49:30  abbott
// Summary: Zero to power zero now gives 1
// Author: JAA
//
// Revision 1.4  2013/05/20 15:48:13  abbott
// Added new fn LogFactorial (placed in IntOperations).
//
// Revision 1.3  2013/05/14 14:23:00  abbott
// Minor improvement to log(BigInt)
//
// Revision 1.2  2012/12/04 09:55:47  abbott
// Added new fns LeastNNegRemainder and SymmRemainder (with doc & new test).
// Some minor corrections to the doc (for operator/ and operator%).
//
// Revision 1.1  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
//
