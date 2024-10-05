#ifndef CoCoA_BigIntOps_H
#define CoCoA_BigIntOps_H

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


#include "CoCoA/BigInt.H"
#include "CoCoA/MachineInt.H"

#include <mutex>

namespace CoCoA
{

  // STANDARD ARITHMETIC OPERATIONS

  int sign(const BigInt& N) noexcept;
  int sign(const MachineInt& n) noexcept;
  BigInt abs(const BigInt& N);
  // unsigned long abs(const MachineInt& n); // see MachineInt.H
  BigInt operator-(const BigInt& N);
  void negate(BigInt& N) noexcept;

  BigInt operator+(const BigInt& N1, const BigInt& N2);
  BigInt operator+(const BigInt& N1, const MachineInt& n2);
  BigInt operator+(const MachineInt& n1, const BigInt& N2);

  BigInt operator-(const BigInt& N1, const BigInt& N2);
  BigInt operator-(const BigInt& N1, const MachineInt& n2);
  BigInt operator-(const MachineInt& n1, const BigInt& N2);

  BigInt operator*(const BigInt& N1, const BigInt& N2);
  BigInt operator*(const BigInt& N1, const MachineInt& n2);
  BigInt operator*(const MachineInt& n1, const BigInt& N2);

  BigInt operator/(const BigInt& N1, const BigInt& N2);
  BigInt operator/(const BigInt& N1, const MachineInt& n2);
  BigInt operator/(const MachineInt& n1, const BigInt& N2);

  BigInt operator%(const BigInt& N1, const BigInt& N2);
  long operator%(const BigInt& N1, const MachineInt& n2);
  BigInt operator%(const MachineInt& n1, const BigInt& N2);
  BigInt LeastNNegRemainder(const BigInt& N1, const BigInt& N2);
  long LeastNNegRemainder(const BigInt& N1, const MachineInt& n2);
  BigInt LeastNNegRemainder(const MachineInt& n1, const BigInt& N2);
  long LeastNNegRemainder(const MachineInt& n1, const MachineInt& n2);
  BigInt SymmRemainder(const BigInt& N1, const BigInt& N2);
  long SymmRemainder(const BigInt& N1, const MachineInt& n2);
  BigInt SymmRemainder(const MachineInt& n1, const BigInt& N2);
  long SymmRemainder(const MachineInt& n1, const MachineInt& n2);

  enum class BigIntPowerOverflowCheck { DISABLED=0, ENABLED=1 };

  BigInt power(const BigInt& base, const BigInt& exponent, BigIntPowerOverflowCheck OverFlowCheck=BigIntPowerOverflowCheck::ENABLED);
  BigInt power(const BigInt& base, const MachineInt& exponent, BigIntPowerOverflowCheck OverFlowCheck=BigIntPowerOverflowCheck::ENABLED);
  BigInt power(const MachineInt& base, const BigInt& exponent, BigIntPowerOverflowCheck OverFlowCheck=BigIntPowerOverflowCheck::ENABLED);
  BigInt power(const MachineInt& base, const MachineInt& exponent, BigIntPowerOverflowCheck OverFlowCheck=BigIntPowerOverflowCheck::ENABLED);
  long SmallPower(const MachineInt& base, const MachineInt& exponent); ///< base^exponent (assuming it fits into a long)


  class SumBigInt
  {
  public:
    SumBigInt() = default;
    SumBigInt(const SumBigInt&) = default;
    SumBigInt(SumBigInt&&) = default;
    SumBigInt& operator=(const SumBigInt&) = default;
    SumBigInt& operator=(SumBigInt&&) = default;
    ~SumBigInt() = default;
    /*inline*/ SumBigInt& operator+=(long n)          // NOT THREADSAFE!
      { if (sign(n) >= 0)
          myPositiveSum += n;
        else
          myNegativeSum += n;
        return *this;
      }
    /*inline*/ SumBigInt& operator+=(const BigInt& N) // NOT THREADSAFE!
      { if (sign(N) >= 0)
          myPositiveSum += N;
        else
          myNegativeSum += N;
        return *this;
      }
///    void myAdd(const BigInt& N); // threadsafe version of op+=
    const BigInt& myTotal() const;
  private: // data members
///    std::mutex myMutex; // causes too much complication :-(
    mutable BigInt myPositiveSum;
    mutable BigInt myNegativeSum;
    friend std::ostream& operator<<(std::ostream& out, const SumBigInt& S);
  };

  inline const BigInt& total(const SumBigInt& S)  { return S.myTotal(); }
  std::ostream& operator<<(std::ostream& out, const SumBigInt& S);


  // Efficient arithmetic procedures...

  void add(BigInt& lhs, const BigInt& N1, const BigInt& N2);
  void add(BigInt& lhs, const BigInt& N1, const MachineInt& n2);
  void add(BigInt& lhs, const MachineInt& n1, const BigInt& N2);

  void sub(BigInt& lhs, const BigInt& N1, const BigInt& N2);
  void sub(BigInt& lhs, const BigInt& N1, const MachineInt& n2);
  void sub(BigInt& lhs, const MachineInt& n1, const BigInt& N2);

  void mul(BigInt& lhs, const BigInt& N1, const BigInt& N2);
  void mul(BigInt& lhs, const BigInt& N1, const MachineInt& n2);
  void mul(BigInt& lhs, const MachineInt& n1, const BigInt& N2);

  void div(BigInt& lhs, const BigInt& N1, const BigInt& N2);
  void div(BigInt& lhs, const BigInt& N1, const MachineInt& n2);
  void div(BigInt& lhs, const MachineInt& n1, const BigInt& N2);

  void mod(BigInt& lhs, const BigInt& N1, const BigInt& N2);
  void mod(BigInt& lhs, const BigInt& N1, const MachineInt& n2);
  void mod(BigInt& lhs, const MachineInt& n1, const BigInt& N2);

  void quorem(BigInt& quo, BigInt& rem, const BigInt& num, const BigInt& den);
  void quorem(BigInt& quo, BigInt& rem, const BigInt& num, const MachineInt& den); // why is rem not of type MachineInt???
  void quorem(BigInt& quo, BigInt& rem, const MachineInt& num, const BigInt& den);

  void divexact(BigInt& lhs, const BigInt& N1, const BigInt& N2);
//???  void divexact(BigInt& lhs, const MachineInt& N1, const BigInt& N2);
//???  void divexact(BigInt& lhs, const BigInt& N1, const MachineInt& N2);


  // SUNDRY PREDICATES.

  bool IsZero(const BigInt& N) noexcept;
  bool IsOne(const BigInt& N) noexcept;
  bool IsMinusOne(const BigInt& N) noexcept;

  bool IsOdd(const BigInt& N) noexcept;
  bool IsOdd(const MachineInt& n) noexcept;
  bool IsEven(const BigInt& N) noexcept;
  bool IsEven(const MachineInt& n) noexcept;
  bool IsPowerOf2(const BigInt& N) noexcept;
  bool IsPowerOf2(const MachineInt& n) noexcept;

  // IsDivisible will throw ERR::DivByZero if D==0
  bool IsDivisible(const MachineInt& N, const MachineInt& D); // is N divisibile by D?
  bool IsDivisible(const MachineInt& N, const BigInt& D);     // is N divisibile by D?
  bool IsDivisible(const BigInt& N, const MachineInt& D);     // is N divisibile by D?
  bool IsDivisible(const BigInt& N, const BigInt& D);         // is N divisibile by D?
//bool IsPerfectSquare(const BigInt& N);
// bool IsPerfectSquare(const MachineInt& n);


  // COMPARISON FUNCTIONS

  int cmp(const BigInt& N1, const BigInt& N2) noexcept;
  int cmp(const BigInt& N1, const MachineInt& n2) noexcept;
  int cmp(const MachineInt& n1, const BigInt& N2) noexcept;
  int cmp(const MachineInt& n1, const MachineInt& n2) noexcept;

  int CmpAbs(const BigInt& N1, const BigInt& N2) noexcept;
  int CmpAbs(const BigInt& N1, const MachineInt& n2) noexcept;
  int CmpAbs(const MachineInt& n1, const BigInt& N2) noexcept;
  int CmpAbs(const MachineInt& n1, const MachineInt& n2) noexcept;

  bool operator==(const BigInt& N1, const BigInt& N2) noexcept;
  bool operator==(const BigInt& N1, const MachineInt& n2) noexcept;
  bool operator==(const MachineInt& n1, const BigInt& N2) noexcept;

  bool operator!=(const BigInt& N1, const BigInt& N2) noexcept;
  bool operator!=(const BigInt& N1, const MachineInt& n2) noexcept;
  bool operator!=(const MachineInt& n1, const BigInt& N2) noexcept;

  bool operator> (const BigInt& N1, const BigInt& N2) noexcept;
  bool operator> (const BigInt& N1, const MachineInt& n2) noexcept;
  bool operator> (const MachineInt& n1, const BigInt& N2) noexcept;

  bool operator>=(const BigInt& N1, const BigInt& N2) noexcept;
  bool operator>=(const BigInt& N1, const MachineInt& n2) noexcept;
  bool operator>=(const MachineInt& n1, const BigInt& N2) noexcept;

  bool operator< (const BigInt& N1, const BigInt& N2) noexcept;
  bool operator< (const BigInt& N1, const MachineInt& n2) noexcept;
  bool operator< (const MachineInt& n1, const BigInt& N2) noexcept;

  bool operator<=(const BigInt& N1, const BigInt& N2) noexcept;
  bool operator<=(const BigInt& N1, const MachineInt& n2) noexcept;
  bool operator<=(const MachineInt& n1, const BigInt& N2) noexcept;
			

  // MISCELLANEOUS FUNCTIONS

  double mantissa(long& exp, const BigInt& N) noexcept;
//???  double mantissa(const BigInt& N) noexcept;
  long exponent(const BigInt& N) noexcept;
  long MultiplicityOf2(const BigInt& N);
  long MultiplicityOf2(const MachineInt& N);

  double log(const BigInt& N);
  double LogAbs(const BigInt& N);
  // std::log works for machine integers
  long FloorLog2(const MachineInt& n);
  long FloorLog2(const BigInt& n);
  long FloorLog10(const MachineInt& n);
  long FloorLog10(const BigInt& n);

  long FloorLogBase(const MachineInt& n, const MachineInt& base);
  long FloorLogBase(const MachineInt& n, const BigInt& base);
  long FloorLogBase(const BigInt& N,     const MachineInt& base);
  long FloorLogBase(const BigInt& N,     const BigInt& base);

  BigInt RangeFactorial(const MachineInt& lo, const MachineInt& hi); // computes lo*(lo+1)*(lo+2)*...*hi
  BigInt factorial(const BigInt& N);
  BigInt factorial(const MachineInt& n);
  double LogFactorial(const MachineInt& n);
  double LogFactorial(const BigInt& n);
  BigInt primorial(const BigInt& N);
  BigInt primorial(const MachineInt& n);

  BigInt binomial(const BigInt& N,     const BigInt& R);
  BigInt binomial(const MachineInt& n, const BigInt& R);
  BigInt binomial(const BigInt& N,     const MachineInt& r);
  BigInt binomial(const MachineInt& n, const MachineInt& r);
  BigInt fibonacci(const BigInt& N);
  BigInt fibonacci(const MachineInt& n);


  BigInt RoundDiv(const BigInt& num,     const BigInt& den);     ///< round num/den to an integer (halves round away from 0)
  BigInt RoundDiv(const MachineInt& num, const BigInt& den);     ///< round num/den to an integer (halves round away from 0)
  BigInt RoundDiv(const BigInt& num,     const MachineInt& den); ///< round num/den to an integer (halves round away from 0)
  long         RoundDiv(const MachineInt& num, const MachineInt& den); ///< round num/den to an integer (halves round away from 0)

  BigInt FloorRoot(const MachineInt& n, const MachineInt& r);
  BigInt FloorRoot(const MachineInt& n, const BigInt& R);
  BigInt FloorRoot(const BigInt& N,     const MachineInt& r);
  BigInt FloorRoot(const BigInt& N,     const BigInt& R);

  bool IsExactFloorRoot(long& ans,   const MachineInt& n, const MachineInt& r);
  bool IsExactFloorRoot(BigInt& ans, const MachineInt& n, const MachineInt& r);
  bool IsExactFloorRoot(long& ans,   const MachineInt& n, const BigInt& R);
  bool IsExactFloorRoot(BigInt& ans, const MachineInt& n, const BigInt& R);
//???  bool IsExactFloorRoot(long& ans, const BigInt& N, const MachineInt& r);
  bool IsExactFloorRoot(BigInt& ans, const BigInt& N,     const MachineInt& r);
//???  bool IsExactFloorRoot(long& ans, const BigInt& N, const BigInt& R);
  bool IsExactFloorRoot(BigInt& ans, const BigInt& N,     const BigInt& R);

  long FloorSqrt(const MachineInt& n);
  BigInt FloorSqrt(const BigInt& N);

  bool IsSquare(const MachineInt& n);
  bool IsSquare(const BigInt& N) noexcept;
  bool IsPower(const MachineInt& n);
  bool IsPower(const BigInt& N) noexcept;


/*
  friend void           SqrtRem(BigInt& sqrt, BigInt& rem, const BigInt& b);
BigInt         Div2Exp(const BigInt& b, const MachineInt& exponent_of_2);
BigInt         Mod2Exp(const BigInt& b, const MachineInt& exponent_of_2);
BigInt         Mul2Exp(const BigInt& b, const MachineInt& exponent_of_2);
*/


  // Inline impls

  inline int sign(const MachineInt& n) noexcept
  { if (IsNegative(n)) return -1; if (IsZero(n)) return 0; return 1; }

  inline int sign(const BigInt& N) noexcept
  { return mpz_sgn(mpzref(N)); }


  inline void negate(BigInt& lhs) noexcept
  { mpz_neg(mpzref(lhs), mpzref(lhs)); }


  inline bool IsOdd(const BigInt& N) noexcept
  { return mpz_odd_p(mpzref(N)); }

  inline bool IsOdd(const MachineInt& n) noexcept
  { return uabs(n)&1; }


  inline bool IsEven(const BigInt& N) noexcept
  { return mpz_even_p(mpzref(N)); }

  inline bool IsEven(const MachineInt& n) noexcept
  { return !IsOdd(n); }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/BigIntOps.H,v 1.14 2024/01/25 21:06:00 abbott Exp $
// $Log: BigIntOps.H,v $
// Revision 1.14  2024/01/25 21:06:00  abbott
// Summary: Removed thread-safety from SumBigInt/SumBigRat (costly/complicated to do it right)
//
// Revision 1.13  2023/12/10 15:31:18  abbott
// Summary: Cleaned and tidied SumBigInt & SumBigRat
//
// Revision 1.12  2023/12/09 21:37:32  abbott
// Summary: Tidied SumBigInt/SumBigRat; added doc
//
// Revision 1.11  2023/12/01 20:25:28  abbott
// Summary: Added flag to disable overflow check in power (for integers)
//
// Revision 1.10  2022/10/25 17:23:49  abbott
// Summary: Added myAdd; threadsafe version of op+=; need mutex now
//
// Revision 1.9  2022/10/11 18:06:36  abbott
// Summary: Added comment about swap, move ctor & mpve assign
//
// Revision 1.8  2022/10/10 09:52:43  abbott
// Summary: Added SumBigInt class
//
// Revision 1.7  2022/08/08 17:53:45  abbott
// Summary: SizeInBase is obsolete (redmine 1693)
//
// Revision 1.6  2022/02/18 14:11:44  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/08/03 12:27:33  abbott
// Summary: Removed "useless" const qualifier on rtn values (redmine 1606)
//
// Revision 1.4  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.3  2021/01/31 10:08:20  abbott
// Summary: Added noexcept; added LogAbs
//
// Revision 1.2  2019/09/16 14:28:43  abbott
// Summary: Added primorial; changed iroot into FloorRoot, and IsExactIRoot into IsExactFloorRoot
//
// Revision 1.1  2018/05/18 09:15:44  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.11  2017/11/29 20:30:17  abbott
// Summary: Added MultiplicityOf2C
//
// Revision 1.10  2016/08/02 12:49:34  abbott
// Summary: Renamed NumDigits to SizeInBase; updated doc.
//
// Revision 1.9  2015/11/24 12:46:47  abbott
// Summary: Renamed "isqrt" --> "FloorSqrt"
//
// Revision 1.8  2015/11/23 18:20:27  abbott
// Summary: Renamed ILogBase -> FloorLogBase; added FloorLog2, FloorLog10
//
// Revision 1.7  2015/11/05 14:20:06  abbott
// Summary: Changed rtn type of LeastNonNegRemainder to long when modulus is MachineInt
//
// Revision 1.6  2015/10/09 18:28:44  abbott
// Summary: Corrected redmine reference
//
// Revision 1.5  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.4  2014/04/08 13:04:13  abbott
// Summary: Added new fn IsPowerOf2
// Author: JAA
//
// Revision 1.3  2013/05/20 15:48:13  abbott
// Added new fn LogFactorial (placed in IntOperations).
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

#endif
