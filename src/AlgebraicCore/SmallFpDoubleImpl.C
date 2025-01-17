//   Copyright (c)  2005,2009,2011-2013  John Abbott and Anna M. Bigatti

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

#include "CoCoA/SmallFpDoubleImpl.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/error.H"

#include <cmath>
using std::floor;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;

namespace CoCoA
{

  // MaxInt gives the largest integer N such that all integers from 1 to N
  // N can be represented exactly using the type SmallFpDoubleImpl::value_t.
  constexpr long SmallFpDoubleImpl::ourMaxInt() noexcept
  {
    constexpr double MAX = 2 / numeric_limits<value_t>::epsilon();
    return static_cast<long>(MAX);
  }




  namespace // anonymous namespace
  {
    // This "double" version is slightly faster than "InvModNoArgCheck<unsigned long>" (in NumTheory.C)
    // on some platforms, and for some inputs; on other platforms it can be slower
    // (JAA observed about 20% slower than InvModNoArgCheck on one platform).
    
    // JAA believes that this routine always returns a reduced value,
    // but does not have a proof of this currently.
    SmallFpDoubleImpl::value_t InvModNoArgCheck(SmallFpDoubleImpl::value_t r, SmallFpDoubleImpl::value_t m)
    {
      CoCoA_ASSERT(r > 0);  CoCoA_ASSERT(m >= 2);
      CoCoA_ASSERT(r < m);
      SmallFpDoubleImpl::value_t p = m;
      SmallFpDoubleImpl::value_t cr = 1;
      SmallFpDoubleImpl::value_t cm = 0; // this is minus what you think it is!
      while (r != 0)
      {
        SmallFpDoubleImpl::value_t q = std::floor(m/r);
	m -= q*r;
	cm += q*cr;
	if (m == 0) break;
	q = std::floor(r/m);
	r -= q*m;
	cr += q*cm;
      }
      if (r+m != 1) CoCoA_THROW_ERROR(ERR::DivByZero, "SmallFpDoubleImpl::myDiv");
      if (r == 0) return p-cm;
      return cr;
    }


  } // end of anonymous namespace


  // These two inline fns are used only in the ctors.
  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::ourCalcResidueUPB(value_t p) noexcept
  {
    return p*(ourMaxInt()/p/2); // Largest (integer) multiple of p <= MaxInt/2.
  }

  inline long SmallFpDoubleImpl::ourCalcIterLimit(value_t p) noexcept
  {
return (ourMaxInt()/(p-1)/(p-1)/2); // Max no. of unreduced products you can sum without exceeding MaxInt/2.
    // const long MaxIters = (ourMaxInt()/(p-1)/(p-1)/2); // Max no. of unreduced products you can sum without exceeding MaxInt/2.
    // constexpr long MaxLong = numeric_limits<long>::max();
    // if (MaxIters > MaxLong) return MaxLong; // JAA reckons this'll never happen.
    // return static_cast<long>(MaxIters);
  }


  SmallFpDoubleImpl::SmallFpDoubleImpl(const MachineInt& n, GlobalSettings::ResidueRepr repr):
      myModulusValue(ourCheckCtorArg(n)),
      myResiduesAreSymm(repr == GlobalSettings::ResidueRepr::symmetric),
      myResidueUPBValue(ourCalcResidueUPB(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue))
  {}

  SmallFpDoubleImpl::SmallFpDoubleImpl(SmallPrime p, GlobalSettings::ResidueRepr repr):
      myModulusValue(ourCheckCtorArg(p)),
      myResiduesAreSymm(repr == GlobalSettings::ResidueRepr::symmetric),
      myResidueUPBValue(ourCalcResidueUPB(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue))
  {}


  bool SmallFpDoubleImpl::IsGoodCtorArg(const MachineInt& n) noexcept
  {
    if (IsNegative(n) || !IsSignedLong(n)) return false;
    const long N = AsSignedLong(n);
    return N <= ourMaxModulus() && IsPrime(N);
  }

  bool SmallFpDoubleImpl::IsGoodCtorArg(SmallPrime p) noexcept
  {
    return p <= ourMaxModulus();
  }

  long SmallFpDoubleImpl::ourMaxModulus() noexcept(!CoCoA_DEBUG_MODE)
  {
    const double SqrtMax2 = std::floor(sqrt(ourMaxInt()/2))-1;
    const long ans = PrevPrime(static_cast<long>(SqrtMax2));
    CoCoA_ASSERT(ourCalcIterLimit(ans) > 0);
    return ans;
  }


  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myReduce(const MachineInt& n) const noexcept
  {
    const value_t ans =  uabs(n)%myModulus();
    if (!IsNegative(n) || ans == 0) return ans;
    return myModulusValue - ans;
  }

  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myReduce(const BigInt& N) const noexcept
  {
    return mpz_fdiv_ui(mpzref(N), myModulus());
  }

  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myReduce(const BigRat& q) const
  {
    const value_t D = mpz_fdiv_ui(mpq_denref(mpqref(q)), myModulus());
    if (D == 0) CoCoA_THROW_ERROR(ERR::DivByZero, "SmallFpDoubleImpl::myReduce");
    const value_t N = mpz_fdiv_ui(mpq_numref(mpqref(q)), myModulus());
    return myDiv(N, D);
  }


  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myDiv(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    CoCoA_ASSERT(0 <= y && y < myModulusValue && y == std::floor(y));
    CoCoA_ASSERT(y != 0);
    // static_cast to ulong is OK because myModulusValue is small enough (see ctor arg check)
    const value_t recip = InvModNoArgCheck(static_cast<unsigned long>(y),
                                           static_cast<unsigned long>(myModulusValue));
    // const value_t recip = InvModNoArgCheck(y, myModulusValue);
    return myMul(x, recip);
  }


  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myPower(value_t x, long n) const noexcept
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    CoCoA_ASSERT(x != 0 || n > 0);  // complain about any non-positive power of 0
    if (x == 0) { return 0; }
    if (x == 1) { return 1; }
    const unsigned long Mminus1 = myModulus() - 1;
    n %= Mminus1; // OK by fermat's little theorem.
    if (n < 0) n += Mminus1;
    if (n == 0) { return 1; }
    if (n == 1) { return x; }
    // Here we know that n >= 2 and x is not 0 or 1.
    // Below is an iterative version of repeated squaring.
    unsigned long mask = 1;
    unsigned long quartern = n/4;
    while (mask <= quartern) mask <<= 1;
    value_t ans = x;
    for (; mask != 0; mask >>= 1)
    {
      ans = myMul(ans, ans);
      if (n & mask) ans = myMul(ans, x);
    }
    return ans;
  }


  // These two fns check that the modulus is a suitably small prime.
  // If the arg is not suitable then an exception is thrown, otherwise
  // the return value is the suitably small prime.
  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::ourCheckCtorArg(const MachineInt& n)
  {
    if (!IsGoodCtorArg(n))
      CoCoA_THROW_ERROR(ERR::BadSmallFpChar, "SmallFpDoubleImpl ctor");
    return AsUnsignedLong(n);
  }

  SmallFpDoubleImpl::value_t SmallFpDoubleImpl::ourCheckCtorArg(SmallPrime p)
  {
    if (!IsGoodCtorArg(p))
      CoCoA_THROW_ERROR(ERR::BadSmallFpChar, "SmallFpDoubleImpl ctor");
    return p;
  }


  std::ostream& operator<<(std::ostream& out, const SmallFpDoubleImpl& arith)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "SmallFpDoubleImpl(" << arith.myModulus() << ")";
  }


} // end of namespace CoCoA
