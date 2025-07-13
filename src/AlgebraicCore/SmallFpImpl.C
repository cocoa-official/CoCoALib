//   Copyright (c)  2004-2007,2009,2011,2013  John Abbott and Anna M. Bigatti

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

#include "CoCoA/SmallFpImpl.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"

#include <cmath>
using std::floor;
using std::sqrt;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;

namespace CoCoA
{

  // These two inline fns are used only in the ctors.
  inline SmallFpImpl::repr_t SmallFpImpl::ourCalcHalfwayPoint(repr_t p) noexcept
  {
    CoCoA_ASSERT(p >= 2);
    constexpr repr_t MAX = numeric_limits<repr_t>::max();
    return p*(MAX/p/2); // Largest multiple of p <= MAX/2; exploits integer division.
  }

  inline long SmallFpImpl::ourCalcIterLimit(repr_t p) noexcept
  {
    CoCoA_ASSERT(p >= 2);
    constexpr repr_t MAX = numeric_limits<repr_t>::max();
    const repr_t MaxIters = (MAX/(p-1))/(p-1)/2; // Max no. of unreduced products you can sum without exceeding MAX/2.
    constexpr unsigned long MaxLong = numeric_limits<long>::max();
    if (MaxIters > MaxLong) return MaxLong; // JAA reckons this'll never happen.
    return MaxIters; // implicit cast is safe
  }


  SmallFpImpl::SmallFpImpl(const MachineInt& n, GlobalSettings::ResidueRepr repr):
      myModulusValue(ourCheckCtorArg(n)),
      myHalfWayPoint(ourCalcHalfwayPoint(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue)),
      myResiduesAreSymm(repr == GlobalSettings::ResidueRepr::symmetric)
  {}

  SmallFpImpl::SmallFpImpl(SmallPrime p, GlobalSettings::ResidueRepr repr):
      myModulusValue(ourCheckCtorArg(p)),
      myHalfWayPoint(ourCalcHalfwayPoint(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue)),
      myResiduesAreSymm(repr == GlobalSettings::ResidueRepr::symmetric)
  {}

  bool SmallFpImpl::IsGoodCtorArg(const MachineInt& n) noexcept
  {
    if (IsNegative(n) || !IsSignedLong(n)) return false;
    const long N = AsSignedLong(n);
    return N <= ourMaxModulus() && IsPrime(N);
  }

  bool SmallFpImpl::IsGoodCtorArg(SmallPrime p) noexcept
  {
    return p <= ourMaxModulus();
  }

  long SmallFpImpl::ourMaxModulus() noexcept
  {
    constexpr double HalfMaxIntVal = static_cast<double>(numeric_limits<repr_t>::max()/2); // may harmlessly round up
    const long ans = PrevPrime(ConvertTo<long>(std::floor(sqrt(HalfMaxIntVal))));
    CoCoA_ASSERT(ourCalcIterLimit(ans) > 0); // check that 2*ans^2 <= MAXLONG
    return ans;
  }


  SmallFpImpl::value SmallFpImpl::myReduce(const MachineInt& n) const noexcept
  {
    const repr_t ans =  uabs(n)%myModulusValue;
    if (!IsNegative(n) || ans == 0) return ans;
    return myModulusValue - ans;
  }

  SmallFpImpl::value SmallFpImpl::myReduce(const BigInt& N) const noexcept
  {
    return mpz_fdiv_ui(mpzref(N), myModulusValue);
  }

  SmallFpImpl::value SmallFpImpl::myReduce(const BigRat& q) const
  {
    if (mpz_cmp_ui(mpq_denref(mpqref(q)),1) == 0) // short-cut if den(q) == 1
      return mpz_fdiv_ui(mpq_numref(mpqref(q)), myModulusValue);
    const repr_t D = mpz_fdiv_ui(mpq_denref(mpqref(q)), myModulusValue);
    if (D == 0)  CoCoA_THROW_ERROR1(ERR::DivByZero);
    const repr_t N = mpz_fdiv_ui(mpq_numref(mpqref(q)), myModulusValue);
    return myDiv(N, D);
  }


  SmallFpImpl::value SmallFpImpl::myRecip(value x) const
  {
    CoCoA_ASSERT(x == myNormalize(x));
    if (IsZero(x))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    return InvModNoArgCheck(x.myVal, myModulusValue);
  }


  SmallFpImpl::value SmallFpImpl::myDiv(value x, value y) const
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    if (IsZero(y))  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (IsZero(x)) return 0;
    const value yrecip = InvModNoArgCheck(y.myVal, myModulusValue);
    return myMul(x, yrecip);
  }


  SmallFpImpl::value SmallFpImpl::myPower(value x, long n) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(!IsZero(x) || n > 0); // complain about any non-positive power of 0
    if (IsZero(x)) { return 0; }
    if (x.myVal == 1) { return 1; }
    n %= (myModulusValue-1); // OK by fermat's little theorem.
    if (n < 0) n += myModulusValue-1;
    if (n == 0) { return 1; }
    if (n == 1) { return x; }
    // Here we know that n >= 2 and x is not 0 or 1.
    // Below is an iterative version of repeated squaring.
    unsigned long mask = 1;
    unsigned long quartern = n/4;
    while (mask <= quartern) mask <<= 1;
    value ans = x;
    for (; mask != 0; mask >>= 1)
    {
      ans = myMul(ans, ans);
      if (n & mask)  ans = myMul(ans, x);
    }
    return ans;
  }


  // If p is a small prime, return p as a repr_t (unsigned integral type).
  // Otherwise throw an exception.
  SmallFpImpl::repr_t SmallFpImpl::ourCheckCtorArg(const MachineInt& n)
  {
    if (!IsGoodCtorArg(n))
      CoCoA_THROW_ERROR2(ERR::BadSmallFpChar, "SmallFpImpl ctor");
    return AsUnsignedLong(n);
  }

  SmallFpImpl::repr_t SmallFpImpl::ourCheckCtorArg(SmallPrime p)
  {
    if (!IsGoodCtorArg(p))
      CoCoA_THROW_ERROR2(ERR::BadSmallFpChar, "SmallFpImpl ctor");
    return p;
  }


  std::ostream& operator<<(std::ostream& out, SmallFpImpl::NonRedValue x)
  {
    if (!out)  return out;  // short-cut for bad ostreams
    return out << "!!" << x.myVal << "!!"; // "!!" to emphasise that internal repr is non-reduced
  }

  std::ostream& operator<<(std::ostream& out, SmallFpImpl::value x)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << '!' << x.myVal << '!';  // '!' to emphasise that internal repr is printed
  }

  std::ostream& operator<<(std::ostream& out, const SmallFpImpl& arith)
  {
    if (!out) return out;  // short-cut for bad ostreams
    out << "SmallFpImpl(" << arith.myModulus() << ", export=";
    if (arith.myResiduesAreSymm) out << "SymmResidues";
    else out << "NonNegResidues";
    out << ")";
    return out;
  }


} // end of namespace CoCoA
