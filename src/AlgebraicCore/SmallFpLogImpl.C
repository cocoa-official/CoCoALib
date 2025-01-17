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

#include "CoCoA/SmallFpLogImpl.H"

#include "CoCoA/BigRat.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/error.H"

#include <algorithm>
using std::min;
#include <cmath>
using std::floor;
using std::sqrt;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits; // only in SmallFpLogImpl ctor

namespace CoCoA
{

  // These two inline fns are used only in the ctors.
  inline SmallFpLogImpl::value_t SmallFpLogImpl::ourCalcResidueUPB(value_t p) noexcept
  {
    const value_t MAX = numeric_limits<value_t>::max();
    return p*(MAX/p/2); // Largest multiple of p <= MAX/2; exploits integer division.
  }

  inline long SmallFpLogImpl::ourCalcIterLimit(value_t p) noexcept
  {
    const value_t MAX = numeric_limits<value_t>::max();
    const value_t MaxIters = MAX/(p-1)/(p-1)/2; // Max no. of unreduced products you can sum without exceeding MAX/2.
    constexpr unsigned long UPB = 30000; /////MaxLong = numeric_limits<long>::max();
    if (MaxIters > UPB) return UPB; // JAA reckons this'll never happen.
    return MaxIters;
  }


  SmallFpLogImpl::SmallFpLogImpl(const MachineInt& p, GlobalSettings::ResidueRepr repr):
      myModulusValue(ourCheckCtorArg(p)),
      myResiduesAreSymm(repr == GlobalSettings::ResidueRepr::symmetric),
      myResidueUPBValue(ourCalcResidueUPB(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue)),
      myRoot(PrimitiveRoot(myModulusValue)),
      myLog(myModulusValue),
      myExp(2*myModulusValue-1)
  {
    myCtorBody();
  }

  SmallFpLogImpl::SmallFpLogImpl(SmallPrime p, GlobalSettings::ResidueRepr repr):
      myModulusValue(ourCheckCtorArg(p)),
      myResiduesAreSymm(repr == GlobalSettings::ResidueRepr::symmetric),
      myResidueUPBValue(ourCalcResidueUPB(myModulusValue)),
      myIterLimit(ourCalcIterLimit(myModulusValue)),
      myRoot(PrimitiveRoot(myModulusValue)),
      myLog(myModulusValue),
      myExp(2*myModulusValue-1)
  {
    myCtorBody();
  }


  bool SmallFpLogImpl::IsGoodCtorArg(const MachineInt& n) noexcept
  {
    if (IsNegative(n) || !IsSignedLong(n)) return false;
    const long N = AsSignedLong(n);
    return N <= ourMaxModulus() && IsPrime(N);
  }

  bool SmallFpLogImpl::IsGoodCtorArg(SmallPrime p) noexcept
  {
    return p <= ourMaxModulus();
  }


  long SmallFpLogImpl::ourMaxModulus() noexcept
  {
    const double HalfMaxIntVal = numeric_limits<value_t>::max()/2;
    const long candidate1 = PrevPrime(min<long>(65535, numeric_limits<FpTableElem>::max()));
    const long candidate2 = PrevPrime(static_cast<long>(std::floor(sqrt(HalfMaxIntVal))));
    const long ans = min(candidate1, candidate2);
    CoCoA_ASSERT(ourCalcIterLimit(ans) > 0); // check that 2*ans^2 < MAXLONG
    return ans;
  }



  SmallFpLogImpl::value_t SmallFpLogImpl::myReduce(const MachineInt& n) const noexcept
  {
    const value_t ans =  uabs(n)%myModulusValue;
    if (!IsNegative(n) || ans == 0) return ans;
    return myModulusValue - ans;
  }

  SmallFpLogImpl::value_t SmallFpLogImpl::myReduce(const BigInt& N) const noexcept
  {
    return mpz_fdiv_ui(mpzref(N), myModulusValue);
  }


  SmallFpLogImpl::value_t SmallFpLogImpl::myReduce(const BigRat& q) const
  {
    const value_t D = mpz_fdiv_ui(mpq_denref(mpqref(q)), myModulusValue);
    if (IsZero(D))  CoCoA_THROW_ERROR(ERR::DivByZero, "SmallFpLogImpl::myReduce");
    const value_t N = mpz_fdiv_ui(mpq_numref(mpqref(q)), myModulusValue);
    return myDiv(N, D);
  }


  SmallFpLogImpl::value_t SmallFpLogImpl::myPower(value_t x, long n) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(x != 0 || n > 0); // complain about any non-positive power of 0
    if (x == 0) { return 0; }
    if (x == 1) { return 1; }
    n %= myModulusValue-1;
    if (n < 0) n += myModulusValue-1;
    if (n == 0) { return 1; }
    if (n == 1) { return x; }
    // At this point 0 <= n < myModulusValue-1 and x != 0 and x != 1
    // (n*myLog[x]) cannot overflow, since square of myModulusValue-1 fits in a value_t (see ctor)
    return myExp[(n*myLog[x])%(myModulusValue-1)];  // no risk of overflow in the product
  }



  // If p is a small prime, return p as a value_t.
  // Otherwise throw an exception.
  SmallFpLogImpl::value_t SmallFpLogImpl::ourCheckCtorArg(const MachineInt& n)
  {
    if (!IsGoodCtorArg(n))
      CoCoA_THROW_ERROR(ERR::BadSmallFpChar, "SmallFpLogImpl ctor");
    return AsUnsignedLong(n);
  }

  SmallFpLogImpl::value_t SmallFpLogImpl::ourCheckCtorArg(SmallPrime p)
  {
    if (!IsGoodCtorArg(p))
      CoCoA_THROW_ERROR(ERR::BadSmallFpChar, "SmallFpLogImpl ctor");
    return p;
  }


  // NOTE: When this is called when the Log/Exp tables are still empty.
  void SmallFpLogImpl::myCtorBody() noexcept
  {
    const value_t p = myModulusValue; // Convenient short-hand.

    // Build the log/exp tables.
    if (p == 2) // simplest to handle this special case separately
    {
      myLog[1] = 0;
      myExp[0] = myExp[1] = 1;
      return;
    }

    // From here on we know p is odd.
    const value_t p1 = p-1;
    const value_t half_p1 = (p-1)/2;
    const value_t ThreeHalves = p1 + half_p1;

    value_t s = 1;
    for (value_t i = 0; i < half_p1; ++i)
    {
      // All assignments to myLog and myExp truncate safely.
      myLog[s] = i;
      myLog[p-s] = i+half_p1;
      myExp[i] = s;
      myExp[i+half_p1] = p-s;
      myExp[i+p1] = s;
      myExp[i+ThreeHalves] = p-s;
      s = (s*myRoot)%p;
    }
  }


  std::ostream& operator<<(std::ostream& out, const SmallFpLogImpl& arith)
  {
    if (!out) return out;  // short-cut for bad ostreams
    return out << "SmallFpLogImpl(" << arith.myModulus() << ")";
  }


} // end of namespace CoCoA
