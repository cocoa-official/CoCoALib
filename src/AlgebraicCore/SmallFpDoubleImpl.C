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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SmallFpDoubleImpl.C,v 1.41 2023/03/27 13:57:33 abbott Exp $
// $Log: SmallFpDoubleImpl.C,v $
// Revision 1.41  2023/03/27 13:57:33  abbott
// Summary: Cleaning: added better defn of ourMaxInt
//
// Revision 1.40  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.39  2021/03/04 21:03:46  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.38  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.37  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.36  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.35  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.34  2019/03/18 11:18:59  abbott
// Summary: Added include for iostream (after changing NumTheory)
//
// Revision 1.33  2018/06/25 12:28:21  abbott
// Summary: Ctors now accept SmallPrime args (and skip primality tests)
//
// Revision 1.32  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.31  2018/02/27 10:54:12  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.30  2018/01/17 10:30:22  abbott
// Summary: Changed InvModNoCheck into InvModNoArgCheck
//
// Revision 1.29  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.28  2016/07/20 08:45:28  abbott
// Summary: Minor change
//
// Revision 1.27  2016/07/04 11:09:38  abbott
// Summary: Minor cleaning (e.g. changed name of internal fn from exgcd to InvMod)
//
// Revision 1.26  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.25  2015/10/09 18:26:36  abbott
// Summary: Corrected redmine reference
//
// Revision 1.24  2015/10/09 18:18:28  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.23  2013/05/27 13:01:24  abbott
// Simplified IsGoodCtorArg (following hint on redmine).
// Some minor cosmetic changes.
//
// Revision 1.22  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.21  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.20  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.19  2012/01/30 23:27:57  abbott
// Added print function.
//
// Revision 1.18  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.17  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/17 11:59:07  abbott
// In myAssign from a QQ, replaced myModulusValue by myModulus()
// [the former is a double, the latter a long]
//
// Revision 1.15  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.14  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.13  2011/05/20 09:44:20  abbott
// Removed some no-longer-useful bits of code.
//
// Revision 1.12  2011/05/19 14:38:27  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.11  2011/03/22 20:06:13  abbott
// Added static mem fn IsGoodCtorArg (called by RingFp pseudo-ctors).
// Commented out ctors which take ZZ arg -- seems useless.
//
// Revision 1.10  2011/03/14 10:28:14  abbott
// Changed unsigned long into long (and unsigned short into short).
//
// Revision 1.9  2010/09/30 14:43:30  abbott
// Added full qualification to calls to std::floor (needed now that GlobalManager.H
// includes FractionField.H).
//
// Revision 1.8  2010/03/23 11:58:04  abbott
// Simplified MaxInt -- removed static data member (troublesome for thread safety).
//
// Revision 1.7  2009/12/23 18:53:52  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.6  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.5  2009/06/05 12:14:55  abbott
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
// Revision 1.4  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.3  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.4  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.3  2006/11/27 13:06:22  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.2  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.3  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.4  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.3  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.2  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.1  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
