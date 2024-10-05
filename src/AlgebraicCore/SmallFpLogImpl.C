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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SmallFpLogImpl.C,v 1.35 2022/02/18 14:11:58 abbott Exp $
// $Log: SmallFpLogImpl.C,v $
// Revision 1.35  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.34  2022/02/05 20:01:33  abbott
// Summary: Minor changes to avoid warnings from clang (redmine 1528)
//
// Revision 1.33  2021/03/04 21:03:46  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.32  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.31  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.30  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.29  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.28  2019/03/18 11:19:35  abbott
// Summary: Added include for iostream (after changing NumTheory)
//
// Revision 1.27  2018/06/25 12:28:21  abbott
// Summary: Ctors now accept SmallPrime args (and skip primality tests)
//
// Revision 1.26  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.25  2018/02/27 10:55:14  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.24  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.23  2015/10/09 18:26:37  abbott
// Summary: Corrected redmine reference
//
// Revision 1.22  2015/10/09 18:18:28  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.21  2014/05/02 13:56:26  abbott
// Summary: Added explanatopry comment about a CoCoA_ASSERT
// Author: JAA
//
// Revision 1.20  2014/04/30 16:13:32  abbott
// Summary: Changed arg of sqrt to be double
// Author: JAA
//
// Revision 1.19  2013/05/27 13:01:23  abbott
// Simplified IsGoodCtorArg (following hint on redmine).
// Some minor cosmetic changes.
//
// Revision 1.18  2013/04/29 09:05:23  abbott
// Changed local variable to unsigned to avoid compiler warning.
//
// Revision 1.17  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.16  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.15  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.14  2012/01/30 23:27:57  abbott
// Added print function.
//
// Revision 1.13  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.12  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.11  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.10  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.9  2011/05/19 14:38:26  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.8  2011/03/23 21:00:46  abbott
// Removed FindPrimRoot from NumTheory.H because it was already
// available as PrimitiveRoot (a better name).
// Updated documentation for NumTheory.
//
// Revision 1.7  2011/03/22 20:06:13  abbott
// Added static mem fn IsGoodCtorArg (called by RingFp pseudo-ctors).
// Commented out ctors which take ZZ arg -- seems useless.
//
// Revision 1.6  2011/03/14 10:28:14  abbott
// Changed unsigned long into long (and unsigned short into short).
//
// Revision 1.5  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
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
// Revision 1.5  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.4  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.3  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.2  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.3  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.2  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.2  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.1  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
//
//
