#ifndef CoCoA_SmallFpDoubleImpl_H
#define CoCoA_SmallFpDoubleImpl_H

//   Copyright (c)  2005,2009,2011-2013  John Abbott, and Anna M. Bigatti

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


// Header file for the class SmallFpDoubleImpl


#include "CoCoA/assert.H"
#include "CoCoA/GlobalManager.H"

#include <cmath>
//using std::floor;
//using std::fmod;


namespace CoCoA
{

  class MachineInt;     // fwd decl -- defined in MachineInt.H
  class BigInt;         // fwd decl -- defined in BigInt.H
  class BigRat;         // fwd decl -- defined in BigRat.H
  class SmallPrime;     // fwd decl -- defined in NumTheory-prime.H

  /*-----------------------------------------------------------------*/
  /** \include SmallFpDoubleImpl.txt  */
  /*-----------------------------------------------------------------*/
  class SmallFpDoubleImpl
  {
  public:
    explicit SmallFpDoubleImpl(const MachineInt& n, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
    explicit SmallFpDoubleImpl(SmallPrime p, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
    static bool IsGoodCtorArg(const MachineInt& n) noexcept;
    static bool IsGoodCtorArg(SmallPrime p) noexcept;
    static long ourMaxModulus() noexcept(!CoCoA_DEBUG_MODE);
  public: // disable copy ctor & assignment
    SmallFpDoubleImpl(const SmallFpDoubleImpl&) = delete;
    SmallFpDoubleImpl& operator=(const SmallFpDoubleImpl&) = delete;

  public:
    typedef double value_t;
    static const int ourDatumSize = sizeof(value_t);

    long myModulus() const noexcept;
    value_t myReduce(const MachineInt& n) const noexcept;  ///< n % myModulus
    value_t myReduce(const BigInt& N) const noexcept;      ///< N % myModulus
    value_t myReduce(const BigRat& q) const;      ///< q % myModulus (error if den(q)%myModulus == 0)
    long myExport(value_t x) const noexcept(!CoCoA_DEBUG_MODE); ///< exports x into a long (according to myResiduesAreSymm)

    value_t myNegate(value_t x) const noexcept(!CoCoA_DEBUG_MODE);
    value_t myAdd(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE);
    value_t mySub(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE);
    value_t myMul(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE);
    value_t myDiv(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE);
    value_t myPower(value_t x, long n) const noexcept;

    bool myIsZeroAddMul(value_t& lhs, value_t y, value_t z) const noexcept(!CoCoA_DEBUG_MODE);///< lhs += y*z, result says whether lhs == 0.
    value_t myNormalize(value_t n) const noexcept(!CoCoA_DEBUG_MODE);
    value_t myHalfNormalize(value_t n) const noexcept(!CoCoA_DEBUG_MODE);
    long myMaxIters() const noexcept { return myIterLimit; }

  private: // Data members
    const value_t myModulusValue;
    const bool myResiduesAreSymm;  // Used only in myExport
    const value_t myResidueUPBValue;
    const long myIterLimit;

  private: // impl details
    static value_t ourCheckCtorArg(const MachineInt& n);
    static value_t ourCheckCtorArg(SmallPrime p);
    static value_t ourCalcResidueUPB(value_t p) noexcept;
    static long ourCalcIterLimit(value_t p) noexcept;
    static constexpr long ourMaxInt() noexcept;
  };  // end of class SmallFpDoubleImpl


  std::ostream& operator<<(std::ostream& out, const SmallFpDoubleImpl& arith);
  bool operator==(const SmallFpDoubleImpl& arith1, const SmallFpDoubleImpl& arith2) noexcept;
  bool operator!=(const SmallFpDoubleImpl& arith1, const SmallFpDoubleImpl& arith2) noexcept;


  //------------------------------------------------------------
  // SmallFpDoubleImpl inline functions
  //------------------------------------------------------------

  inline long SmallFpDoubleImpl::myModulus() const noexcept
  {
    return static_cast<long>(myModulusValue); // cast is safe, cannot overflow.
  }


  inline long SmallFpDoubleImpl::myExport(value_t x) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    if (!myResiduesAreSymm || x <= myModulusValue/2.0)
      return static_cast<long>(x); // cast cannot overflow
    return static_cast<long>(x - myModulusValue); // cast cannot overflow
  }


  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myNormalize(value_t n) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(n >= 0 && n == std::floor(n));
    return std::fmod(n,myModulusValue);
  }


  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myNegate(value_t x) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    if (x == 0) return x;
    return myModulusValue - x;
  }


  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myAdd(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    CoCoA_ASSERT(0 <= y && y < myModulusValue && y == std::floor(y));
    const value_t ans = x+y;
    if (ans < myModulusValue) return ans;
    return ans - myModulusValue;
  }


  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::mySub(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    CoCoA_ASSERT(0 <= y && y < myModulusValue && y == std::floor(y));
    if (x >= y) return x-y;
    return x + (myModulusValue-y);  // avoid trying to create a negative value
  }


  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myMul(value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    CoCoA_ASSERT(0 <= y && y < myModulusValue && y == std::floor(y));
    return std::fmod(x*y, myModulusValue);
  }


  inline bool SmallFpDoubleImpl::myIsZeroAddMul(value_t& lhs, value_t x, value_t y) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= x && x < myModulusValue && x == std::floor(x));
    CoCoA_ASSERT(0 <= y && y < myModulusValue && y == std::floor(y));
    lhs = fmod(lhs+x*y, myModulusValue);
    CoCoA_ASSERT(0 <= lhs && lhs < myModulusValue && lhs == std::floor(lhs));
    return (lhs == 0);
  }


  inline SmallFpDoubleImpl::value_t SmallFpDoubleImpl::myHalfNormalize(value_t n) const noexcept(!CoCoA_DEBUG_MODE)
  {
    CoCoA_ASSERT(0 <= n && n < 2*myResidueUPBValue && n == std::floor(n));
//    if (n < myResidueUPBValue) return n;
//    return n-myResidueUPBValue;
    if (n >= myResidueUPBValue) return n-myResidueUPBValue;
    return n;
  }


  inline bool operator==(const SmallFpDoubleImpl& arith1, const SmallFpDoubleImpl& arith2) noexcept
  {
    return (arith1.myModulus() == arith2.myModulus());
  }

  inline bool operator!=(const SmallFpDoubleImpl& arith1, const SmallFpDoubleImpl& arith2) noexcept
  {
    return !(arith1 == arith2);
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SmallFpDoubleImpl.H,v 1.27 2023/03/27 13:52:45 abbott Exp $
// $Log: SmallFpDoubleImpl.H,v $
// Revision 1.27  2023/03/27 13:52:45  abbott
// Summary: Cleaning: noexcept depends on CoCoA_DEBUG_MODE macro
//
// Revision 1.26  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.25  2021/10/30 17:17:24  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.24  2021/03/04 21:03:45  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.23  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.22  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.21  2020/01/26 14:17:54  abbott
// Summary: Changed include of MachineInt into fwd decl
//
// Revision 1.20  2018/06/25 12:28:20  abbott
// Summary: Ctors now accept SmallPrime args (and skip primality tests)
//
// Revision 1.19  2014/06/17 10:03:15  abbott
// Summary: Changed type of myIterLimit (double -> long); makes more sense.
// Author: JAA
//
// Revision 1.18  2013/05/27 13:10:14  abbott
// Correct param type from long to MachineInt.
//
// Revision 1.17  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.16  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.15  2012/01/30 11:01:11  abbott
// Added naive printing fn (just to help during debugging).
//
// Revision 1.14  2011/11/09 13:50:01  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.13  2011/08/24 10:21:09  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.12  2011/08/17 11:56:12  abbott
// Added two static_casts to keep compiler quiet.
//
// Revision 1.11  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.10  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.9  2011/05/20 09:44:20  abbott
// Removed some no-longer-useful bits of code.
//
// Revision 1.8  2011/05/19 14:38:27  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.7  2011/03/22 20:06:13  abbott
// Added static mem fn IsGoodCtorArg (called by RingFp pseudo-ctors).
// Commented out ctors which take ZZ arg -- seems useless.
//
// Revision 1.6  2011/03/16 15:30:55  abbott
// Changed arg type from MachineInt into const MachineInt&.
//
// Revision 1.5  2011/03/14 10:28:15  abbott
// Changed unsigned long into long (and unsigned short into short).
//
// Revision 1.4  2009/05/14 09:39:29  abbott
// Added possibility to specify "symmetric" or "non-negative" residues
// in quotients of ZZ.  Affects printing of elements in quotients of ZZ
// (also changed printing of elements in general quotient rings).
// Consequent changes in several tests.
//
// Revision 1.3  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.7  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.6  2007/01/11 14:07:42  cocoa
// -- changed names to arguments called "rsh"
//
// Revision 1.5  2006/12/06 17:27:29  cocoa
// -- removed #include "config.H"
//
// Revision 1.4  2006/11/24 17:22:05  cocoa
// -- removed OpenMathFwd.H
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.4  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.3  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
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

#endif
