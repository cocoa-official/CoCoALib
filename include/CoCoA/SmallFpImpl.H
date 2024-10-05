#ifndef CoCoA_SmallFpImpl_H
#define CoCoA_SmallFpImpl_H

//   Copyright (c)  2005,2009,2011-2015,2018  John Abbott, and Anna M. Bigatti

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


// Header file for the class SmallFpImpl


#include "CoCoA/assert.H"
#include "CoCoA/GlobalManager.H"


namespace CoCoA
{

  class MachineInt;     // fwd decl -- defined in MachineInt.H
  class BigInt;         // fwd decl -- defined in BigInt.H
  class BigRat;         // fwd decl -- defined in BigRat.H
  class SmallPrime;     // fwd decl -- defined in NumTheory-prime.H

  // Enum to give arg type for fns zero and one
  enum class SmallFpEnum { marker };
  constexpr SmallFpEnum SmallFp = SmallFpEnum::marker; // safely inject into namespace CoCoA

  class SmallFpImpl
  {
  public: // Type used to repr the values in the finite field (semi-private, but C++ needs it to be public)
    using repr_t = unsigned long;  // ***MUST*** be an unsigned integral type
  public:
    class value;
    class NonRedValue // an unreduced/unnormalized value (as produced by unnorm add/mul).
    {
    public:
      NonRedValue() noexcept: myVal(0) {}
      NonRedValue(value x) noexcept: myVal(x.myVal) {}
      // default copy ctor, assignment, dtor all OK
      friend class SmallFpImpl; // for myNormalize and myHalfNormalize (see below)
      friend NonRedValue operator+(NonRedValue, NonRedValue) noexcept;
      friend NonRedValue operator*(NonRedValue, NonRedValue) noexcept;
      friend NonRedValue& operator+=(NonRedValue&, NonRedValue) noexcept;
      friend NonRedValue& operator*=(NonRedValue&, NonRedValue) noexcept;
      friend std::ostream& operator<<(std::ostream& out, NonRedValue x);
//      friend repr_t GetRepr(NonRedValue x);
    private:
      /*implicit*/ NonRedValue(repr_t x) noexcept: myVal(x) {}
    private: // data member -- must be same as data member of value!!
      repr_t myVal;
    };

    class value // a reduced/normalized value
    {
    public:
      value() noexcept: myVal(0) {}
      // default copy ctor, assignment, dtor all OK
      friend class SmallFpImpl;
      friend class NonRedValue;
      friend value zero(SmallFpEnum) noexcept;
      friend value one(SmallFpEnum) noexcept;
      friend bool IsZero(value x) noexcept;
      friend bool IsOne(value x) noexcept;
      friend bool operator==(value x, value y) noexcept;
      friend std::ostream& operator<<(std::ostream& out, value x);
    private:
      /*implicit*/ value(repr_t n) noexcept: myVal(n) {}
    private: // data member -- must be same as data member of NonRedValue!!
      repr_t myVal;
    };

  public:
    explicit SmallFpImpl(const MachineInt& n, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
    explicit SmallFpImpl(SmallPrime p, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
    static bool IsGoodCtorArg(const MachineInt& n) noexcept;
    static bool IsGoodCtorArg(SmallPrime p) noexcept;
    static long ourMaxModulus() noexcept;
  public: // disable copy ctor & assignment
    SmallFpImpl(const SmallFpImpl&) = delete;
    SmallFpImpl& operator=(const SmallFpImpl&) = delete;

  public:
    static const int ourDatumSize = sizeof(value);

    long myModulus() const noexcept;
    long myMaxIters() const noexcept { return myIterLimit; }
    value myReduce(const MachineInt& n) const noexcept;  ///< n % myModulus
    value myReduce(const BigInt& N) const noexcept;      ///< N % myModulus
    value myReduce(const BigRat& q) const;      ///< q % myModulus (error if den(q)%myModulus == 0)

    long myExportNonNeg(value x) const noexcept;         ///< exports x into a long (least non neg residue)
    long myExportSymm(value x) const noexcept;           ///< exports x into a long (symm residue)
    long myExport(value x) const noexcept;               ///< exports x into a long (according to setting myResiduesAreSymm)

    value myNegate(value x) const noexcept;
    value myRecip(value x) const;
    value myAdd(value x, value y) const noexcept;
    value mySub(value x, value y) const noexcept;
    value myMul(value x, value y) const noexcept;
    value myDiv(value x, value y) const;
    value myPower(value x, long n) const noexcept;
    value myAddMul(value x, value y, value z) const noexcept; /// x+y*z
    bool myIsZeroAddMul(value& lhs, value y, value z) const noexcept;///< lhs += y*z, result says whether lhs == 0.
    value myNormalize(NonRedValue n) const noexcept;
    NonRedValue myHalfNormalize(NonRedValue n) const noexcept; // 2 HalfNormalized values may be added without causing overflow

  private: // Data members
    const repr_t myModulusValue;
    const repr_t myHalfWayPoint;
    const long myIterLimit;
    const bool myResiduesAreSymm;   // Used only in myExport

    friend std::ostream& operator<<(std::ostream& out, const SmallFpImpl& arith); // for access to myResiduesAreSymm

  private: // impl details
    static repr_t ourCheckCtorArg(const MachineInt& n);
    static repr_t ourCheckCtorArg(SmallPrime p);
    static repr_t ourCalcHalfwayPoint(repr_t p) noexcept;
    static long ourCalcIterLimit(repr_t p) noexcept;
  };  // end of class SmallFpImpl


  std::ostream& operator<<(std::ostream& out, const SmallFpImpl& arith);
  bool operator==(const SmallFpImpl& arith1, const SmallFpImpl& arith2) noexcept;
  bool operator!=(const SmallFpImpl& arith1, const SmallFpImpl& arith2) noexcept;

  std::ostream& operator<<(std::ostream& out, SmallFpImpl::NonRedValue x);
  std::ostream& operator<<(std::ostream& out, SmallFpImpl::value x);

  inline SmallFpImpl::value zero(SmallFpEnum) noexcept { return 0; }
  inline SmallFpImpl::value one(SmallFpEnum) noexcept { return 1; }

  //------------------------------------------------------------
  // SmallFpImpl inline functions
  //------------------------------------------------------------

  inline SmallFpImpl::NonRedValue operator+(SmallFpImpl::NonRedValue x, SmallFpImpl::NonRedValue y) noexcept
  {
    CoCoA_ASSERT(x.myVal <= std::numeric_limits<decltype(x.myVal)>::max()-y.myVal);
    return x.myVal + y.myVal;  // NOTE: implicit conversion using private ctor
  }

  inline SmallFpImpl::NonRedValue operator*(SmallFpImpl::NonRedValue x, SmallFpImpl::NonRedValue y) noexcept
  {
    CoCoA_ASSERT(y.myVal == 0 || x.myVal <= std::numeric_limits<decltype(x.myVal)>::max()/y.myVal);
    return x.myVal * y.myVal;  // NOTE: implicit conversion using private ctor
  }

  inline SmallFpImpl::NonRedValue& operator+=(SmallFpImpl::NonRedValue& lhs, SmallFpImpl::NonRedValue x) noexcept
  {
    CoCoA_ASSERT(lhs.myVal <= std::numeric_limits<decltype(x.myVal)>::max()-x.myVal);
    lhs.myVal += x.myVal;
    return lhs;
  }

  inline SmallFpImpl::NonRedValue& operator*=(SmallFpImpl::NonRedValue& lhs, SmallFpImpl::NonRedValue x) noexcept
  {
    CoCoA_ASSERT(x.myVal == 0 || lhs.myVal <= std::numeric_limits<decltype(x.myVal)>::max()/x.myVal);
    lhs.myVal *= x.myVal;
    return lhs;
  }


  inline SmallFpImpl::value SmallFpImpl::myAddMul(SmallFpImpl::value x, SmallFpImpl::value y, SmallFpImpl::value z) const noexcept
  {
    CoCoA_ASSERT(y.myVal == 0 || z.myVal <= std::numeric_limits<repr_t>::max()/y.myVal);
    CoCoA_ASSERT(y.myVal*z.myVal <= std::numeric_limits<repr_t>::max()-x.myVal);
    return myNormalize(x.myVal + y.myVal * z.myVal);
  }


  inline bool IsZero(SmallFpImpl::value x) noexcept
  {
    return x.myVal == 0;
  }

  inline bool IsOne(SmallFpImpl::value x) noexcept
  {
    return x.myVal == 1;
  }


  inline bool operator==(SmallFpImpl::value x, SmallFpImpl::value y) noexcept
  {
    return x.myVal == y.myVal;
  }

  inline bool operator!=(SmallFpImpl::value x, SmallFpImpl::value y) noexcept
  {
    return !(x == y);
  }




  inline long SmallFpImpl::myModulus() const noexcept
  {
    return myModulusValue; // implicit cast is safe
  }


  inline long SmallFpImpl::myExportNonNeg(value x) const noexcept
  {
    return x.myVal; // implicit cast is safe
  }

  inline long SmallFpImpl::myExportSymm(value x) const noexcept
  {
    if (x.myVal <= myModulusValue/2) return x.myVal; // implicit cast is safe
    return -static_cast<long>(myModulusValue - x.myVal); // cast is safe
  }

  inline long SmallFpImpl::myExport(value x) const noexcept
  {
    if (myResiduesAreSymm) return myExportSymm(x);
    return myExportNonNeg(x);
  }


  inline SmallFpImpl::value SmallFpImpl::myNormalize(NonRedValue n) const noexcept
  {
    return n.myVal%myModulusValue;
  }

  inline SmallFpImpl::NonRedValue SmallFpImpl::myHalfNormalize(NonRedValue n) const noexcept
  {
    if (n.myVal >= myHalfWayPoint)
      n.myVal -= myHalfWayPoint;
    return n;
  }


  inline SmallFpImpl::value SmallFpImpl::myNegate(value x) const noexcept /*assert*/
  {
    CoCoA_ASSERT(x == myNormalize(x));
    if (IsZero(x)) return x;
    return myModulusValue - x.myVal;
  }


  inline SmallFpImpl::value SmallFpImpl::myAdd(value x, value y) const noexcept /*assert*/
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    const repr_t ans = x.myVal + y.myVal;
    if (ans < myModulusValue) return ans;
    return ans - myModulusValue;
  }


  inline SmallFpImpl::value SmallFpImpl::mySub(value x, value y) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    if (x.myVal >= y.myVal) return x.myVal - y.myVal;
    return x.myVal + (myModulusValue-y.myVal);  // avoid creating a negative value
  }


  inline SmallFpImpl::value SmallFpImpl::myMul(value x, value y) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    return (x.myVal*y.myVal)%myModulusValue; // overflow cannot occur!
  }


  inline bool SmallFpImpl::myIsZeroAddMul(value& lhs, value y, value z) const noexcept
  {
    CoCoA_ASSERT(y == myNormalize(y));
    CoCoA_ASSERT(z == myNormalize(z));
    CoCoA_ASSERT(lhs == myNormalize(lhs));
    lhs = (lhs.myVal + y.myVal * z.myVal)%myModulusValue;  // overflow cannot occur!
    return IsZero(lhs);
  }


  inline bool operator==(const SmallFpImpl& arith1, const SmallFpImpl& arith2) noexcept
  {
    return (arith1.myModulus() == arith2.myModulus());
  }

  inline bool operator!=(const SmallFpImpl& arith1, const SmallFpImpl& arith2) noexcept
  {
    return !(arith1 == arith2);
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SmallFpImpl.H,v 1.41 2023/03/27 13:51:33 abbott Exp $
// $Log: SmallFpImpl.H,v $
// Revision 1.41  2023/03/27 13:51:33  abbott
// Summary: Cleaning: updated some assertions to C++11
//
// Revision 1.40  2022/11/16 19:41:21  abbott
// Summary: Replaced typedef by using
//
// Revision 1.39  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.38  2021/10/30 17:17:39  abbott
// Summary: Used keyword delete (redmine 1627)
//
// Revision 1.37  2021/03/04 21:03:45  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.36  2021/03/03 22:13:25  abbott
// Summary: enum class (redmine 894)
//
// Revision 1.35  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.34  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.33  2020/01/26 14:17:55  abbott
// Summary: Changed include of MachineInt into fwd decl
//
// Revision 1.32  2018/06/25 12:28:20  abbott
// Summary: Ctors now accept SmallPrime args (and skip primality tests)
//
// Revision 1.31  2018/02/27 10:52:42  abbott
// Summary: Updated to offer special ctor for SmallPrime arg (no need to check primality)
//
// Revision 1.30  2018/02/15 17:26:11  abbott
// Summary: Added EratosthenesRange, and PrimeSeq
//
// Revision 1.29  2018/01/18 14:50:05  abbott
// Summary: Corrected a check in operator* or NonRedValue
//
// Revision 1.28  2018/01/17 10:27:25  abbott
// Summary: Corrected private implicit ctor for value; added useful comments
//
// Revision 1.27  2018/01/11 17:16:19  abbott
// Summary: Removed embarassing bug in private ctor for NonRedValue (arg type was uint)
//
// Revision 1.26  2016/03/25 20:22:48  abbott
// Summary: Minor reordering of source code
//
// Revision 1.25  2015/12/11 13:08:11  abbott
// Summary: Added myExportNonNeg and myExportSymm
//
// Revision 1.24  2015/12/02 15:30:49  abbott
// Summary: Added op!=
//
// Revision 1.23  2015/11/04 10:09:47  abbott
// Summary: MAJOR REVISION to handle cleanly non-reduced values
//
// Revision 1.22  2015/10/09 18:28:44  abbott
// Summary: Corrected redmine reference
//
// Revision 1.21  2015/10/09 18:18:27  abbott
// Summary: Renamed "abs" to "uabs" for MachineInt; new fn "negate"; see redmine 783
//
// Revision 1.20  2014/11/27 11:28:51  abbott
// Summary: Corrected signature of mul
// Author: JAA
//
// Revision 1.19  2014/11/18 15:46:48  abbott
// Summary: intermediate check in (new and old impls together)
// Author: JAA
//
// Revision 1.18  2013/05/27 14:05:44  abbott
// Removed typedefs SmallFpElem_t and SmallFpLogElem_t from config.H;
// they are now nest typedefs SmallFpImpl::value_t and SmallFpLogImpl::value_t.
//
// Revision 1.17  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.16  2013/02/21 13:00:37  abbott
// Added new fn operator!=
//
// Revision 1.15  2012/11/23 17:26:45  abbott
// Added mem fn myResidueUPB.
//
// Revision 1.14  2012/10/04 08:57:28  abbott
// Added keyword const to a local variable.
//
// Revision 1.13  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.12  2012/01/30 11:01:11  abbott
// Added naive printing fn (just to help during debugging).
//
// Revision 1.11  2011/11/09 13:50:01  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.10  2011/08/24 10:21:09  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.9  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.8  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.7  2011/05/19 14:38:27  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.6  2011/03/22 20:06:13  abbott
// Added static mem fn IsGoodCtorArg (called by RingFp pseudo-ctors).
// Commented out ctors which take ZZ arg -- seems useless.
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
// Revision 1.5  2006/12/06 17:30:02  cocoa
// -- rearranged #include
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
// Revision 1.4  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/02/28 15:58:56  cocoa
// Resynch after some minor changes.
//
// Revision 1.2  2005/02/11 16:45:24  cocoa
// Removed the useless and misleading functions myInit and myKill
// from the SmallFp*Impl classes; various consequential changes.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.7  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.6  2004/11/11 13:05:49  cocoa
// -- added \include *.txt for doxygen
//
// Revision 1.5  2004/11/04 18:47:42  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.4  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.3  2004/07/16 15:45:12  cocoa
// First stage of new RingElem implementation completed.
//
// Revision 1.2  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
// Revision 1.1  2004/07/13 16:32:26  cocoa
// First stage of major revamp of ring elements.
// Implementation of RingFp has been split into "ring interface"
// and "algorithms plus data structures".
//
//


#endif
