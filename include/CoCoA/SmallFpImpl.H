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

#endif
