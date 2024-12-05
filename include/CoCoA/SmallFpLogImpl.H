#ifndef CoCoA_SmallFpLogImpl_H
#define CoCoA_SmallFpLogImpl_H

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


// Header file for the class SmallFpLogImpl


#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/GlobalManager.H"


#include <vector>
// using std::vector


namespace CoCoA
{

  class MachineInt;     // fwd decl -- defined in MachineInt.H
  class BigInt;         // fwd decl -- defined in BigInt.H
  class BigRat;         // fwd decl -- defined in BigRat.H
  class SmallPrime;     // fwd decl -- defined in NumTheory-prime.H

  class SmallFpLogImpl
  {
  public:
    explicit SmallFpLogImpl(const MachineInt& n, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
    explicit SmallFpLogImpl(SmallPrime p, GlobalSettings::ResidueRepr repr = DefaultResidueRepr());
    static bool IsGoodCtorArg(const MachineInt& n) noexcept;
    static bool IsGoodCtorArg(SmallPrime p) noexcept;
    static long ourMaxModulus() noexcept;
  public: // disable copy ctor & assignment
    SmallFpLogImpl(const SmallFpLogImpl&) = delete;
    SmallFpLogImpl& operator=(const SmallFpLogImpl&) = delete;

  public:
    typedef unsigned int value_t; // data type for internal repr of mod p values
    static const int ourDatumSize = sizeof(value_t);

    long myModulus() const noexcept;
    value_t myReduce(const MachineInt& n) const noexcept;  ///< n % myModulus
    value_t myReduce(const BigInt& N) const noexcept;      ///< N % myModulus
    value_t myReduce(const BigRat& q) const;               ///< q % myModulus (error if den(q)%myModulus == 0)
    long myExport(value_t x) const noexcept;               ///< exports x into a long (according to myResiduesAreSymm)

    value_t myNegate(value_t x) const noexcept;
    value_t myAdd(value_t x, value_t y) const noexcept;
    value_t mySub(value_t x, value_t y) const noexcept;
    value_t myMul(value_t x, value_t y) const noexcept;
    value_t myDiv(value_t x, value_t y) const;
    value_t myPower(value_t x, long n) const noexcept;

    bool myIsZeroAddMul(value_t& lhs, value_t y, value_t z) const noexcept;///< lhs += y*z, result says whether lhs == 0.
    value_t myNormalize(value_t n) const noexcept;
    value_t myHalfNormalize(value_t n) const noexcept;
    long myMaxIters() const noexcept { return myIterLimit; }

  private: // data members
    const value_t myModulusValue;
    const bool myResiduesAreSymm;   // Used only in myExport
    const value_t myResidueUPBValue;
    const value_t myIterLimit;
    const value_t myRoot;
    typedef unsigned short FpTableElem; // the log/exp table elements are of this type
    std::vector<FpTableElem> myLog;
    std::vector<FpTableElem> myExp;

  private: // impl details
    static value_t ourCheckCtorArg(const MachineInt& n);
    static value_t ourCheckCtorArg(SmallPrime p);
    static value_t ourCalcResidueUPB(value_t p) noexcept;
    static long ourCalcIterLimit(value_t p) noexcept;
    void myCtorBody() noexcept;
  };  // end of class SmallFpLogImpl


  std::ostream& operator<<(std::ostream& out, const SmallFpLogImpl& arith);
  bool operator==(const SmallFpLogImpl& arith1, const SmallFpLogImpl& arith2) noexcept;
  bool operator!=(const SmallFpLogImpl& arith1, const SmallFpLogImpl& arith2) noexcept;


  //------------------------------------------------------------
  // SmallFpLogImpl inline functions
  //------------------------------------------------------------

  inline long SmallFpLogImpl::myModulus() const noexcept
  {
    return myModulusValue; // implicit cast is safe
  }


  inline long SmallFpLogImpl::myExport(value_t x) const noexcept
  {
    if (!myResiduesAreSymm || x <= myModulusValue/2) return x; // implicit cast is safe
    return -static_cast<long>(myModulusValue - x); // cast is safe
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::myNormalize(value_t n) const noexcept
  {
    return n%myModulusValue;
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::myNegate(value_t x) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    if (x == 0) return x;
    return myModulusValue - x;
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::myAdd(value_t x, value_t y) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    const value_t ans = x+y;
    if (ans < myModulusValue) return ans;
    return ans - myModulusValue;
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::mySub(value_t x, value_t y) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    if (x >= y) return x-y;
    return x + (myModulusValue-y);  // avoid trying to create a negative value
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::myMul(value_t x, value_t y) const noexcept
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    if (x == 0 || y == 0)  return 0;
    else return myExp[myLog[x] + myLog[y]];
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::myDiv(value_t x, value_t y) const
  {
    CoCoA_ASSERT(x == myNormalize(x));
    CoCoA_ASSERT(y == myNormalize(y));
    if (y == 0)  CoCoA_THROW_ERROR1(ERR::DivByZero);
    if (x == 0)  return 0;
    return myExp[myModulusValue-1 + myLog[x] - myLog[y]];
  }


  inline bool SmallFpLogImpl::myIsZeroAddMul(value_t& lhs, value_t y, value_t z) const noexcept
  {
    CoCoA_ASSERT(y == myNormalize(y));
    CoCoA_ASSERT(z == myNormalize(z));
    CoCoA_ASSERT(lhs == myNormalize(lhs));
    lhs += myMul(y,z);
    if (lhs >= myModulusValue) lhs -= myModulusValue;
    return (lhs == 0);
  }


  inline SmallFpLogImpl::value_t SmallFpLogImpl::myHalfNormalize(value_t n) const noexcept
  {
//    if (n < myResidueUPBValue) return n;
//    return n-myResidueUPBValue;
    if (n >= myResidueUPBValue) return n-myResidueUPBValue;
    return n;
  }


  inline bool operator==(const SmallFpLogImpl& arith1, const SmallFpLogImpl& arith2) noexcept
  {
    return (arith1.myModulus() == arith2.myModulus());
  }

  inline bool operator!=(const SmallFpLogImpl& arith1, const SmallFpLogImpl& arith2) noexcept
  {
    return !(arith1 == arith2);
  }


} // end of namespace CoCoA

#endif
