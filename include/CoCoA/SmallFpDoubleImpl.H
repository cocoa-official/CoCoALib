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

#endif
