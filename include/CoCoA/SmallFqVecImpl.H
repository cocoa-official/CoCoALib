#ifndef CoCoA_SmallFqVecImpl_H
#define CoCoA_SmallFqVecImpl_H

//   Copyright (c)  2015  John Abbott, and Anna M. Bigatti

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

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class FFqImpl_vec
  {
  public: typedef SmallFpImpl::value FpElem;
  public:
    FFqImpl_vec(long p, const std::vector<FpElem>& m);
    FFqImpl_vec(const ideal& I);
//???  FFqImpl_vec(const SmallFpImpl& FFp, const vector<FpElem_t>& m);

//  typedef std::vector<FpElem_t> value_t;
    typedef const FpElem* value_t;
    typedef FpElem* OUTvalue_t;
  private:
    void myMulByGen(std::vector<FpElem>& v) const noexcept;
    long compress(const std::vector<FpElem>& v) const noexcept;
    void MakeLogExpTbls() const;
  public:
    // value myNegate(value_t& x) const;
    // value myRecip(value_t& x) const;
    void myGen(OUTvalue_t ans) const noexcept; // sets ans to the canonical generator
    void myAdd(OUTvalue_t ans, value_t x, value_t y) const noexcept;
    void mySub(OUTvalue_t ans, value_t x, value_t y) const noexcept;
    void myMul(OUTvalue_t ans, value_t x, value_t y) const;
    void myDiv(OUTvalue_t ans, value_t x, value_t y) const;
    void myPower(OUTvalue_t ans, value_t x, long n) const;  // assumes n >= 1
    long myExtnDeg() const noexcept { return myDeg; }
    // value myNormalize(uvalue n) const;
    // uvalue myHalfNormalize(uvalue n) const;
    long myModulus() const noexcept { return myFFp.myModulus(); }
    const SmallFpImpl& myFpArith() const noexcept { return myFFp; }

    std::vector<long> myExport(value_t x) const;               ///< exports x into a vector<long>: k-th entry is coeff of alpha^k (according to myResiduesAreSymm)

  private: // impl details
    void myBinaryPowerLoop(OUTvalue_t ans, value_t x, long n) const;  // assumes n >= 1

  private:
    const SmallFpImpl myFFp; // includes modulus
    const std::vector<FpElem> myMinPoly;
    int myDeg;
    int myTailDeg; // tail degree = deg(MinPoly-LM(MinPoly))
  };



} // end of namespace CoCoA

#endif
