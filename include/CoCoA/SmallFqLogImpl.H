#ifndef CoCoA_SmallFqLogImpl_H
#define CoCoA_SmallFqLogImpl_H

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

#include <vector>
//using std::vector;

namespace CoCoA
{

  class FFqLogImpl
  {
  public:
//???  FFqLogImpl(long p, const vector<SmallFpImpl::value>& m);
    FFqLogImpl(const ideal& I);
//???  FFqImpl_vec(const SmallFpImpl& FFp, const vector<SmallFpImpl::value_t>& m);

    using repr_t = unsigned int;
////  typedef int repr_t; // unsigned?  long?
  private: // impl details
//  void myMulByGen(std::vector<SmallFpImpl::value>& v) const;
//  long compress(const std::vector<SmallFpImpl::value>& v) const;
    void MakeArithTbls(const ideal& I);
    void myMulByGen(const SmallFpImpl& FFp, std::vector<SmallFpImpl::value>& v, std::vector<SmallFpImpl::value>& MonicMinPoly) const;
    long compress(const SmallFpImpl& FFp, const std::vector<SmallFpImpl::value>& v) const noexcept;
  public:
    repr_t myZero() const noexcept { return 0; }
    repr_t myOne() const noexcept { return 1; }
    bool myIsZero(repr_t x) const noexcept { return (x == 0); }
    bool myIsOne(repr_t x) const noexcept { return (x == 1); }
    bool myIsMinusOne(repr_t x) const noexcept { if (myModulus == 2) return myIsOne(x); return (x == 1+myCard/2); }
    repr_t myReduce(long n) const noexcept;
    repr_t myNegate(repr_t x) const noexcept;
    // value myRecip(value_t& x) const;
    repr_t myGen() const noexcept; // sets to the canonical generator
    repr_t myAdd(repr_t x, repr_t y) const noexcept;
    repr_t mySub(repr_t x, repr_t y) const noexcept;
    repr_t myMul(repr_t x, repr_t y) const noexcept;
    repr_t myDiv(repr_t x, repr_t y) const;
    repr_t myPower(repr_t x, long n) const noexcept;  // assumes n >= 1
    long myExtnDeg() const { return myDeg; }
    // value myNormalize(uvalue n) const;
    // uvalue myHalfNormalize(uvalue n) const;


    long myExport(repr_t x) const noexcept;   ///< exports x into a long (using a compressed repr)

  public://????????private: // data members
    const repr_t myModulus;
//  const SmallFpImpl myFFp;
//  const std::vector<SmallFpImpl::value> myMinPoly;
    int myDeg;
    repr_t myCard;
//  int myTailDeg; // tail degree = deg(MinPoly-LM(MinPoly))
    std::vector<repr_t> myAdd1Tbl;
    std::vector<repr_t> myExpTbl;
    std::vector<repr_t> mySmallLogTbl;
  };

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SmallFqLogImpl.H,v 1.5 2022/10/25 17:17:48 abbott Exp $
// $Log: SmallFqLogImpl.H,v $
// Revision 1.5  2022/10/25 17:17:48  abbott
// Summary: Improved indentation
//
// Revision 1.4  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.2  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//

#endif
