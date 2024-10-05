//   Copyright (c)  2007-2017  John Abbott and Anna Bigatti
//   Author:  2007  Anna Bigatti

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


// Source code for class DenseUPolyClean

#include "CoCoA/DenseUPolyClean.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/error.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::max;
using std::copy; // copy ctor
//using std::swap;
#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;


namespace CoCoA
{

  namespace  {    //anonymous
  long deg(const DenseUPolyClean& f)
  {
    if (IsZero(f)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "deg(f)");
    return f.myDegPlus1()-1;
  }
  }
  
  ring CoeffRing(const DenseUPolyClean& f)
  { return f.myCoeffRingValue; }


  DenseUPolyClean::DenseUPolyClean(const ring& R, long MinCapacity):
      myCoeffRingValue(R)
  {
    CoCoA_ASSERT(MinCapacity > 0);
    myCoeffsValue.reserve(MinCapacity);
    myDegPlus1Value = 0;
    mySizeValue = 0;
  }


  DenseUPolyClean::~DenseUPolyClean()
  {
  }


  DenseUPolyClean::DenseUPolyClean(const DenseUPolyClean& copy, long MinCapacity):
      myCoeffRingValue(copy.myCoeffRingValue),
      myDegPlus1Value(copy.myDegPlus1()),
      mySizeValue(copy.mySize())
  {
    CoCoA_ASSERT(MinCapacity > 0);
    myCoeffsValue.reserve(max(MinCapacity, copy.myDegPlus1()));
    myCoeffsValue = copy.myCoeffsValue;
    //    copy(rawcopy.begin(), rawcopy.end(), rawlhs.back_inserter());
  }

  // ANNA: this is commented out to check whether it is really needed
//   DenseUPolyClean::DenseUPolyClean(const DenseUPolyClean& copy):
//       myCoeffRingValue(copy.myCoeffRingValue),
//       myDegPlus1Value(copy.myDegPlus1Value),
//       mySizeValue(copy.mySize())
//   {
//     myCoeffsValue.reserve(copy.myDegPlus1Value);
//     myCoeffsValue = copy.myCoeffsValue;
//     //    copy(rawcopy.begin(), rawcopy.end(), rawlhs.back_inserter());
//   }


  DenseUPolyClean& DenseUPolyClean::operator=(const DenseUPolyClean& rhs)
  {
    if (this == &rhs) return *this;
    if (IsZero(rhs))
      myAssignZero();
    else
    {
      DenseUPolyClean copy(rhs, rhs.myDegPlus1());
      ourSwap(*this, copy);
    }
    return *this;
  }


  DenseUPolyClean& DenseUPolyClean::operator=(const MachineInt& rhs)
  { // exception safe
    if (IsZero(rhs))
      myAssignZero();
    else
    {
      RingElem x(myCoeffRingValue, rhs);
      myResize(1);
      myDegPlus1Value = 1;
      swap(myCoeffsValue[0], x);
    }
    return *this;
  }


  DenseUPolyClean& DenseUPolyClean::operator=(const BigInt& rhs)
  { // exception safe
    //    myAssignZero();  // this would probably save memory
    if (IsZero(rhs))
      myAssignZero();
    else
    {
      RingElem x(myCoeffRingValue, rhs);
      myResize(1);
      myDegPlus1Value = 1;
      swap(myCoeffsValue[0], x);
    }
    return *this;
  }

  DenseUPolyClean& DenseUPolyClean::operator=(const BigRat& rhs)
  { // exception safe
    //    myAssignZero();  // this would probably save memory
    if (IsZero(rhs))
    {
      myAssignZero();
      return *this;
    }
    RingElem x(myCoeffRingValue, rhs);
    myResize(1);
    myDegPlus1Value = 1;
    swap(myCoeffsValue[0], x);
    return *this;
  }


  bool IsCompatible(const DenseUPolyClean& f, const DenseUPolyClean& g)
  {
    return (f.myCoeffRingValue == g.myCoeffRingValue);
  }


  // internal
  inline bool IsZero(const DenseUPolyClean& f)
  {
    return f.myDegPlus1() == 0;
  }

//----------------------------------------------------------------------

  void DenseUPolyClean::ourSwap(DenseUPolyClean& f, DenseUPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    std::swap(f.myCoeffsValue, g.myCoeffsValue);
    std::swap(f.myDegPlus1Value, g.myDegPlus1Value);
    std::swap(f.mySizeValue, g.mySizeValue);
  }


  void DenseUPolyClean::myAssignZero()
  {
    for (long i=0; i < myDegPlus1(); ++i)
      myCoeffsValue[i] = 0;
    myDegPlus1Value = 0;
  }


  void DenseUPolyClean::myResize(long NewSize)
  {
    myCoeffsValue.resize(NewSize, zero(myCoeffRingValue));
    mySizeValue = NewSize;
  }


  void DenseUPolyClean::myResetDeg()
  {
    long i = myDegPlus1();
    for (; i!=0; --i)
      if (!IsZero(myCoeffsValue[i-1])) break;
    myDegPlus1Value = i;  // either i==0 or myCoeffsValue[i-1]!=0
  }


  void DenseUPolyClean::myAssignZeroCoeff(long d)
  { // exception safe?
    CoCoA_ASSERT(d >= 0);
    if (d >= myDegPlus1())  return;
    // d <= deg
    myCoeffsValue[d] = 0;
    if (d == myDegPlus1()-1) myResetDeg();
  }


  void DenseUPolyClean::myAssignNonZeroCoeff(ConstRefRingElem c, long d)
  { // exception safe
    CoCoA_ASSERT(d >= 0);
    CoCoA_ASSERT(d <= mySizeValue-1);
    RingElem x(c);
    swap(myCoeffsValue[d], x);
    if (d >= myDegPlus1()) myDegPlus1Value = d+1;
  }


  long NumTerms(const DenseUPolyClean& f)
  {
    if (IsZero(f)) return 0;
    const long D = deg(f);
    long nterms = 0;
    for (long d=0; d <= D; ++d)
      if (!IsZero(f.myCoeff(d)))
        ++nterms;
    return nterms;
  }


  ConstRefRingElem LC(const DenseUPolyClean& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    CoCoA_ASSERT(!IsZero(f.myCoeffsValue[f.myDegPlus1()]));
    return f.myCoeffsValue[f.myDegPlus1()];
  }


  void DenseUPolyClean::myNegate()
  {
    for (long i=0 ; i < myDegPlus1() ; ++i)
      myCoeffRingValue->myNegate(raw(myCoeffsValue[i]), raw(myCoeffsValue[i]));   // MIGHT THROW???
  }


  void add(DenseUPolyClean& lhs, const DenseUPolyClean& f, const DenseUPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(f, g));
    const ring& R = lhs.myCoeffRingValue;
    DenseUPolyClean ans(R, max(f.myCoeffsValue.capacity(), g.myCoeffsValue.capacity())); // ans.myDegPlus1() is 0
    ans.myResize(max(g.myDegPlus1(), f.myDegPlus1())+1);

    long i = 0;
    if (f.myDegPlus1() >= g.myDegPlus1())
    {
      for (; i<g.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = f.myCoeffsValue[i] + g.myCoeffsValue[i];
      for (; i<f.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = f.myCoeffsValue[i];
    }
    else
    {
      for (; i<f.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = f.myCoeffsValue[i] + g.myCoeffsValue[i];
      for (; i<g.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = g.myCoeffsValue[i];
    }
    ans.mySizeValue = len(ans.myCoeffsValue);
    ans.myDegPlus1Value = max(f.myDegPlus1(), g.myDegPlus1());
    if (f.myDegPlus1() == g.myDegPlus1())
      ans.myResetDeg();
    swap(lhs, ans);
  }


  void sub(DenseUPolyClean& lhs, const DenseUPolyClean& f, const DenseUPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(f, g));
    const ring& R = lhs.myCoeffRingValue;
    DenseUPolyClean ans(R, max(f.myCoeffsValue.capacity(), g.myCoeffsValue.capacity())); // ans.myDegPlus1() is 0
    ans.myResize(max(g.myDegPlus1(), f.myDegPlus1())+1);

    long i = 0;
    if (f.myDegPlus1() >= g.myDegPlus1())
    {
      for (; i<g.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = f.myCoeffsValue[i] - g.myCoeffsValue[i];
      for (; i<f.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = f.myCoeffsValue[i];
    }
    else
    {
      for (; i<f.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = f.myCoeffsValue[i] - g.myCoeffsValue[i];
      for (; i<g.myDegPlus1(); ++i)
        ans.myCoeffsValue[i] = - g.myCoeffsValue[i];
    }
    ans.mySizeValue = len(ans.myCoeffsValue);
    ans.myDegPlus1Value = max(f.myDegPlus1(), g.myDegPlus1());
    if (f.myDegPlus1() == g.myDegPlus1())
      ans.myResetDeg();
    swap(lhs, ans);
  }


  void DenseUPolyClean::myAddMulLM(ConstRefRingElem c, long d, const DenseUPolyClean& g)
  {
    myResize(max(myDegPlus1(), g.myDegPlus1() + d));

    if (IsOne(c)) // special case for Hilbert Poincare Series
      for (long i=0; i<g.myDegPlus1(); ++i)
      {
        if (!IsZero(g.myCoeffsValue[i]))
          myCoeffsValue[i+d] += g.myCoeffsValue[i];
      }
    else
      for (long i=0; i<g.myDegPlus1(); ++i)
      {
        if (!IsZero(g.myCoeffsValue[i]))
          myCoeffsValue[i+d] += c * g.myCoeffsValue[i];
      }
    myDegPlus1Value = max(myDegPlus1(), g.myDegPlus1()+d);
    if (myDegPlus1() == g.myDegPlus1()+d ||
        !IsIntegralDomain(myCoeffRingValue))
      myResetDeg();
  }


  void DenseUPolyClean::myMulByCoeff(ConstRefRingElem c)
  {
    for (long i=0 ; i < myDegPlus1() ; ++i)
      myCoeffsValue[i] *= c;
    if (!IsIntegralDomain(myCoeffRingValue))  myResetDeg();
  }


  void DenseUPolyClean::myDivByCoeff(ConstRefRingElem c)
  {
    for (long i=0 ; i < myDegPlus1() ; ++i)
      myCoeffsValue[i] /= c;
    if (!IsIntegralDomain(myCoeffRingValue))  myResetDeg();
  }


  void DenseUPolyClean::myMulByXExp(long n)
  { // EXCEPTION SAFE
    myResize(myDegPlus1() + n); // EXCEPTION SAFE ???
    for (long i=myDegPlus1(); i != 0; )
    {
      --i;
      swap(myCoeffsValue[i+n], myCoeffsValue[i]);
    }
    myDegPlus1Value = myDegPlus1() + n;
  }


  void DenseUPolyClean::myMulBy1MinusXExp(long n)
  { // NOT EXCEPTION SAFE
    myResize(myDegPlus1() + n);
    for (long i=myDegPlus1(); i != 0; )
    {
      --i;
      myCoeffsValue[i+n] -= myCoeffsValue[i];
    }
    myDegPlus1Value = myDegPlus1() + n;
  }


  void output(ostream& out, const DenseUPolyClean& f)  // for debugging only
  {
    if (!out) return;  // short-cut for bad ostreams
    if (IsZero(f)) { out << "0"; return; }
    for (long i=f.myDegPlus1()-1; i >= 0; --i)
      out << " +(" << f.myCoeff(i) << ")*indet^" << i;
  }


  bool IsMonomial(const DenseUPolyClean& f)
  { // is it useful?
    return NumTerms(f) == 1;
  }


  bool IsEqual(const DenseUPolyClean& f, const DenseUPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    if (&f == &g) return true;
    if (f.myDegPlus1() != g.myDegPlus1()) return false;
    for (long i=0; i < f.myDegPlus1(); ++i)
      if (f.myCoeff(i) != g.myCoeff(i)) return false;
    return true;
  }


  void deriv(DenseUPolyClean& lhs, const DenseUPolyClean& f)
  {
    if (IsZero(f) || deg(f) == 0) { lhs.myAssignZero(); return; }
    const long degf = deg(f);
    DenseUPolyClean ans(CoeffRing(f), (degf-1)+1);
    ans.myResize(degf);
    for (long d=1; d <= degf; ++d)
    {
      const RingElem c = d*f.myCoeff(d);
      if (IsZero(c)) continue;
      ans.myAssignNonZeroCoeff(c, d-1);
    }
    swap(lhs, ans);
  }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DenseUPolyClean.C,v 1.35 2024/07/02 09:57:57 bigatti Exp $
// $Log: DenseUPolyClean.C,v $
// Revision 1.35  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.34  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.33  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.32  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.31  2020/02/11 16:12:17  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.30  2019/11/14 17:55:43  abbott
// Summary: Removec commented out use of clog
//
// Revision 1.29  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.28  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.27  2017/03/28 12:56:29  bigatti
// -- changed copyright
//
// Revision 1.26  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.25  2014/04/30 16:04:32  abbott
// Summary: Replaced X.size() by len(X)
// Author: JAA
//
// Revision 1.24  2013/05/28 17:08:43  abbott
// A fair amount of cleaning (needs more, but it's not urgent).
//
// Revision 1.23  2013/05/22 13:35:37  bigatti
// -- fix for deriv
// -- StdDeg --> deg  fixed and moved inside anonymous namespace
//
// Revision 1.22  2013/05/22 10:07:48  bigatti
// -- modified implementation of "=" for x = zero(P)
//
// Revision 1.21  2013/05/20 15:56:10  abbott
// Added new fns CoeffRing and StdDeg.
// Added impls for deriv, NumTerms, IsMonomial.
//
// Revision 1.20  2013/04/17 16:30:28  abbott
// Replaced size_t by long.
//
// Revision 1.19  2013/03/15 17:50:39  abbott
// Changed int parameter into long.
//
// Revision 1.18  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.17  2011/11/09 14:03:40  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.16  2011/08/24 10:24:17  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.15  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.14  2011/04/27 08:21:48  bigatti
// -- added gcd with coefficients in GCDDomain
//
// Revision 1.13  2009/10/30 10:15:41  abbott
// Cleaned up the include directives: removed unnecessary ones,
// put them into alphabetical order.
//
// Revision 1.12  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.11  2008/05/30 12:48:32  abbott
// Commented out some unused formal parameters (to avoid compiler warnings).
//
// Revision 1.10  2008/01/29 16:16:48  bigatti
// -- use of member functions instead of member fields for rvalues
//
// Revision 1.9  2007/12/21 12:29:08  bigatti
// -- abstract implementation in DenseUPolyRing of myDiv, myIsDivisible, myIsInvertible, myGcd
// -- abstract implementation in DenseUPolyRing of ideal operations
// -- some cleaning
//
// Revision 1.8  2007/11/29 17:43:58  bigatti
// -- fixed  operator=  which used automatic copy-ctor
// -- copy-ctor declared but not defined for checking
//
// Revision 1.7  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/10/19 10:04:23  bigatti
// -- RingDenseUPolyClean now allow to specify the MinCapacity for all
//    coeff vectors (to avoid too many reallocations)
//
// Revision 1.5  2007/10/19 09:03:46  bigatti
// -- the two contructors reserve 150 slots for coefficients.
//
// Revision 1.4  2007/10/15 13:43:24  bigatti
// -- improved myAddMul for c==1 (called from RecPoincare)
//
// Revision 1.3  2007/10/15 12:40:31  bigatti
// -- new implementation for myAddMul
// -- new: coeffs between deg+1 and size are now guaranteed to be 0
//
// Revision 1.2  2007/10/10 14:02:37  bigatti
// -- added myMulBy1MinusXExp
// -- fixed a few little bugs
//
// Revision 1.1  2007/10/05 15:28:56  bigatti
// -- added abstract class DenseUPolyRing for representing dense
//    univariate polynomials
// -- added concrete class RingDenseUPolyClean, the cleanest
//    implementation
// -- just for testing, still horribly inefficient and incomplete
//
