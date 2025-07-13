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

  namespace   //anonymous
  {

    long deg(const DenseUPolyClean& f)
    {
      if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
      return f.myDegPlus1()-1;
    }

  }  // end of anonymous namespace
  
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
