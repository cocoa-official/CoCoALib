//   Copyright (c)  2001-2017  John Abbott and Anna M. Bigatti
//   Author:  2001-2009  John Abbott

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

#include "CoCoA/PPMonoid.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoidEvOv.H" //  for NewPPMonoid
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H" // for OrdMat
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::sort;
#include <iostream>
using std::ostream;
#include <numeric>
using std::accumulate;
#include <set>
using std::set;
// #include <vector>
using std::vector;

namespace CoCoA
{

  PPMonoidBase::PPMonoidBase(const PPOrdering& ord, const std::vector<symbol>& IndetNames):
    IntrusiveReferenceCount(),
    myOrd(ord),
    myIndetSymbols(IndetNames),
    myNumIndets(len(myIndetSymbols))
  {
    CoCoA_ASSERT(NumIndets(ord) == myNumIndets);
    CoCoA_ASSERT(AreDistinct(myIndetSymbols));
    CoCoA_ASSERT(AreArityConsistent(myIndetSymbols));
  }


  namespace // for functions local to this file/compilation unit.
  {
    inline void CheckCompatible(const ConstRefPPMonoidElem& pp1, const ConstRefPPMonoidElem& pp2, const ErrorContext& ErrCtx)
    {
      if (owner(pp1) != owner(pp2))
        CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::MixedPPMs, ErrCtx);
    }
  }


  inline PPMonoidElemRawPtr PPMonoidBase::myNewCheckVecSize(const std::vector<long>& v) const
  {
    if (len(v) != myNumIndets)  CoCoA_THROW_ERROR1(ERR::BadArraySize);
    return myNew(v);
  }

  inline PPMonoidElemRawPtr PPMonoidBase::myNewCheckVecSize(const std::vector<BigInt>& v) const
  {
    if (len(v) != myNumIndets)  CoCoA_THROW_ERROR1(ERR::BadArraySize);
    return myNew(v);
  }


  // Default definition of a virtual function.
  PPMonoidElemRawPtr PPMonoidBase::myNew(const std::vector<BigInt>& EXPV) const
  {
    // Attempt to convert to vector<long> then call myNew on the result
    CoCoA_ASSERT(len(EXPV) == myNumIndets);
    vector<long> expv(myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
      expv[i] = ConvertTo<long>(EXPV[i]);
    return myNew(expv);
  }


  // Default definition of a virtual function.
  bool PPMonoidBase::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myCmp(rawpp1, rawpp2) == 0;
  }


  std::ostream& operator<<(std::ostream& out, const PPMonoid& PPM)
  {
    if (!out) return out;  // short-cut for bad ostreams
    PPM->myOutputSelf(out);
    return out;
  }


  // Default impl, uses small exp version or throws ERR::ArgTooBig
  void PPMonoidBase::myMulIndetPower(RawPtr rawpp, long var, const BigInt& EXP) const
  {
    CoCoA_ASSERT(EXP >= 0);
    long exp;
    if (!IsConvertible(exp, EXP))  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    myMulIndetPower(rawpp, var, exp);
  }


  // Assume inputs are mathematically valid (i.e. exp >= 0).
  // Deal with all trivial cases; pass other cases to myPowerSmallExp.
  void PPMonoidBase::myPower(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const
  {
    CoCoA_ASSERT(exp >= 0);
    if (exp==0 || myIsOne(rawpp1))
    {
      myAssignOne(rawpp);
      return;
    }
    myPowerSmallExp(rawpp, rawpp1, exp);
  }

  // Assume inputs are mathematically valid (i.e. exp >= 0).
  // Deal with all trivial cases; pass other cases to myPowerSmallExp or myPowerBigExp.
  void PPMonoidBase::myPower(RawPtr rawpp, ConstRawPtr rawpp1, const BigInt& EXP) const
  {
    CoCoA_ASSERT(EXP >= 0);
    if (IsZero(EXP) || myIsOne(rawpp1))
    {
      myAssignOne(rawpp);
      return;
    }
    long exp;
    if (IsConvertible(exp, EXP))
      myPowerSmallExp(rawpp, rawpp1, exp);
    else
      myPowerBigExp(rawpp, rawpp1, EXP);
  }


  // Default defn of virtual fn: just gives an error.
  void PPMonoidBase::myPowerBigExp(RawPtr /*rawpp*/, ConstRawPtr /*rawpp1*/, const BigInt& /*EXP*/) const
  { CoCoA_THROW_ERROR1(ERR::ExpTooBig); }


  // Generic implementation -- default defn of virtual function.
  bool PPMonoidBase::myIsIndetPosPower(long& indet, BigInt& EXP, ConstRawPtr rawpp) const
  {
    if (myIsOne(rawpp)) return false; // changed in 0.99700 (previously was error)
    vector<long> expv(myNumIndets);  // SLUG wasteful new+delete ???
    myExponents(expv, rawpp);
    // ??? it should be myBigExponents 2010-02-02
    long TmpIndet = 0;
    for (long i = 1; i < myNumIndets; ++i)
    {
      if (expv[i] == 0) continue;
      if (expv[TmpIndet] != 0) return false;
      TmpIndet = i;
    }
    indet = TmpIndet;
    EXP = expv[indet];
    return true;
  }


  bool PPMonoidBase::myIsIndetPosPower(long& indet, long& pow, ConstRawPtr rawpp) const
  {
    if (myIsOne(rawpp)) return false; // changed in 0.99700 (previous gave error)
    vector<long> expv(myNumIndets);  // SLUG wasteful new+delete ???
    myExponents(expv, rawpp);
    long TmpIndet = 0;
    for (long i = 1; i < myNumIndets; ++i)
    {
      if (expv[i] == 0) continue;
      if (expv[TmpIndet] != 0) return false;
      TmpIndet = i;
    }
    indet = TmpIndet;
    pow = expv[indet];
    return true;
  }


  bool PPMonoidBase::myIsIndetPosPower(ConstRawPtr rawpp) const
  {
    if (myIsOne(rawpp)) return false; // changedin 0.99700 (previously gave error)
    vector<long> expv(myNumIndets);  // SLUG wasteful new+delete ???
    myExponents(expv, rawpp);
    long j = myNumIndets;
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (expv[i] == 0) continue;
      if (j != myNumIndets) return false;
      j = i;
    }
    return j != myNumIndets;
  }


  // Generic PP printing routine.
  void PPMonoidBase::myOutput(std::ostream& out, ConstRawPtr rawpp) const
  {
    if (!out) return;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));

    bool all0 = true;
///    vector<BigInt> expv(myNumIndets);
///    myBigExponents(expv, rawpp);
    BigInt d;
    for (long indet=0; indet < myNumIndets; ++indet)
    {
      myBigExponent(d, rawpp, indet); // Genericity is more important than efficiency here.
///      const BigInt& d = expv[indet];
      if (IsZero(d)) continue;
      if (!all0) out << "*";
      all0 = false;
      out << myIndetSymbol(indet);
      if (d > 1) out << "^" << d;
      else if (d < 0) out << "^(" << d << ")"; // ...in case we ever allow negative exponents.
    }
    if (all0) out << "1";
  }


  void PPMonoidBase::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawpp) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "PPExponents");
    BigInt d;
    for (long indet=0; indet < myNumIndets; ++indet)
    {
      myBigExponent(d, rawpp, indet); // Genericity is more important than efficiency here.
      OMOut << d;
    }
    OMOut->mySendApplyEnd();
  }


  const symbol& IndetSymbol(const PPMonoid& PPM, long idx)
  { 
    if (idx < 0 || idx >= NumIndets(PPM))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    return PPM->myIndetSymbol(idx);
  }


  const std::vector<symbol>& symbols(const PPMonoid& PPM)
  {
    return PPM->mySymbols();
  }


  ConstRefPPMonoidElem PPMonoidBase::mySymbolValue(const symbol& s) const
  {
    for (long i=0; i < myNumIndets; ++i)
      if ( s == myIndetSymbol(i) )
        return myIndets()[i];
    CoCoA_THROW_ERROR2(ERR::BadArg, "unknown symbol in PPM");
    return myOne(); // Just to keep the compiler quiet
  }


  //----------------------------------------------------------------------
  // For fairly obvious reasons assignment (i.e. operator=) does not permit
  // polymorphism on the destination type (e.g. to avoid slicing); so I must
  // define these operators in addition to those for PPMonoidElem.
  // All four operators should have identical definitions: bear this in mind
  // should you ever want to change one of them!

  RefPPMonoidElem& RefPPMonoidElem::operator=(const RefPPMonoidElem& rhs)
  {
    CheckCompatible(*this, rhs, CoCoA_ERROR_CONTEXT);
    // if (this == &rhs) return *this;  // trivial check not really needed
    myPPM->myAssign(myPPPtr, raw(rhs));
    return *this;
  }


  RefPPMonoidElem& RefPPMonoidElem::operator=(const ConstRefPPMonoidElem& rhs)
  {
    CheckCompatible(*this, rhs, CoCoA_ERROR_CONTEXT);
    // if (this == &rhs) return *this;  // trivial check not really needed
    myPPM->myAssign(myPPPtr, raw(rhs));
    return *this;
  }


  ConstMatrixView OrdMat(const PPMonoid& PPM)
  { return OrdMat(ordering(PPM)); }

  ConstMatrixView GradingMat(const PPMonoid& PPM)
  { return GradingMat(ordering(PPM)); }


  //----------------------------------------------------------------------
  // Implementation of class PPMonoidElem below...

  //-------------------- constructors & destructor --------------------//


  PPMonoidElem::PPMonoidElem(const PPMonoid& PPM):
      RefPPMonoidElem(PPM, PPM->myNew())
  {}


  PPMonoidElem::PPMonoidElem(const PPMonoid& PPM, const std::vector<long>& v):
      RefPPMonoidElem(PPM, PPM->myNewCheckVecSize(v))
  {}

  PPMonoidElem::PPMonoidElem(const PPMonoid& PPM, const std::vector<BigInt>& v):
      RefPPMonoidElem(PPM, PPM->myNewCheckVecSize(v))
  {}



  // This function assumes ownership of the value pointed to by ToBeOwned
  PPMonoidElem::PPMonoidElem(const PPMonoid& PPM, PPMonoidElemRawPtr rawToBeOwned):
      RefPPMonoidElem(PPM, rawToBeOwned)
  {}


  // This is a copy ctor.
  PPMonoidElem::PPMonoidElem(const PPMonoidElem& copy):
      RefPPMonoidElem(owner(copy), owner(copy)->myNew(raw(copy)))
  {}


  // This is a copy ctor.
  PPMonoidElem::PPMonoidElem(const ConstRefPPMonoidElem& copy):
      RefPPMonoidElem(owner(copy), owner(copy)->myNew(raw(copy)))
  {}


  // These assignment operators are replicated for RefPPMonoidElem because
  // C++ cannot safely allow polymorphism in the lhs (e.g. slicing troubles).
  // All four operators should have identical definitions: bear this in mind
  // should you ever want to change one of them!

  PPMonoidElem& PPMonoidElem::operator=(const PPMonoidElem& rhs)
  {
    CheckCompatible(*this, rhs, CoCoA_ERROR_CONTEXT);
    // if (this == &rhs) return *this;  // trivial check not really needed
    myPPM->myAssign(myPPPtr, raw(rhs));
    return *this;
  }

  PPMonoidElem& PPMonoidElem::operator=(const ConstRefPPMonoidElem& rhs)
  {
    CheckCompatible(*this, rhs, CoCoA_ERROR_CONTEXT);
    // if (this == &rhs) return *this;  // trivial check not really needed
    myPPM->myAssign(myPPPtr, raw(rhs));
    return *this;
  }


  PPMonoidElem::~PPMonoidElem()
  {
    myPPM->myDelete(myPPPtr);
  }



  long StdDeg(ConstRefPPMonoidElem pp)
  {
    return owner(pp)->myStdDeg(raw(pp));
  }

  long deg(ConstRefPPMonoidElem pp)
  {
    return owner(pp)->myStdDeg(raw(pp));
  }


  degree wdeg(ConstRefPPMonoidElem pp)
  {
    degree ans(GradingDim(ordering(owner(pp))));
    owner(pp)->myWDeg(ans, raw(pp));
    return ans;
  }


  int CmpWDeg(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myCmpWDeg(raw(pp1), raw(pp2));
  }


  int CmpWDegPartial(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2, long PartialGradingDim)
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    if ( PartialGradingDim < 0 || PartialGradingDim > GradingDim(owner(pp1)) )
      CoCoA_THROW_ERROR2(ERR::BadArg, "PartialGradingDim > GradingDim");
    return owner(pp1)->myCmpWDegPartial(raw(pp1), raw(pp2), PartialGradingDim);
  }


  long exponent(ConstRefPPMonoidElem pp, long indet) // degree in indet
  {
    if (indet >= NumIndets(owner(pp)))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    return owner(pp)->myExponent(raw(pp), indet);
  }


  BigInt BigExponent(ConstRefPPMonoidElem pp, long indet)
  {
    if (indet >= NumIndets(owner(pp)))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    BigInt ans;
    owner(pp)->myBigExponent(ans, raw(pp), indet);
    return ans;
  }


  const std::vector<long>& exponents(std::vector<long>& expv, ConstRefPPMonoidElem pp)
  {
    expv.resize(NumIndets(owner(pp)));
    owner(pp)->myExponents(expv, raw(pp));
    return expv;
  }


  std::vector<long> exponents(ConstRefPPMonoidElem pp)
  {
    std::vector<long> expv;
    return exponents(expv, pp);
  }
  

  const std::vector<BigInt>& BigExponents(std::vector<BigInt>& expv, ConstRefPPMonoidElem pp)
  {
    expv.resize(NumIndets(owner(pp)));
    owner(pp)->myBigExponents(expv, raw(pp));
    return expv;
  }


  std::vector<BigInt> BigExponents(ConstRefPPMonoidElem pp)
  {
    std::vector<BigInt> expv;
    BigExponents(expv, pp);
    return expv;
  }


  void swap(RefPPMonoidElem pp1, RefPPMonoidElem pp2)
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    owner(pp1)->mySwap(raw(pp1), raw(pp2));
  }


  bool IsOne(ConstRefPPMonoidElem pp)
  {
    return owner(pp)->myIsOne(raw(pp));
  }


  bool IsIndet(long& index, ConstRefPPMonoidElem pp)
  {
    return owner(pp)->myIsIndet(index, raw(pp));
  }


  bool IsIndet(ConstRefPPMonoidElem pp)
  {
    long NoUse;
    return IsIndet(NoUse, pp);
  }


  bool IsIndetPosPower(long& index, BigInt& EXP, ConstRefPPMonoidElem pp)
  {
    if (IsOne(pp)) return false;
    return owner(pp)->myIsIndetPosPower(index, EXP, raw(pp));
  }


  bool IsIndetPosPower(long& index, long& pow, ConstRefPPMonoidElem pp)
  {
    if (IsOne(pp)) return false;
    return owner(pp)->myIsIndetPosPower(index, pow, raw(pp));
  }


  bool IsIndetPosPower(ConstRefPPMonoidElem pp)
  {
    if (IsOne(pp)) return false;
    return owner(pp)->myIsIndetPosPower(raw(pp));
  }


  // ------------------------------ Comparisons ------------------------------

  int cmp(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)                  // <0, =0, >0 as pp1 < = > pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myCmp(raw(pp1), raw(pp2));
  }


  bool operator==(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 == pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myIsEqual(raw(pp1), raw(pp2));
  }


  bool operator!=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 != pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return !owner(pp1)->myIsEqual(raw(pp1), raw(pp2));
  }


  bool operator<(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)           // pp1 < pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) < 0;
  }


  bool operator<=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 <= pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) <= 0;
  }


  bool operator>(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)           // pp1 > pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) > 0;
  }


  bool operator>=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 >= pp2
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) >= 0;
  }


// ------------------------------ Arithmetic ------------------------------

  PPMonoidElem operator*(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)   //  pp1*pp2;
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myMul(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem operator/(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)   // pp1/pp2;
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    if (!owner(pp1)->myIsDivisible(raw(pp1), raw(pp2)))  CoCoA_THROW_ERROR1(ERR::BadQuot);
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myDiv(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem colon(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)       //  pp1:pp2;
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myColon(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem gcd(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)         // gcd(pp1,pp2);
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myGcd(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem lcm(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)         // lcm(pp1,pp2);
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myLcm(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem radical(ConstRefPPMonoidElem pp)         // radical(pp)
  {
    PPMonoidElem ans(owner(pp));
    owner(pp)->myRadical(raw(ans), raw(pp));
    return ans;
  }


  std::vector<long> IndetsIn(ConstRefPPMonoidElem pp)
  {
    // Simple & universal rather than efficient.
    const long n = NumIndets(owner(pp));
    vector<BigInt> v; v.reserve(n);
    BigExponents(v, pp);
    vector<long> ans;
    for (long i=0; i < n; ++i)
      if (v[i] != 0)
        ans.push_back(i);
    return ans;
  }


  PPMonoidElem power(ConstRefPPMonoidElem pp, long exp)                        // pp^exp
  {
    if (exp < 0)  CoCoA_THROW_ERROR1(ERR::NegExp);
    PPMonoidElem ans(owner(pp));
    owner(pp)->myPower(raw(ans), raw(pp), exp);
    return ans;
  }


  PPMonoidElem power(ConstRefPPMonoidElem pp, const BigInt& EXP)                 // pp^EXP
  {
    if (EXP < 0)  CoCoA_THROW_ERROR1(ERR::NegExp);
    PPMonoidElem ans(owner(pp));
    owner(pp)->myPower(raw(ans), raw(pp), EXP);
    return ans;
  }


  void PowerOverflowCheck(ConstRefPPMonoidElem pp, long exp)
  {
    if (exp < 0)  CoCoA_THROW_ERROR1(ERR::NegExp);
    owner(pp)->myPowerOverflowCheck(raw(pp), exp);
  }


  PPMonoidElem root(ConstRefPPMonoidElem pp, const MachineInt& exp)
  {
    if (IsNegative(exp) || IsZero(exp))  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    const long n = AsSignedLong(exp);
    if (n == 1) return pp;
    vector<long> exps = exponents(pp);
    const long nvars = NumIndets(owner(pp));
    for (long i=0; i < nvars; ++i)
    {
      if (exps[i]%n != 0) CoCoA_THROW_ERROR2(ERR::BadArg, "non an nth-power");
      exps[i] /= n;
    }
    return PPMonoidElem(owner(pp), exps);
  }


  bool IsPower(ConstRefPPMonoidElem pp, const MachineInt& exp)
  {
    if (IsNegative(exp))  CoCoA_THROW_ERROR1(ERR::NegExp);
    if (IsZero(exp)) return IsOne(pp);
    const long n = AsSignedLong(exp);
    if (n == 1) return true;
    //    vector<long> exps;    exponents(exps, pp);
    vector<long> exps = exponents(pp);
    const long nvars = NumIndets(owner(pp));
    for (long i=0; i < nvars; ++i)
      if (exps[i]%n != 0) return false;
    return true;
  }


  bool IsCoprime(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)  // are pp1, pp2 coprime?
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myIsCoprime(raw(pp1), raw(pp2));
  }


  bool IsDivisible(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)
  {
    CheckCompatible(pp1, pp2, CoCoA_ERROR_CONTEXT);
    return owner(pp1)->myIsDivisible(raw(pp1), raw(pp2));
  }


  bool IsSqFree(ConstRefPPMonoidElem pp)
  {
    return owner(pp)->myIsSqFree(raw(pp));
  }


//////////////////////////////////////////////////////////////////
// Next fns are for IsFactorClosed

  namespace // anonymous for local fns
  {
    bool CheckPredecessors(const std::set< std::vector<long> >& PrevDeg, const std::set< std::vector<long> >& CurrDeg)
    {
      CoCoA_ASSERT(!CurrDeg.empty());
      typedef set< vector<long> > SetPP;

      const int n = CurrDeg.begin()->size();
      vector<long> predecessor;
      for (SetPP::const_iterator it=CurrDeg.begin(); it != CurrDeg.end(); ++it)
      {
        predecessor = *it;
        for (int i=0; i < n; ++i)
        {
          if (predecessor[i] == 0) continue;
          --predecessor[i];
          if (PrevDeg.find(predecessor) == PrevDeg.end()) return false;
          ++predecessor[i];
        }
      }
      return true;
    }

    long deg(const vector<long>& a)
    {
      return accumulate(a.begin(), a.end(), 0l);
    }


    bool DegLex(const vector<long>& a, const vector<long>& b)
    {
      const long dega = deg(a);
      const long degb = deg(b);
      if (dega < degb) return true;
      if (dega > degb) return false;
      return lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
    }
  } // end of anonymous namespace


  bool IsFactorClosed(const std::vector<PPMonoidElem>& SetPP)
  {
    if (SetPP.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    PPMonoid PPM = owner(SetPP.front());
    const long NumPPs = len(SetPP);
    vector<long> expv(NumIndets(PPM));
    vector< vector<long> > S; S.reserve(NumPPs);
    for (long i=0; i < NumPPs; ++i)
      S.push_back(exponents(expv, SetPP[i]));
    sort(S.begin(), S.end(), DegLex);
    S.erase(unique(S.begin(), S.end()), S.end());
    if (deg(S.front()) != 0) return false;
    if (deg(S.back()) <= 1) return true;
    set< vector<long> > PrevDeg;
    long i = 1;
    long d = 1;
    while (i < NumPPs && deg(S[i]) == d)
    {
      PrevDeg.insert(S[i]);
      ++i;
    }

    while (i < NumPPs)
    {
      ++d;
      set< vector<long> > CurrDeg;
      while (i < NumPPs && deg(S[i]) == d)
      {
        CurrDeg.insert(S[i]);
        ++i;
      }

      if (!CheckPredecessors(PrevDeg, CurrDeg))
        return false;
      swap(PrevDeg, CurrDeg);
    }
    return true;
  }


  void AssignOne(RefPPMonoidElem dest)
  {
    owner(dest)->myAssignOne(raw(dest));
  }


  RefPPMonoidElem operator*=(RefPPMonoidElem pp, ConstRefPPMonoidElem pp1)
  {
    CheckCompatible(pp, pp1, CoCoA_ERROR_CONTEXT);
    owner(pp)->myMul(raw(pp), raw(pp), raw(pp1));
    return pp;
  }


  RefPPMonoidElem operator/=(RefPPMonoidElem pp, ConstRefPPMonoidElem pp1)
  {
    CheckCompatible(pp, pp1, CoCoA_ERROR_CONTEXT);
    owner(pp)->myDiv(raw(pp), raw(pp), raw(pp1));
    return pp;
  }


  const PPMonoidElem& indet(const PPMonoid& M, long index)
  {
    if (index >= NumIndets(M))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    return indets(M)[index];
  }


  const PPMonoidElem& indet(const PPMonoid& M, const BigInt& index)
  {
    CoCoA_ASSERT(index >= 0);
    long i;
    if (!IsConvertible(i, index))  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return indet(M, i);
  }


  PPMonoidElem IndetPower(const PPMonoid& M, long indet, long exp)
  {
    if (exp < 0)  CoCoA_THROW_ERROR1(ERR::NegExp);
    if (indet >= NumIndets(M))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    PPMonoidElem ans(M);
    M->myMulIndetPower(raw(ans), indet, exp);
    return ans;
  }

  PPMonoidElem IndetPower(const PPMonoid& M, long indet, const BigInt& EXP)
  {
    if (EXP < 0)  CoCoA_THROW_ERROR1(ERR::NegExp);
    if (indet >= NumIndets(M))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    PPMonoidElem ans(M);
    M->myMulIndetPower(raw(ans), indet, EXP);
    return ans;
  }


//   // Possible default implementation -- would wastefully alloc/free though
//   void PPMonoidElem::myComputeDivMask(DivMask::value& dm, const DivMask::base& DivMaskImpl, const RawPtr& pp) const
//   {
//     std::vector<long> v(myNumIndets);

//     myExponents(v, pp);
//     DivMaskImpl.myAssignFromExpv(dm, &v[0], myNumIndets);
//   }

// ------------------------------ Input/Output ------------------------------

  std::ostream& operator<<(std::ostream& out, ConstRefPPMonoidElem pp)
  {
    if (!out) return out;  // short-cut for bad ostreams
    owner(pp)->myOutput(out, raw(pp));
    return out;
  }


  OpenMathOutput& operator<<(OpenMathOutput& OMOut, ConstRefPPMonoidElem pp)
  {
    owner(pp)->myOutput_OM(OMOut, raw(pp));
    return OMOut;
  }


  PPMonoid NewPPMonoid(const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
    return NewPPMonoidEvOv(IndetNames, ord); // use EvOv by default
  }

  PPMonoid NewPPMonoid(const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  {
    return NewPPMonoidEvOv(IndetNames, OrdCtor); // use EvOv by default
  }


} // end of namespace CoCoA
