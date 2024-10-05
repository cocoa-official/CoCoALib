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
    inline void CheckCompatible(const ConstRefPPMonoidElem& pp1, const ConstRefPPMonoidElem& pp2, const char* const FnName)
    {
      if (owner(pp1) != owner(pp2))
        CoCoA_THROW_ERROR(ERR::MixedPPMs, FnName);
    }
  }


  inline PPMonoidElemRawPtr PPMonoidBase::myNewCheckVecSize(const std::vector<long>& v) const
  {
    if (len(v) != myNumIndets)
      CoCoA_THROW_ERROR(ERR::BadArraySize, "PPMonoidElem(PPM, ExpVec)  or  PPMonoidBase::myNewCheckVecSize");
    return myNew(v);
  }

  inline PPMonoidElemRawPtr PPMonoidBase::myNewCheckVecSize(const std::vector<BigInt>& v) const
  {
    if (len(v) != myNumIndets)
      CoCoA_THROW_ERROR(ERR::BadArraySize, "PPMonoidElem(PPM, ExpVec)  or  PPMonoidBase::myNewCheckVecSize");
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
    if (!IsConvertible(exp, EXP))
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "IndetPower/myMulIndetPower");
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
  {
    CoCoA_THROW_ERROR(ERR::ExpTooBig, "power(pp,N)");
  }


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
    if (myIsOne(rawpp)) return false; // changed in 0.99700 (previos gave error)
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
    if (idx < 0 || idx >= NumIndets(PPM)) CoCoA_THROW_ERROR(ERR::BadIndex, "IndetSymbol");
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
    CoCoA_THROW_ERROR("unknown symbol in PPM", "PPMonoid::mySymbolValue");
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
    CheckCompatible(*this, rhs, "PPMonoidElem::operator=");
    // if (this == &rhs) return *this;  // trivial check not really needed
    myPPM->myAssign(myPPPtr, raw(rhs));
    return *this;
  }


  RefPPMonoidElem& RefPPMonoidElem::operator=(const ConstRefPPMonoidElem& rhs)
  {
    CheckCompatible(*this, rhs, "PPMonoidElem::operator=");
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



  // This function assumes onwership of the value pointed to by ToBeOwned
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
    CheckCompatible(*this, rhs, "PPMonoidElem::operator=");
    // if (this == &rhs) return *this;  // trivial check not really needed
    myPPM->myAssign(myPPPtr, raw(rhs));
    return *this;
  }

  PPMonoidElem& PPMonoidElem::operator=(const ConstRefPPMonoidElem& rhs)
  {
    CheckCompatible(*this, rhs, "PPMonoidElem::operator=");
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
    CheckCompatible(pp1, pp2, "CmpWDeg(PP,PP)");
    return owner(pp1)->myCmpWDeg(raw(pp1), raw(pp2));
  }


  int CmpWDegPartial(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2, long PartialGradingDim)
  {
    CheckCompatible(pp1, pp2, "CmpWDeg(PP,PP)");
    if ( PartialGradingDim < 0 || PartialGradingDim > GradingDim(owner(pp1)) )
      CoCoA_THROW_ERROR("PartialGradingDim > GradingDim", "CmpWDegPartial");
    return owner(pp1)->myCmpWDegPartial(raw(pp1), raw(pp2), PartialGradingDim);
  }


  long exponent(ConstRefPPMonoidElem pp, long indet) // degree in indet
  {
    if (indet >= NumIndets(owner(pp)))
      CoCoA_THROW_ERROR(ERR::BadIndex, "exponent(PP)");
    return owner(pp)->myExponent(raw(pp), indet);
  }


  BigInt BigExponent(ConstRefPPMonoidElem pp, long indet)
  {
    if (indet >= NumIndets(owner(pp)))
      CoCoA_THROW_ERROR(ERR::BadIndex, "BigExponent(PP)");
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
    CheckCompatible(pp1, pp2, "swap(PP,PP)");
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
    CheckCompatible(pp1, pp2, "cmp(PP,PP)");
    return owner(pp1)->myCmp(raw(pp1), raw(pp2));
  }


  bool operator==(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 == pp2
  {
    CheckCompatible(pp1, pp2, "PP == PP");
    return owner(pp1)->myIsEqual(raw(pp1), raw(pp2));
  }


  bool operator!=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 != pp2
  {
    CheckCompatible(pp1, pp2, "PP != PP");
    return !owner(pp1)->myIsEqual(raw(pp1), raw(pp2));
  }


  bool operator<(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)           // pp1 < pp2
  {
    CheckCompatible(pp1, pp2, "PP < PP");
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) < 0;
  }


  bool operator<=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 <= pp2
  {
    CheckCompatible(pp1, pp2, "PP <= PP");
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) <= 0;
  }


  bool operator>(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)           // pp1 > pp2
  {
    CheckCompatible(pp1, pp2, "PP > PP");
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) > 0;
  }


  bool operator>=(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)          // pp1 >= pp2
  {
    CheckCompatible(pp1, pp2, "PP >= PP");
    return owner(pp1)->myCmp(raw(pp1), raw(pp2)) >= 0;
  }


// ------------------------------ Arithmetic ------------------------------

  PPMonoidElem operator*(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)   //  pp1*pp2;
  {
    CheckCompatible(pp1, pp2, "PP * PP");
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myMul(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem operator/(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)   // pp1/pp2;
  {
    CheckCompatible(pp1, pp2, "PP / PP");
    if (!owner(pp1)->myIsDivisible(raw(pp1), raw(pp2)))
      CoCoA_THROW_ERROR(ERR::BadQuot, "PP / PP");
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myDiv(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem colon(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)       //  pp1:pp2;
  {
    CheckCompatible(pp1, pp2, "colon(PP, PP)");
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myColon(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem gcd(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)         // gcd(pp1,pp2);
  {
    CheckCompatible(pp1, pp2, "gcd(PP, PP)");
    PPMonoidElem ans(owner(pp1));
    owner(pp1)->myGcd(raw(ans), raw(pp1), raw(pp2));
    return ans;
  }


  PPMonoidElem lcm(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)         // lcm(pp1,pp2);
  {
    CheckCompatible(pp1, pp2, "lcm(PP, PP)");
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
    if (exp < 0)
      CoCoA_THROW_ERROR(ERR::NegExp, "power(pp, n) and n < 0");
    PPMonoidElem ans(owner(pp));
    owner(pp)->myPower(raw(ans), raw(pp), exp);
    return ans;
  }


  PPMonoidElem power(ConstRefPPMonoidElem pp, const BigInt& EXP)                 // pp^EXP
  {
    if (EXP < 0)
      CoCoA_THROW_ERROR(ERR::NegExp, "power(pp, N) and N < 0");
    PPMonoidElem ans(owner(pp));
    owner(pp)->myPower(raw(ans), raw(pp), EXP);
    return ans;
  }


  void PowerOverflowCheck(ConstRefPPMonoidElem pp, long exp)
  {
    if (exp < 0)
      CoCoA_THROW_ERROR(ERR::NegExp, "PowerOverflowCheck(pp, n) and n < 0");
    owner(pp)->myPowerOverflowCheck(raw(pp), exp);
  }


  PPMonoidElem root(ConstRefPPMonoidElem pp, const MachineInt& exp)
  {
    if (IsNegative(exp) || IsZero(exp)) CoCoA_THROW_ERROR(ERR::ReqPositive, "root(PP, exp)");
    const long n = AsSignedLong(exp);
    if (n == 1) return pp;
    vector<long> exps = exponents(pp);
    const long nvars = NumIndets(owner(pp));
    for (long i=0; i < nvars; ++i)
    {
      if (exps[i]%n != 0) CoCoA_THROW_ERROR(ERR::BadArg, "root(PP,exp)");
      exps[i] /= n;
    }
    return PPMonoidElem(owner(pp), exps);
  }


  bool IsPower(ConstRefPPMonoidElem pp, const MachineInt& exp)
  {
    if (IsNegative(exp)) CoCoA_THROW_ERROR(ERR::ReqNonNegative, "IsPower(PP, exp)");
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
    CheckCompatible(pp1, pp2, "IsCoprime(PP, PP)");
    return owner(pp1)->myIsCoprime(raw(pp1), raw(pp2));
  }


  bool IsDivisible(ConstRefPPMonoidElem pp1, ConstRefPPMonoidElem pp2)
  {
    CheckCompatible(pp1, pp2, "IsDivisible(PP, PP)");
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
    if (SetPP.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
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
    CheckCompatible(pp, pp1, "PP *= PP");
    owner(pp)->myMul(raw(pp), raw(pp), raw(pp1));
    return pp;
  }


  RefPPMonoidElem operator/=(RefPPMonoidElem pp, ConstRefPPMonoidElem pp1)
  {
    CheckCompatible(pp, pp1, "PP /= PP");
    owner(pp)->myDiv(raw(pp), raw(pp), raw(pp1));
    return pp;
  }


  const PPMonoidElem& indet(const PPMonoid& M, long index)
  {
    if (index >= NumIndets(M)) CoCoA_THROW_ERROR(ERR::BadIndex, "indet(PPMonoid, index)");
    return indets(M)[index];
  }


  const PPMonoidElem& indet(const PPMonoid& M, const BigInt& index)
  {
    CoCoA_ASSERT(index >= 0);
    long i;
    if (!IsConvertible(i, index)) CoCoA_THROW_ERROR(ERR::ArgTooBig, "indet(PPMonoid, index)");
    return indet(M, i);
  }


  PPMonoidElem IndetPower(const PPMonoid& M, long indet, long exp)
  {
    if (exp < 0) CoCoA_THROW_ERROR(ERR::NegExp, "IndetPower(PPMonoid, indet, exp)");
    if (indet >= NumIndets(M)) CoCoA_THROW_ERROR(ERR::BadIndex, "IndetPower(PPMonoid, indet, exp)");
    PPMonoidElem ans(M);
    M->myMulIndetPower(raw(ans), indet, exp);
    return ans;
  }

  PPMonoidElem IndetPower(const PPMonoid& M, long indet, const BigInt& EXP)
  {
    if (EXP < 0) CoCoA_THROW_ERROR(ERR::NegExp, "IndetPower(PPMonoid, indet, EXP)");
    if (indet >= NumIndets(M)) CoCoA_THROW_ERROR(ERR::BadIndex, "IndetPower(PPMonoid, indet, EXP)");
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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPMonoid.C,v 1.55 2024/07/15 17:10:00 bigatti Exp $
// $Log: PPMonoid.C,v $
// Revision 1.55  2024/07/15 17:10:00  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.54  2024/07/03 15:35:22  bigatti
// Summary: Specialized some BadArg with explicit codes
//
// Revision 1.53  2024/07/03 09:35:13  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.52  2022/10/11 15:27:41  bigatti
// Summary: using new interface of exponents returning vector
//
// Revision 1.51  2022/09/14 19:15:04  abbott
// Summary: Revised interface to exponents (redmine 1694)
//
// Revision 1.50  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.49  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.48  2021/02/17 17:50:36  abbott
// Summary: Now asserts ostream is in decimal mode (redmine 1547)
//
// Revision 1.47  2020/06/17 15:49:24  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.46  2020/02/15 11:10:22  abbott
// Summary: IsIndetPosPower now gives false instead of error for input 1
//
// Revision 1.45  2020/02/14 12:21:51  abbott
// Summary: Renames indets to IndetsIn
//
// Revision 1.44  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.43  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.42  2019/11/14 17:51:32  abbott
// Summary: Removed some commented out cruft
//
// Revision 1.41  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.40  2017/12/01 17:16:16  bigatti
// -- just Copyright
//
// Revision 1.39  2016/11/11 14:15:32  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.38  2016/11/03 12:25:24  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.37  2015/12/01 13:11:00  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.36  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.35  2015/04/13 14:42:07  abbott
// Summary: Added myPowerOverflowCheck (1st version)
// Author: JAA
//
// Revision 1.34  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.33  2014/07/03 15:36:35  abbott
// Summary: Cleaned up impl of PPMonoids: moved myIndetSymbols & myNumIndets to base class
// Author: JAA
//
// Revision 1.32  2014/06/17 10:11:34  abbott
// Summary: Changed name of a param
// Author: JAA
//
// Revision 1.31  2014/03/05 13:46:11  abbott
// Summary: IsFactorClosed now gives error when called on empty list
// Author: JAA
//
// Revision 1.30  2014/02/25 16:27:57  abbott
// Summary: Added new fn IsFactorClosed
// Author: JAA
//
// Revision 1.29  2012/10/22 16:48:56  bigatti
// -- modified IsIndetPosPower so that N unchanged when result is false
//
// Revision 1.28  2012/05/29 15:26:45  bigatti
// -- fixed capitalization of IntOperations.H
//
// Revision 1.27  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.26  2012/04/17 19:55:40  abbott
// Added 2 new fns: root and IsPower.
//
// Revision 1.25  2012/02/10 17:05:54  abbott
// Added default defn of pseudo-ctor myNew(vector<BigInt>).
// Added new fn indets.
//
// Revision 1.24  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.23  2011/03/22 16:44:39  bigatti
// -- removed debugging printout
//
// Revision 1.22  2011/03/22 15:26:46  bigatti
// -- added myIsIndetPosPower/IsIndetPosPower
//
// Revision 1.21  2011/03/16 15:19:25  bigatti
// -- added IsIndet(pp), IsIndetPosPower(pp)
//
// Revision 1.20  2011/03/16 13:23:17  abbott
// Corrected arg sanity check in CmpWDegPartial.
//
// Revision 1.19  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.18  2010/11/30 11:18:11  bigatti
// -- renamed IndetName --> IndetSymbol
//
// Revision 1.17  2010/11/05 16:21:08  bigatti
// -- added ZZExponents
//
// Revision 1.16  2010/11/02 16:03:12  bigatti
// -- added indet(*, ZZ) [for CoCoA-5]
// -- indet(PPM, i) now returns "const PPMonoidElem&"
//
// Revision 1.15  2010/10/01 15:45:17  bigatti
// -- added mySymbolValue
//
// Revision 1.14  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.13  2010/02/03 14:12:34  bigatti
// -- added IsIndetPosPower with exponent of type long
//
// Revision 1.12  2010/02/02 16:44:31  abbott
// Added radical & IsRadical (via mem fns myRadical & myIsRadical)
// for PPMonoidElems.
//
// Revision 1.11  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.10  2009/10/29 18:31:50  abbott
// Changed order of include directives (now alphabetical by file name).
//
// Revision 1.9  2009/09/22 14:01:33  bigatti
// -- added myCmpWDegPartial (ugly name, I know....)
// -- cleaned up and realigned code in PPMonoid*.C files
//
// Revision 1.8  2008/10/07 12:23:43  abbott
// Added deg function (same as StdDeg).
//
// Revision 1.7  2008/09/16 14:03:39  bigatti
// -- added comment
//
// Revision 1.6  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.5  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/24 14:21:29  abbott
// Added default defn of myMulIndetPower.
// Added IndetPower.
// Renamed myIsIndetPower to myIsIndetPosPower.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.8  2007/03/07 13:44:49  bigatti
// -- minor cleanup for -Wextra
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/06 17:34:51  cocoa
// -- style: RawPtr args are now called "raw.."
//
// Revision 1.5  2006/11/24 17:04:32  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/11/23 17:35:48  cocoa
// -- changed: PPMonoid is now a class (instead of typedef)
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.10  2006/04/21 14:59:22  cocoa
// Changed return type of indet(PPM, i) to ConstRefPPMonoidElem
// instead of PPMonoidElem (which made a copy).
//
// Revision 1.9  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.8  2006/03/14 17:21:18  cocoa
// Moved concrete PPMonoid impls entirely into their respective .C files.
// Now the corresponding .H files are very compact.
//
// Revision 1.7  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.6  2006/03/07 10:08:48  cocoa
// -- fixed: PPMonoidElem(const ConstRefPPMonoidElem& copy);  [added const]
//
// Revision 1.5  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.4  2006/02/13 12:08:04  cocoa
// Fixed a problem with some missing assignment ops for certain PPMonoidElems.
// Fixed a bug in RingDistrMPoly::myIndetPower.
//
// Revision 1.3  2006/02/01 16:56:13  cocoa
// Added some missing assignment operators for (Ref)PPMonoidElems.
//
// Revision 1.2  2005/11/17 15:43:04  cocoa
// -- added ctor PPMonoidElem(ConstRefPPMonoidElem& copy);
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.7  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.6  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.5  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.4  2005/06/27 14:55:24  cocoa
// Cleaned up some more PPMonoid code.
//
// Revision 1.3  2005/06/23 15:42:41  cocoa
// Fixed typo in GNU fdl -- all doc/*.txt files affected.
// Minor corrections to PPMonoid (discovered while writing doc).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.3  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.11  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.10  2004/11/11 13:40:39  cocoa
// -- minor changes for doxygen
//
// Revision 1.9  2004/11/02 14:56:33  cocoa
// -- changed *Print* into *Output* (myPrint --> myOutput)
// -- changed *Var* into *Indet* (myPrintVarName --> myOutputIndetName)
// -- removed suffix "IgnoreDivMask"
// -- added myComputeDivMask
// -- improved storing of IndetNames
// -- changed ExpvElem into SmallExponent_t
//
// Revision 1.8  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.7  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.6  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.5  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.4  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.3  2004/01/28 15:20:37  cocoa
// Added CmpDeg function and myCmpDeg member function.
// Aligned name of data member of PPMonoid with the convention used
// for the data member of ring.
//
// Revision 1.2  2003/10/01 10:35:31  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.17  2003/09/23 17:01:51  bigatti
// - added #include "CoCoA/PPMonoidDiv.H"
// - changed default monoid to PPMonoidDiv
//
// Revision 1.16  2003/09/22 13:04:24  abbott
// Improved definition of operator!= and tidied a few comments.
// Something odd had happened: the corrected version was write
// protected but not checked in.
//
// Revision 1.15  2003/06/23 16:57:08  abbott
// Minor cleaning prior to public release.
//
// Revision 1.14  2003/05/14 16:41:41  abbott
// Added new functions: IsOne and PPMonoidBase::IsEqualIgnoreDivMask
// (default definition of a virtual functions).
// Revitalized the function exponents.
//
// Revision 1.13  2003/04/29 13:59:37  abbott
// Uh oh.  Lots of changes.  Checking in prior to restructuring
// to align with new PPOrdering code.
//
// Revision 1.12  2003/01/15 16:03:50  abbott
// Checking in prior to major change.
// Minor changes:
//   new definition of deg (to go with new signature),
//   better error messages in fns on PPMonoid::elems,
//   better exception safety in fns on PPMonoid::elems.
//
// Revision 1.11  2002/11/15 17:07:21  bigatti
// - added check of myDivMask in operators == and !=
//
// Revision 1.10  2002/11/15 16:25:41  abbott
// MAJOR CHANGE reflecting a new "philosophy".
// The new philosophy is that PPMonoid::elems will be "slow and safe" since
// they are now higly decoupled from the internal representation of
// polynomials (which use "order vectors").
// Eliminated HalfPP, alias and temp subclasses.
// Every RawPP now contains its own "divisibility mask" because Anna needs
// the DivMask check to be inline.
//
// Revision 1.9  2002/06/22 17:16:16  abbott
// Numerous changes to accommodate the new idea of PPOrdering.
// And some other changes (e.g. NumVars is now NumIndets).
//
// Revision 1.8  2002/03/06 11:24:19  abbott
// Un-inlined the destructor for PPMonoid::elem (but shouldn't
// it be inline in the header file?).
//
// Revision 1.7  2002/02/08 12:07:48  abbott
// MAJOR CHANGE: added two new proxy classes "alias" and "temp" (see PP.txt)
// Numerous consequential changes.  Changed order of some definitions.
//
// Revision 1.6  2002/01/30 14:21:06  abbott
// Tidied up use of "using" declaration.
// Added definition of PPMonoid::print -- all concrete derived classes had
// effectively the same implementation.  This is the only non-pure virtual
// function in the class PPMonoid.
//
// Revision 1.5  2001/12/04 19:53:32  abbott
// Names changed in accordance with the new coding conventions.
//
// Revision 1.4  2001/11/19 15:25:47  abbott
// Corrected an imcompatibility between PPmonoid::cmp and
// bool operator<(const PPmonoid::elem&, const PPmonoid::elem&).
// The latter expected cmp to return one of -1,0 or +1 whereas
// cmp actually returned negative, zero, or positive.
//
// Revision 1.3  2001/11/16 18:41:58  bigatti
// added:   using namespace std;
// for compatibility with gcc-3
//
// Revision 1.2  2001/10/31 20:42:57  abbott
// Pseudo copy constructor for PPmonoid::elem now takes a const copy as
// argument (previous const was omitted).
//
// Revision 1.1  2001/10/25 17:54:55  abbott
// Initial revision
//

