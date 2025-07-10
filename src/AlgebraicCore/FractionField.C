//   Copyright (c)  2005-2009,2014,2021  John Abbott and Anna M. Bigatti

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


// Source code for classes FractionField & FractionFieldImpl

#include "CoCoA/FractionField.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/GlobalManager.H" // needed only by NewFractionField
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"         // needed only by NewFractionField
#include "CoCoA/RingZZ.H"         // needed only by NewFractionField
#include "CoCoA/assert.H"
///#include "CoCoA/error.H"
#include "CoCoA/ideal.H"

#include <memory>
using std::unique_ptr;
#include <iostream>
using std::ostream;
#include <new>
//for placement new


namespace CoCoA
{

  const FractionFieldBase* FractionFieldPtr(const ring& R)
  {
    return dynamic_cast<const FractionFieldBase*>(R.myRawPtr());
  }

  // const FractionFieldBase* FractionFieldPtr(const ring& R, const ErrorContext& FnName)
  // {
  //   const FractionFieldBase* ptr = FractionFieldPtr(R);
  //   if (ptr == nullptr)
  //     CoCoA_THROW_ERROR1(ERR::NotFracField, FnName);
  //   return ptr;
  // }

  const FractionFieldBase* FractionFieldPtr(const ring& R, const ErrorContext& ErrCtx)
  {
    const FractionFieldBase* ptr = FractionFieldPtr(R);
    if (ptr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::NotFracField, ErrCtx);
    return ptr;
  }


  BigInt FractionFieldBase::myCharacteristic() const
  {
    return characteristic(myBaseRingValue);
///???    return myBaseRingValue->myCharacteristic();
  }


  bool FractionFieldBase::IamCommutative() const
  {
    return true;
  }


  bool3 FractionFieldBase::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool FractionFieldBase::IamOrderedDomain() const
  {
    return myBaseRingValue->IamOrderedDomain();
  }


  bool FractionFieldBase::IamField() const
  {
    return true;
  }


  bool FractionFieldBase::IamFiniteField() const
  {
    return false;
  }


  bool FractionFieldBase::IamExact() const
  {
    return IsExact(myBaseRing());
  }



  // This ctor is not inline since it probably won't be called very much.
  FractionField::FractionField(const FractionFieldBase* RingPtr):
      ring(RingPtr)
  {}


  RingElem FractionFieldBase::mySymbolValue(const symbol& sym) const
  {
    return myEmbeddingHomCtor()(myBaseRing()->mySymbolValue(sym));
  }


  // This fn is called only by DerivFrF (see file PolyRing.C)
  void FractionFieldBase::myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(IsPolyRing(myBaseRing()));
    const PolyRing Rx = myBaseRing();
    CoCoA_ASSERT(Rx->myIsOne(myRawDen(rawx)));
    const RingElemAlias x(Rx, myRawNum(rawx));
    if (!IsIndet(x))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
    // Code below *ASSUMES* x is an indet (and not a power product)
    const RingElemAlias N(Rx, myRawNum(rawf));
    const RingElemAlias D(Rx, myRawDen(rawf));
    const RingElem derivN = deriv(N, x);
    const RingElem derivD = deriv(D, x);
    const RingHom phi = myEmbeddingHomCtor();
    const RingElem ans = phi(derivN*D - N*derivD)/phi(D*D);
    myAssign(rawlhs, raw(ans)); ///??? swap???
  }


  //----------------------------------------------------------------------
  // Below is a generic implementation of a fraction field.

  class FractionFieldImpl; // fwd decl for friend decl
  // This class defines the data representation used for elements of a FractionFieldImpl.
  // Instantiations of this class are not true objects: they cannot delete themselves.
  class FractionFieldElem
  {
  public:
    FractionFieldElem(RingElemRawPtr rawN, RingElemRawPtr rawD): myNumerator(rawN), myDenominator(rawD) {}
  private: // data members
    RingElemRawPtr myNumerator;
    RingElemRawPtr myDenominator;
    friend class FractionFieldImpl;
  };


  class FractionFieldImpl: public FractionFieldBase
  {
  private: // data members
    mutable MemPool myMemMgr;           // MemPool must come before myZeroPtr and myOnePtr
    unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.

  private:
    typedef FractionFieldElem value_t; // FractionFieldElem is the actual type of the values in a FractionFieldImpl
    static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(RingElemConstRawPtr rawx);

    static RawPtr RefNum(RawPtr rawq);           // result belongs to BaseRing
    static RawPtr RefDen(RawPtr rawq);           // result belongs to BaseRing
    static ConstRawPtr RefNum(ConstRawPtr rawq); // result belongs to BaseRing
    static ConstRawPtr RefDen(ConstRawPtr rawq); // result belongs to BaseRing
    void CancelFactors(RawPtr rawx) const;

  private:
    friend FractionField NewFractionField(const ring& R); // the only fn which calls the ctor
    FractionFieldImpl(const ring& R);
    ~FractionFieldImpl();

  public: // functions that every ring must implement (see also FractionFieldBase)
    ConstRefRingElem myZero() const override;
    ConstRefRingElem myOne() const override;
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& N) const override;
    RingElemRawPtr myNew(ConstRawPtr rawt) const override;
    void myDelete(RawPtr rawx) const override;                        // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;             // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;    // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override; // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;     // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                  // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;    // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; // lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; // lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                          // true iff x is invertible
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // lhs = gcd(x,y) in a field
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;// lhs = x^n, n>1, x not -1,0,1
    void mySymbols(std::vector<symbol>& SymList) const override;                   // appends ring's symbols to SymList
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;             // out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;
    void myOutputSelf(std::ostream& out) const override;                           // out << R
    void myOutputSelfLong(std::ostream& out) const override; // out << R (descr)
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                       // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;         // OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                                // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                 // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                            // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                  // true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                     // true iff x is rational
//???    bool myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const;   // lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;                // x == y

    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))

    bool myImageLiesInSubfield(const RingHom& phi) const override;

    // functions every FractionField must implement
    ConstRawPtr myRawNum(ConstRawPtr rawq) const override;  ///< result belongs to BaseRing!!
    ConstRawPtr myRawDen(ConstRawPtr rawq) const override;  ///< result belongs to BaseRing!!
    RingHom myEmbeddingHomCtor() const override;
    RingHom myInducedHomCtor(const RingHom& phi) const override;

  private:
    class InducedHomImpl: public RingHomInducedBase
    {
    public:
      InducedHomImpl(const FractionField& FrF, const RingHom& InducingHom);
      // Default copy ctor & assignment disabled in RingHomBase.
      // Default dtor works fine.
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return IsPartial(myInducingHom) || !IsField(myCodomain); } /// BUG not strictly correct, image may be a subfield of myCodomain
    };

  private:
    class EmbeddingHomImpl: public RingHomEmbeddingBase
    {
    public:
      EmbeddingHomImpl(const FractionField& FrF);
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return false; }
    };
  };


  //----------------------------------------------------------------------

  inline FractionFieldImpl::value_t& FractionFieldImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const FractionFieldImpl::value_t& FractionFieldImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }


  // Four handy accessor functions -- really just a shorthand.
  // These are more useful than the "import" functions above.
  inline RingElemRawPtr FractionFieldImpl::RefNum(RawPtr rawq)
  {
    return import(rawq).myNumerator;
  }


  inline RingElemRawPtr FractionFieldImpl::RefDen(RawPtr rawq)
  {
    return import(rawq).myDenominator;
  }


  inline RingElemConstRawPtr FractionFieldImpl::RefNum(ConstRawPtr rawq)
  {
    return import(rawq).myNumerator;
  }


  inline RingElemConstRawPtr FractionFieldImpl::RefDen(ConstRawPtr rawq)
  {
    return import(rawq).myDenominator;
  }



  // This ctor is called only from NewFractionField.
  FractionFieldImpl::FractionFieldImpl(const ring& R):
    FractionFieldBase(R),
    myMemMgr(sizeof(FractionFieldElem), "FractionFieldImpl.myMemMgr")
  {
    // NewFractionField has already checked that R is commutative GCD domain, not a field
    CoCoA_ASSERT(IsCommutative(R) && IsTrueGCDDomain(R));
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myRefCountZero();
  }


  FractionFieldImpl::~FractionFieldImpl()
  {}


  ConstRefRingElem FractionFieldImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem FractionFieldImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemConstRawPtr FractionFieldImpl::myRawNum(ConstRawPtr rawq) const
  {
    return RefNum(rawq);
  }

  RingElemConstRawPtr FractionFieldImpl::myRawDen(ConstRawPtr rawq) const
  {
    return RefDen(rawq);
  }


  /////////////////////////////////////////////////////////////////////////////

  RingElemRawPtr FractionFieldImpl::myNew() const
  {
    AutoRingElem num(myBaseRingValue, myBaseRingValue->myNew());  //
    AutoRingElem den(myBaseRingValue, myBaseRingValue->myNew(1)); // Any of these 3 lines could throw.
    void* ptr = myMemMgr.alloc();         //
    new(ptr) FractionFieldElem(release(num), release(den));
    return RingElemRawPtr(ptr);
//    return new(myMemMgr.alloc()) FractionFieldElem(release(N), release(D));
//    return new(myMemMgr.alloc()) FractionFieldElem(myBaseRingValue->myNew(), myBaseRingValue->myNew(1));
//     unique_ptr<value_t> ans = static_cast<value_t*>(myMemMgr.alloc());
//     ans->myNumerator = myBaseRingValue->myNew();
//     ans->myDenominator = myBaseRingValue->myNew(1);
//     return ans.release();
  }


  RingElemRawPtr FractionFieldImpl::myNew(const MachineInt& n) const
  {
    AutoRingElem num(myBaseRingValue, myBaseRingValue->myNew(n)); //
    AutoRingElem den(myBaseRingValue, myBaseRingValue->myNew(1)); // Any of these 3 lines could throw.
    void* ptr = myMemMgr.alloc();         //
    new(ptr) FractionFieldElem(release(num), release(den));
    return RingElemRawPtr(ptr);
//    return new(myMemMgr.alloc()) FractionFieldElem(myBaseRingValue->myNew(n), myBaseRingValue->myNew(1));
//     unique_ptr<value_t> ans = static_cast<value_t*>(myMemMgr.alloc());
//     ans->myNumerator = myBaseRingValue->myNew(n);
//     ans->myDenominator = myBaseRingValue->myNew(1);
//     return ans.release();
  }


  RingElemRawPtr FractionFieldImpl::myNew(const BigInt& N) const
  {
    AutoRingElem num(myBaseRingValue, myBaseRingValue->myNew(N)); //
    AutoRingElem den(myBaseRingValue, myBaseRingValue->myNew(1)); // Any of these 3 lines could throw.
    void* ptr = myMemMgr.alloc();         //
    new(ptr) FractionFieldElem(release(num), release(den));
    return RingElemRawPtr(ptr);
//    return new(myMemMgr.alloc()) FractionFieldElem(myBaseRingValue->myNew(N), myBaseRingValue->myNew(1));
//     unique_ptr<value_t> ans = static_cast<value_t*>(myMemMgr.alloc());
//     ans->myNumerator = myBaseRingValue->myNew(N);
//     ans->myDenominator = myBaseRingValue->myNew(1);
//     return ans.release();
  }


  RingElemRawPtr FractionFieldImpl::myNew(ConstRawPtr rawCopyMe) const
  {
    AutoRingElem num(myBaseRingValue, myBaseRingValue->myNew(RefNum(rawCopyMe))); //
    AutoRingElem den(myBaseRingValue, myBaseRingValue->myNew(RefDen(rawCopyMe))); // Any of these 3 lines could throw.
    void* ptr = myMemMgr.alloc();                      //
    new(ptr) FractionFieldElem(release(num), release(den));
    return RingElemRawPtr(ptr);
//    return new(myMemMgr.alloc()) FractionFieldElem(myBaseRingValue->myNew(RefNum(rawCopyMe)), myBaseRingValue->myNew(RefDen(rawCopyMe)));
//     unique_ptr<value_t> ans = static_cast<value_t*>(myMemMgr.alloc());
//     ans->myNumerator = myBaseRingValue->myNew(RefNum(rawCopyMe));
//     ans->myDenominator = myBaseRingValue->myNew(RefDen(rawCopyMe));
//     return ans.release();
  }


  void FractionFieldImpl::myDelete(RawPtr rawx) const
  {
    myBaseRingValue->myDelete(RefDen(rawx));
    myBaseRingValue->myDelete(RefNum(rawx));
    myMemMgr.free(rawx.myRawPtr());
  }


  void FractionFieldImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    myBaseRingValue->mySwap(RefNum(rawx), RefNum(rawy));
    myBaseRingValue->mySwap(RefDen(rawx), RefDen(rawy));
  }


  void FractionFieldImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myBaseRingValue->myAssign(RefNum(rawlhs), RefNum(rawx));
    myBaseRingValue->myAssign(RefDen(rawlhs), RefDen(rawx));
  }


  void FractionFieldImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    myBaseRingValue->myAssign(RefNum(rawlhs), n);
    myBaseRingValue->myAssign(RefDen(rawlhs), 1);
  }


  void FractionFieldImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    myBaseRingValue->myAssign(RefNum(rawlhs), N);
    myBaseRingValue->myAssign(RefDen(rawlhs), 1);
  }


  void FractionFieldImpl::myAssignZero(RawPtr rawlhs) const
  {
    myBaseRingValue->myAssignZero(RefNum(rawlhs));
    myBaseRingValue->myAssign(RefDen(rawlhs), 1);
  }


  void FractionFieldImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myBaseRingValue));
    myBaseRingValue->myRecvTwinFloat(RefNum(rawlhs), rawx);
    myBaseRingValue->myAssign(RefDen(rawlhs), 1);
  }


  void FractionFieldImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myBaseRingValue->myNegate(RefNum(rawlhs), RefNum(rawx));
    myBaseRingValue->myAssign(RefDen(rawlhs), RefDen(rawx));
  }


  void FractionFieldImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myBaseRingValue->myIsOne(RefDen(rawx)) &&
	myBaseRingValue->myIsOne(RefDen(rawy)))
    {
      myBaseRingValue->myAdd(RefNum(rawlhs), RefNum(rawx), RefNum(rawy));
      myBaseRingValue->myAssign(RefDen(rawlhs), 1);
      return;
    }
    RingElem prod1(myBaseRingValue);
    myBaseRingValue->myMul(raw(prod1), RefNum(rawx), RefDen(rawy));
    RingElem prod2(myBaseRingValue);
    myBaseRingValue->myMul(raw(prod2), RefNum(rawy), RefDen(rawx));
    myBaseRingValue->myMul(RefDen(rawlhs), RefDen(rawx), RefDen(rawy));
    myBaseRingValue->myAdd(RefNum(rawlhs), raw(prod1), raw(prod2));
    myBaseRingValue->myNormalizeFrac(RefNum(rawlhs),RefDen(rawlhs));
//    CancelFactors(rawlhs);
  }


  void FractionFieldImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myBaseRingValue->myIsOne(RefDen(rawx)) &&
	myBaseRingValue->myIsOne(RefDen(rawy)))
    {
      myBaseRingValue->mySub(RefNum(rawlhs), RefNum(rawx), RefNum(rawy));
      myBaseRingValue->myAssign(RefDen(rawlhs), 1);
      return;
    }
    RingElem prod1(myBaseRingValue);
    myBaseRingValue->myMul(raw(prod1), RefNum(rawx), RefDen(rawy));
    RingElem prod2(myBaseRingValue);
    myBaseRingValue->myMul(raw(prod2), RefNum(rawy), RefDen(rawx));
    myBaseRingValue->myMul(RefDen(rawlhs), RefDen(rawx), RefDen(rawy));
    myBaseRingValue->mySub(RefNum(rawlhs), raw(prod1), raw(prod2));
    myBaseRingValue->myNormalizeFrac(RefNum(rawlhs),RefDen(rawlhs));
//    CancelFactors(rawlhs);
  }


  void FractionFieldImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    // A "sophisticated" implementation which does "minimal" gcd computation
    RingElem g1(myBaseRingValue), g2(myBaseRingValue), q1(myBaseRingValue), q2(myBaseRingValue);
    myBaseRingValue->myGcd(raw(g1), RefNum(rawx), RefDen(rawy));   // g1 = gcd(num(x), den(y))
    myBaseRingValue->myGcd(raw(g2), RefNum(rawy), RefDen(rawx));   // g2 = gcd(num(y), den(x))
    myBaseRingValue->myDiv(raw(q1), RefNum(rawx), raw(g1));        // q1 = num(x)/g1
    myBaseRingValue->myDiv(raw(q2), RefNum(rawy), raw(g2));        // q2 = num(y)/q2
    myBaseRingValue->myMul(RefNum(rawlhs), raw(q1), raw(q2));      // num(ans) = q1*q2

    myBaseRingValue->myDiv(raw(q1), RefDen(rawx), raw(g2));        // q1 = den(x)/g2
    myBaseRingValue->myDiv(raw(q2), RefDen(rawy), raw(g1));        // q2 = den(y)/q1
    myBaseRingValue->myMul(RefDen(rawlhs), raw(q1), raw(q2));      // den(ans) = q1*q2
    myBaseRingValue->myNormalizeFracNoGcd(RefNum(rawlhs),RefDen(rawlhs));
    // No need to call CancelFactors here!
  }


  void FractionFieldImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myBaseRingValue->myIsZero(RefNum(rawy)));
    // A "sophisticated" implementation which does "minimal" gcd computation
    RingElem g1(myBaseRingValue), g2(myBaseRingValue), q1(myBaseRingValue), q2(myBaseRingValue), NumAns(myBaseRingValue);
    myBaseRingValue->myGcd(raw(g1), RefNum(rawx), RefNum(rawy));   // g1 = gcd(num(x), num(y))
    myBaseRingValue->myGcd(raw(g2), RefDen(rawx), RefDen(rawy));   // g2 = gcd(den(x), den(y))
    myBaseRingValue->myDiv(raw(q1), RefNum(rawx), raw(g1));        // q1 = num(x)/g1
    myBaseRingValue->myDiv(raw(q2), RefDen(rawy), raw(g2));        // q2 = den(y)/q2
    myBaseRingValue->myMul(raw(NumAns), raw(q1), raw(q2));         // num(ans) = q1*q2

    myBaseRingValue->myDiv(raw(q1), RefDen(rawx), raw(g2));        // q1 = den(x)/g2
    myBaseRingValue->myDiv(raw(q2), RefNum(rawy), raw(g1));        // q2 = den(y)/q1
    myBaseRingValue->myMul(RefDen(rawlhs), raw(q1), raw(q2));      // den(ans) = q1*q2
    myBaseRingValue->mySwap(RefNum(rawlhs), raw(NumAns));          // In case of aliasing between lhs and y
    myBaseRingValue->myNormalizeFracNoGcd(RefNum(rawlhs),RefDen(rawlhs));
    // No need to call CancelFactors here!
  }


  bool FractionFieldImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myBaseRingValue->myIsZero(RefNum(rawy))) return false;
    myDiv(rawlhs, rawx, rawy);
    return true;
  }


  bool FractionFieldImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void FractionFieldImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void FractionFieldImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    // The result is naturally reduced (if the input is).
    myBaseRingValue->myPower(RefNum(rawlhs), RefNum(rawx), n);  // call myPower because RefNum(rawx) could be 1 or -1
    myBaseRingValue->myPower(RefDen(rawlhs), RefDen(rawx), n);  // call myPower because RefDen(rawx) could be 1 or -1
    myBaseRingValue->myNormalizeFracNoGcd(RefNum(rawlhs),RefDen(rawlhs)); // ???ever useful???
  }


  void FractionFieldImpl::mySymbols(std::vector<symbol>& SymList) const
  {
    myBaseRingValue->mySymbols(SymList);
  }


  bool FractionFieldImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    if (!myBaseRingValue->myIsOne(RefDen(rawx))) return false;
    return myBaseRingValue->myIsPrintAtom(RefNum(rawx));
  }


  bool FractionFieldImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    if (myBaseRingValue->myIsPrintedWithMinus(RefNum(rawx)))
    {
      if (myBaseRingValue->myIsMinusOne(RefNum(rawx))) return true;
      RingElem tmp(myBaseRingValue);
      myBaseRingValue->myNegate(raw(tmp), RefNum(rawx));
      if (myBaseRingValue->myIsPrintAtom(raw(tmp))) return true;
    }
    if (!myBaseRingValue->myIsOne(RefDen(rawx))) return false;
    return myBaseRingValue->myIsPrintedWithMinus(RefNum(rawx));
  }


  void FractionFieldImpl::myOutput(std::ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams

    // should we have a function to normalize the denominator?
    // should we modify the rawx in that case?

    if (myBaseRingValue->myIsOne(RefDen(rawx)))
    {
      myBaseRingValue->myOutput(out, RefNum(rawx));
      return;
    }
    if (myBaseRingValue->myIsMinusOne(RefDen(rawx)))
    {
      out << -RingElemAlias(myBaseRing(), RefNum(rawx));
      return;
    }
    // Denom is not 1, so print both numer and denom, perhaps either or both between brackets.
    bool UseBrackets;
    // numerator
    UseBrackets = (!myBaseRingValue->myIsPrintAtom(RefNum(rawx))) &&
      //      (!myBaseRingValue->myIsMinusOne(RefNum(rawx)));
      (!myIsPrintedWithMinus(rawx));
    if (UseBrackets) out << "(";
    myBaseRingValue->myOutput(out, RefNum(rawx));
    if (UseBrackets) out << ")";
    out << "/";
    // denominator
    UseBrackets = !myBaseRingValue->myIsPrintAtom(RefDen(rawx));
    if (UseBrackets) out << "(";
    myBaseRingValue->myOutput(out, RefDen(rawx));
    if (UseBrackets) out << ")";
  }


  void FractionFieldImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    //    out << "FractionField(" << myBaseRingValue << ")";
    out << "RingWithID(" << myID 
        << ", \"FractionField(RingWithID(" << RingID(myBaseRing()) << "))\")";
  }


  void FractionFieldImpl::myOutputSelfLong(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    myOutputSelf(out);
    out << "\n  with BaseRing  ";
    myBaseRing()->myOutputSelfLong(out);
  }


  void FractionFieldImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("setname2", "QuotientField");
    OMOut << myBaseRingValue;
    OMOut->mySendApplyEnd();
  }


  void FractionFieldImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    myBaseRingValue->myOutput_OM(OMOut, RefNum(rawx));
    myBaseRingValue->myOutput_OM(OMOut, RefDen(rawx));
  }


  bool FractionFieldImpl::myIsZero(ConstRawPtr rawx) const
  {
    return myBaseRingValue->myIsZero(RefNum(rawx));
  }


  bool FractionFieldImpl::myIsOne(ConstRawPtr rawx) const
  {
    // NB This definition is valid even if x is not in reduced form.
    return myBaseRingValue->myIsEqual(RefNum(rawx), RefDen(rawx));
  }


  bool FractionFieldImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    RingElem tmp(myBaseRingValue);
    myBaseRingValue->myAdd(raw(tmp), RefNum(rawx), RefDen(rawx));
    return IsZero(tmp);
  }


  bool FractionFieldImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    // Deal with two easy cases first
    if (myBaseRingValue->myIsOne(RefDen(rawx))) return myBaseRingValue->myIsInteger(N, RefNum(rawx));
    if (myBaseRingValue->myIsMinusOne(RefDen(rawx)))
    {
      if (!myBaseRingValue->myIsInteger(N, RefNum(rawx))) return false;
      N = -N;
      return true;
    }
    // General case, must allow for a non-trivial unit in the denominator.
    if (myBaseRingValue->myIsInvertible(RefDen(rawx)))
    {
      RingElem tmp(myBaseRingValue);
      myBaseRingValue->myDiv(raw(tmp), RefNum(rawx), RefDen(rawx));
      return myBaseRingValue->myIsInteger(N, raw(tmp));
    }
    // The lines below work even if FrF elements are not normalized.
    RingElem tmp(myBaseRingValue);
    if (!myBaseRingValue->myIsDivisible(raw(tmp), RefNum(rawx), RefDen(rawx))) return false;
    return myBaseRingValue->myIsInteger(N, raw(tmp));
  }


  // BUG BUG BUG this impl does not work properly if the ring has units other than 1 and -1.
  bool FractionFieldImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { Q = 0; return true; }
    BigInt N,D;
    if (!myBaseRingValue->myIsInteger(D, RefDen(rawx))) return false;
    if (!myBaseRingValue->myIsInteger(N, RefNum(rawx))) return false;
    Q = BigRat(N,D);
    return true;
  }


  // USE DEFAULT myIsZeroAddMul -- see ring.C


  bool FractionFieldImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    // Fractions may not have normalized denominators, so we cannot simply compare nums and dens.
    // Anyway, we speculatively check for equal nums and dens -- it's cheap, and works often???
    if (myBaseRingValue->myIsEqual(RefNum(rawx), RefNum(rawy)) &&
	myBaseRingValue->myIsEqual(RefDen(rawx), RefDen(rawy))) return true;
    RingElem tmp(myBaseRingValue);
    myBaseRingValue->myMul(raw(tmp), RefNum(rawx), RefDen(rawy));
    myBaseRingValue->myNegate(raw(tmp), raw(tmp));
    return myBaseRingValue->myIsZeroAddMul(raw(tmp), RefNum(rawy), RefDen(rawx));
  }


  void FractionFieldImpl::CancelFactors(RawPtr rawx) const
  {
    if (myBaseRingValue->myIsZero(RefNum(rawx))) { myBaseRingValue->myAssign(RefDen(rawx), 1); return; }
    RingElem h(myBaseRingValue);
    myBaseRingValue->myGcd(raw(h), RefNum(rawx), RefDen(rawx));
    myBaseRingValue->myDiv(RefNum(rawx), RefNum(rawx), raw(h));
    myBaseRingValue->myDiv(RefDen(rawx), RefDen(rawx), raw(h));
    /// what if the denominator is negative???  Cannot test this in general!!!
  }


  ideal FractionFieldImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom FractionFieldImpl::myInducedHomCtor(const RingHom& phi) const
  {
    // DO NOT check that the kernel is ideal(0)  -- see documentation!!!
    return RingHom(new InducedHomImpl(FractionField(this), phi));
  }


  RingHom FractionFieldImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    return myInducedHomCtor(phi(theta(myEmbeddingHomCtor())));
  }


  bool FractionFieldImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }


  FractionFieldImpl::InducedHomImpl::InducedHomImpl(const FractionField& FrF, const RingHom& InducingHom):
      RingHomInducedBase(FrF, InducingHom)
  {}


  RingHom FractionFieldImpl::myEmbeddingHomCtor() const
  {
    return RingHom(new EmbeddingHomImpl(FractionField(this)));
  }


  void FractionFieldImpl::InducedHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    RingElem ImDen(myCodomain);
    RingElem ImNum(myCodomain);
    myInducingHom->myApply(raw(ImDen), RefDen(rawarg));
    if (IsZero(ImDen))
      CoCoA_THROW_ERROR1(ERR::BadPartialRingHomArg);
    myInducingHom->myApply(raw(ImNum), RefNum(rawarg));
    if (!myCodomain->myIsDivisible(rawimage, raw(ImNum), raw(ImDen)))
      CoCoA_THROW_ERROR1(ERR::BadPartialRingHomArg);
  }


  FractionFieldImpl::EmbeddingHomImpl::EmbeddingHomImpl(const FractionField& FrF):
      RingHomEmbeddingBase(BaseRing(FrF), FrF)
  {}


  void FractionFieldImpl::EmbeddingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    // myDomain is same as myBaseRing
    myDomain->myAssign(RefNum(rawimage), rawarg);
    myDomain->myAssign(RefDen(rawimage), 1);
  }


  //----------------------------------------------------------------------

  FractionField NewFractionField(const ring& R)
  {
    if (!IsCommutative(R))
      CoCoA_THROW_ERROR1(ERR::NotCommutative);
    if (!IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR1(ERR::NotTrueGCDDomain);
    if (IsZZ(R))  return RingQQ();
    return FractionField(new FractionFieldImpl(R)); // just make a general fraction field for now
  }


  bool IsFractionFieldOfGCDDomain(const ring& R)
  {
    return IsFractionField(R) && IsTrueGCDDomain(BaseRing(R));
  }


  const ring& BaseRing(const FractionField& FrF)
  {
    return FrF->myBaseRing();
  }


  RingHom EmbeddingHom(const FractionField& FrF)
  {
    return FrF->myEmbeddingHomCtor();
  }


  RingHom InducedHom(const FractionField& FrF, const RingHom& InducingHom)
  {
    if (domain(InducingHom) != BaseRing(FrF))
      CoCoA_THROW_ERROR1(ERR::BadInducingHom);
    return FrF->myInducedHomCtor(InducingHom);
  }


  RingElem num(ConstRefRingElem q)
  {
    if (!IsFractionField(owner(q)))
      CoCoA_THROW_ERROR1(ERR::NotElemFrF);
    const FractionField F = owner(q);
    return RingElemAlias(BaseRing(F), F->myRawNum(raw(q)));
  }


  RingElem den(ConstRefRingElem q)
  {
    if (!IsFractionField(owner(q)))
      CoCoA_THROW_ERROR1(ERR::NotElemFrF);
    const FractionField F = owner(q);
    return RingElemAlias(BaseRing(F), F->myRawDen(raw(q)));
  }



} // end of namespace CoCoA
