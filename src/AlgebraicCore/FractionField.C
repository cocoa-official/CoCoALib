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

  // const FractionFieldBase* FractionFieldPtr(const ring& R, const char* const FnName)
  // {
  //   const FractionFieldBase* ptr = FractionFieldPtr(R);
  //   if (ptr == nullptr)
  //     CoCoA_THROW_ERROR(ERR::NotFracField, FnName);
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
    if (!IsIndet(x)) CoCoA_THROW_ERROR(ERR::ReqIndet, "FrF::myDeriv");
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
      CoCoA_THROW_ERROR(ERR::BadPartialRingHomArg, "FractionFieldImpl::InducedHomImpl::myApply");
    myInducingHom->myApply(raw(ImNum), RefNum(rawarg));
    if (!myCodomain->myIsDivisible(rawimage, raw(ImNum), raw(ImDen)))
      CoCoA_THROW_ERROR(ERR::BadPartialRingHomArg, "FractionFieldImpl::InducedHomImpl::myApply");
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
      CoCoA_THROW_ERROR(ERR::NotCommutative, "NewFractionField (pseudo ctor)");
    if (!IsTrueGCDDomain(R))
      CoCoA_THROW_ERROR(ERR::NotTrueGCDDomain, "NewFractionField (pseudo ctor)");
    if (IsZZ(R)) return RingQQ();
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
      CoCoA_THROW_ERROR(ERR::BadInducingHom, "InducedHom(FractionField,ring,RingHom)");
    return FrF->myInducedHomCtor(InducingHom);
  }


  RingElem num(ConstRefRingElem q)
  {
    if (!IsFractionField(owner(q)))
      CoCoA_THROW_ERROR(ERR::NotElemFrF, "num(RingElem)");
    const FractionField F = owner(q);
    return RingElemAlias(BaseRing(F), F->myRawNum(raw(q)));
  }


  RingElem den(ConstRefRingElem q)
  {
    if (!IsFractionField(owner(q)))
      CoCoA_THROW_ERROR(ERR::NotElemFrF, "den(RingElem)");
    const FractionField F = owner(q);
    return RingElemAlias(BaseRing(F), F->myRawDen(raw(q)));
  }



} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/FractionField.C,v 1.62 2024/07/22 15:43:57 abbott Exp $
// $Log: FractionField.C,v $
// Revision 1.62  2024/07/22 15:43:57  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.61  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.60  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.59  2023/03/08 19:43:56  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.58  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.57  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.56  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.55  2021/10/29 14:51:59  abbott
// Summary: Removed keyword virtual (when there is already override)
//
// Revision 1.54  2021/10/29 14:24:28  abbott
// Summary: Inserted keyword override
//
// Revision 1.53  2021/01/07 15:07:03  abbott
// Summary: Corrected copyright
//
// Revision 1.52  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.51  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.50  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.49  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.48  2019/03/04 10:28:33  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.47  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.46  2017/05/09 13:51:30  bigatti
// -- changed BadRingHomArg2 --> BadPartialRingHomArg
//
// Revision 1.45  2015/07/29 13:50:45  bigatti
// -- renamed ID into RingID
//
// Revision 1.44  2015/07/29 11:04:55  bigatti
// -- added space after comma in printing rings
//
// Revision 1.43  2014/07/28 16:04:00  abbott
// Summary: Renamed myEmbeddingHom to myEmbeddingHomCtor
// Author: JAA
//
// Revision 1.42  2014/07/28 15:43:33  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble); Changed (my)EmbeddingHom:
// Author: JAA
//
// Revision 1.41  2014/07/11 15:44:01  bigatti
// -- added  myOutputSelfLong
//
// Revision 1.40  2014/07/11 10:13:39  abbott
// Summary: Added an assertion, removed a pointless IsField in NewFractionField
// Author: JAA
//
// Revision 1.39  2014/07/09 11:36:28  abbott
// Summary: Tidying: moved several simple fns to FractionalFieldBase
// Author: JAA
//
// Revision 1.38  2014/07/08 15:26:28  abbott
// Summary: Updated to fit in with new FractionFieldBase (still needs tidying)
// Author: JAA
//
// Revision 1.37  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.36  2014/07/08 08:43:48  abbott
// Summary: Removed AsFractionField, added FractionFieldPtr
// Author: JAA
//
// Revision 1.35  2014/07/07 12:15:21  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.34  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.33  2014/07/02 16:53:16  bigatti
// -- new way of printing ring with ID
//
// Revision 1.32  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.31  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.30  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.29  2014/01/29 13:05:09  abbott
// Summary: Improved impl of myIsInteger (will work with non-normalized fractions)
// Author: John Abbott
//
// Revision 1.28  2013/05/21 11:32:09  abbott
// Replaced mem fns FractionFieldBase::myGetNum (and Den) by FractionFieldBase::myRawNum (and Den).
//
// Revision 1.27  2013/05/14 14:21:31  abbott
// Revised/improved impl of derivative of ratfns.
//
// Revision 1.26  2013/02/21 14:14:42  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.25  2012/10/24 13:33:34  abbott
// Changed ConstRefRingElem into RingElemAlias in ctor call.
//
// Revision 1.24  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.23  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.22  2012/04/27 15:04:45  abbott
// Added mem fn IamFiniteField
//
// Revision 1.21  2012/02/10 10:26:40  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.20  2012/02/08 17:09:56  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.19  2011/11/09 14:04:30  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.18  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.17  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.16  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.15  2011/05/30 16:16:16  bigatti
// -- fewer parentheses in printing out fractions
// -- more printing tests
//
// Revision 1.14  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.13  2011/02/28 14:11:21  bigatti
// -- no longer prints with denominator (-1)
//
// Revision 1.12  2011/01/28 11:11:38  bigatti
// -- added IsPrintedWithMinus
//
// Revision 1.11  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.10  2010/10/01 15:45:17  bigatti
// -- added mySymbolValue
//
// Revision 1.9  2009/09/24 14:36:46  abbott
// Removed some unnecessary "std::" prefixes.
//
// Revision 1.8  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.7  2009/05/22 10:19:45  bigatti
// -- removed extra space after printing denominator (myOutput)
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2008/04/22 13:09:16  abbott
// Removed IsPositive and IsNegative functions for ZZs.
// Comparison between RingElem and 0 is now handled specially (specially fast?).
//
// Revision 1.4  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.2  2007/05/21 12:57:28  abbott
// New class for passing machine integers as args; includes some simple
// operations on machine integers (cmp, gcd, IsNegative,...).  Operations
// between ZZ and machine integers modified to use the new class.  Inexact
// integer division (of a ZZ) by a negative value now triggers an error;
// new error for reporting inexact integer division by a negative value.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.19  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.18  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.17  2007/03/03 15:59:05  cocoa
// -- added  #include "CoCoA/RingQ.H"
//
// Revision 1.16  2007/03/03 14:07:23  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.15  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.14  2007/01/15 15:47:57  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.13  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.12  2006/12/07 11:57:53  cocoa
// -- style: RawPtr args are now called "raw.."
//
// Revision 1.11  2006/11/27 14:21:13  cocoa
// -- reorganised #include files
//
// Revision 1.10  2006/11/09 17:36:46  cocoa
// -- fixed: myEmbeddingHom::myApply
//
// Revision 1.9  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.8  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.7  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.6  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.5  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/10/06 10:15:52  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.7  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.6  2006/04/21 15:01:36  cocoa
// Changed default implementation of RingBase::myGcd -- it now gives a SERIOUS
// error.  All fields must now handle a call to gcd explicitly: they can use
// the new myGcdInField function.  It's now cleaner than it was.
//
// Revision 1.5  2006/04/05 13:36:32  cocoa
// -- fixed: myOutput
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.5  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
//
// Revision 1.4  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.3  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.6  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/03/31 16:59:42  cocoa
// Made special matrix ctors private, so a user has to pass via the
// pseudo-ctors (which do all the arg sanity checking).
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.15  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.14  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.13  2004/11/05 15:34:33  cocoa
// Consequential change following from the renaming of
// FieldIdealImpl and the introduction of the new pseudo-ctor.
//
// Revision 1.12  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.11  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.10  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
// Revision 1.9  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.8  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.7  2004/03/20 17:46:11  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.6  2004/02/03 16:16:21  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.5  2004/01/30 14:07:10  cocoa
// Tidied RingRawValue union: now it contains just two fields,
// and has no need of forward declarations of types used internally
// by the concrete rings -- it uses explicitly a void* instead.
//
// I have tidied the "import" functions used by most concrete rings.
//
// I have moved the choice of representation type for RingFp and RingFpLog
// into a typedef in config.H -- this is to recognise that different
// choices may work best on different platforms.
//
// Revision 1.4  2003/11/14 13:06:05  cocoa
// -- New function "myIsPrintAtom" for printing polynomials and fractions
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:16:39  cocoa
// New coding convention for rings.
//
// Revision 1.8  2003/06/23 16:46:11  abbott
// Minor cleaning prior to public release.
// Mostly consequential name changes.
//
// Revision 1.7  2003/05/27 15:51:52  abbott
// Changed name from FractionField to GeneralFractionField.
// Added EmbeddingHom and EmbeddingHom_t (to map from R into FrF(R)).
// Improved the printing function.
//
// Revision 1.6  2003/05/15 14:46:33  abbott
// Corrected two embarrassing mistakes in add and sub.
//
// Revision 1.5  2003/05/15 10:21:05  abbott
// Consequential changes from ring.H: old changes to names and some
// member functions.  Added homomorphisms and ideals (reqd by ring.H).
// Uses FieldIdeal for ideals.  Added implementation of homomorphism.
//
// Revision 1.4  2002/11/15 10:42:28  abbott
// Revised according to the renaming in ring.H.
// A little tidying up, and better error messages.
//
// Revision 1.3  2002/07/05 15:31:59  abbott
// Modified dtor since we now keep myZero in an auto_ptr.
//
// Revision 1.2  2002/07/05 15:13:10  abbott
// Added IsDivisible member function, and "improved" the two IsEqual member functions.
//
// Revision 1.1  2002/06/22 16:47:05  abbott
// Initial revision
//
//
