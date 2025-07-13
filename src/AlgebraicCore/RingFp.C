//   Copyright (c)  2001-2011,2014,2021  John Abbott and Anna M. Bigatti

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


#include "CoCoA/RingFp.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"

#include <algorithm>
//using std::swap;         // only in mySwap
#include <iostream>
using std::ostream;        // only in myOutput
// #include <limits>  ---  included in MachineInt.H (included via BigRat.H)
using std::numeric_limits; // only in ctor
// #include <memory>  ---  included in MemPool.H
using std::unique_ptr;
// #include <vector>  ---  included in ideal.H
using std::vector;


namespace CoCoA
{

  class RingFpImpl: public QuotientRingBase
  {
  private: // data members
    typedef SmallFpImpl::value value_t;
    const long myModulus;
    const SmallFpImpl myImpl;
    mutable MemPool myMemMgr;       // MemPool must come *BEFORE* myZeroPtr and myOnePtr
    unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    unique_ptr<RingElem> myFieldGenPtr;   ///< Finite field has primitive elem

  private: // auxiliary functions
    static long PrincipalGen(const ideal& I); // used for arg checking in ctor
    static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(RingElemConstRawPtr rawx);

  private: // ctor and pseudo-ctors
    explicit RingFpImpl(const ideal& I, GlobalSettings::ResidueRepr repr = DefaultResidueRepr()); // called only by NewRingFp
    ~RingFpImpl();
    friend QuotientRing NewRingFp(const MachineInt& p, GlobalSettings::ResidueRepr repr);
    friend QuotientRing NewRingFp(const BigInt& P);
    friend QuotientRing NewRingFp(const ideal& I);
  public: // disable copy ctor & assignment
    RingFpImpl(const RingFpImpl&) = delete;
    RingFpImpl& operator=(const RingFpImpl&) = delete;

  public: // functions every ring must implement
    BigInt myCharacteristic() const override;
    long myLogCardinality() const override;
    bool IamCommutative() const override;
    bool3 IamIntegralDomain3(bool) const override;
    bool IamField() const override;
    bool IamFiniteField() const override;
    bool IamExact() const override;

    RingElem myFieldGen() const override;
    ConstRefRingElem myZero() const override;
    ConstRefRingElem myOne() const override;
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& N) const override;
    RingElemRawPtr myNew(ConstRawPtr rawt) const override;
    void myDelete(RawPtr rawx) const override;                                      // destroys x (incl all resources)

    void mySwap(RawPtr rawx, RawPtr rawy) const override;                           // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;                  // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;               // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;                   // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                                // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;

    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;                  // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                           // true iff x is invertible
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;   // lhs = gcd(x,y) in a field

    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;              // out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;
    void myOutputSelfShort(std::ostream& out) const override;                       // out << R
    void myOutputSelf(std::ostream& out) const override;                            // out << R
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                     // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;       // OMOut << x

    bool myIsZero(ConstRawPtr rawx) const override;                                 // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                  // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                             // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                   // always true
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                  // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                    // false iff x overflows
    bool myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const override;// lhs += y*z, result says whether lhs == 0.
    bool myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawy, ConstRawPtr rawz) const override;// lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;              // x == y

    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))

    bool myImageLiesInSubfield(const RingHom& phi) const override;

  protected:
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;        // lhs = x^n, n>1, x not -1,0,1
    void myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const override; // lhs = x^N, N big, x not -1,0,1

  public: // functions every QuotientRing must implement
    RingElem myCanonicalRepr(ConstRawPtr rawx) const override; // result is element of myReprRing
    void myReduction(RawPtr rawimage, ConstRawPtr rawarg) const override;
    RingHom myInducedHomCtor(const RingHom& InducingHom) const override;

    // Special fns for RingFp only
    const SmallFpImpl& myModularArith() const;
  private: // Impl detail: hom class
    class InducedHomImpl: public RingHomBase
    {
    public:
      InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom);
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return false; }
    };

  };



  inline RingFpImpl::value_t& RingFpImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const RingFpImpl::value_t& RingFpImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }


  // Returns generator of I as a long; returns 0 if value is too large to fit.
  long RingFpImpl::PrincipalGen(const ideal& I)
  {
    if (IsZero(I)) return 0;
    const BigInt GenI = ConvertTo<BigInt>(TidyGens(I)[0]);
    long p;
    if (!IsConvertible(p, GenI))  // check that the value of the principal generator will fit
      return 0; // 0 will trigger an error in ctor for SmallFpImpl
    return p;
  }


  RingFpImpl::RingFpImpl(const ideal& I, GlobalSettings::ResidueRepr repr):
      QuotientRingBase(RingZZ(), I),  // confirms that I is an ideal of ZZ
      myModulus(PrincipalGen(I)),
      myImpl(myModulus, repr),  // also checks that myModulus is a small prime
      myMemMgr(SmallFpImpl::ourDatumSize, "RingFpImpl.myMemMgr")
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myFieldGenPtr.reset(new RingElem(ring(this), PrimitiveRoot(myImpl.myModulus())));
    myRefCountZero();
  }


  RingFpImpl::~RingFpImpl()
  {}


  BigInt RingFpImpl::myCharacteristic() const
  {
    return BigInt(myModulus);
  }


  long RingFpImpl::myLogCardinality() const
  {
    return 1;
  }


  bool RingFpImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingFpImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingFpImpl::IamField() const
  {
    return true;
  }


  bool RingFpImpl::IamFiniteField() const
  {
    return true;
  }


  bool RingFpImpl::IamExact() const
  {
    return true;
  }


  RingElem RingFpImpl::myFieldGen() const
  {
    return *myFieldGenPtr;
  }


  ConstRefRingElem RingFpImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingFpImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingFpImpl::myNew() const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    new(ans) value_t(); // init to 0
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpImpl::myNew(const MachineInt& n) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    new(ans) value_t();
    *ans = myImpl.myReduce(n);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpImpl::myNew(const BigInt& N) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    new(ans) value_t();
    *ans = myImpl.myReduce(N);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpImpl::myNew(ConstRawPtr rawy) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    new(ans) value_t(import(rawy));
    return RingElemRawPtr(ans);
  }


  void RingFpImpl::myDelete(RawPtr rawx) const
  {
    import(rawx).~value(); // dtor actually does nothing (currently)
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingFpImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    std::swap(import(rawx), import(rawy));
  }


  void RingFpImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingFpImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = myImpl.myReduce(n);
  }


  void RingFpImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = myImpl.myReduce(N);
  }


  void RingFpImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs) = zero(SmallFp);
  }


  void RingFpImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
  }


  void RingFpImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = myImpl.myNegate(import(rawx));
  }


  void RingFpImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myAdd(import(rawx), import(rawy));
  }


  void RingFpImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.mySub(import(rawx), import(rawy));
  }


  void RingFpImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myMul(import(rawx), import(rawy));
  }


  void RingFpImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    import(rawlhs) = myImpl.myDiv(import(rawx), import(rawy));
  }


  bool RingFpImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (IsZero(import(rawy))) return false;
    import(rawlhs) = myImpl.myDiv(import(rawx), import(rawy));
    return true;
  }


  bool RingFpImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingFpImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void RingFpImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    import(rawlhs) = myImpl.myPower(import(rawx), n);
  }

  void RingFpImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(N > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    // Use Fermat's Little Theorem to reduce exponent...
    import(rawlhs) = myImpl.myPower(import(rawx), N%(myModulus-1));
  }


  void RingFpImpl::myOutput(ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));
    out << myImpl.myExport(import(rawx));
  }


  bool RingFpImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return myImpl.myExport(import(rawx)) >= 0;
  }


  bool RingFpImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return myImpl.myExport(import(rawx)) < 0;
  }


  void RingFpImpl::myOutputSelfShort(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "ZZ/(" << myModulus << ")";
  }

  void RingFpImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "RingWithID(" << myID << ", \"ZZ/(" << myModulus << ")\")";
  }


  void RingFpImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("setname2", "GFp");
    OMOut << myModulus;
    OMOut->mySendApplyEnd();
  }


  void RingFpImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut << myImpl.myExport(import(rawx));
  }


  bool RingFpImpl::myIsZero(ConstRawPtr rawx) const
  {
    return IsZero(import(rawx));
  }


  bool RingFpImpl::myIsOne(ConstRawPtr rawx) const
  {
    return IsOne(import(rawx));
  }


  bool RingFpImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return import(rawx) == myImpl.myReduce(-1);
  }


  bool RingFpImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    N = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    Q = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    d = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpImpl::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  {
    return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  }


  bool RingFpImpl::myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  { // same as above: just to avoid calling RingBase::myIsZeroAddMul with 4 args
    return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  }


  bool RingFpImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return import(rawx) == import(rawy);
  }




  ideal RingFpImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom RingFpImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    // No need to check compatibility -- it was checked when theta and phi were built
    return RingHom(new InducedHomImpl(QuotientRing(this), phi(theta(myQuotientingHomCtor()))));
  }


  bool RingFpImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }



  RingElem RingFpImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    return RingElem(myReprRing, myImpl.myExport(import(rawx)));
  }


  void RingFpImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    BigInt tmp;
    CoCoA_ASSERT(myReprRing->myIsInteger(tmp, rawarg));
    myReprRing->myIsInteger(tmp, rawarg);
    import(rawimage) = myImpl.myReduce(tmp);
  }


  RingHom RingFpImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    // Compatibility has already been checked (see InducedHom in QuotientRing.C)
    CoCoA_ASSERT(IsZero(InducingHom(myModulus)));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  const SmallFpImpl& RingFpImpl::myModularArith() const
  {
    return myImpl;
  }


  //---------------------------------------------------------------------------
  // Functions to do with ring homomorphisms


  RingFpImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom):
      RingHomBase(domain, codomain(InducingHom))
  { /* Compatibility already checked in InducedHom in QuotientRing.C */  }


  void RingFpImpl::InducedHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    BigInt tmp;  //??? wasteful new/delete
    CoCoA_ASSERT(myDomain->myIsInteger(tmp, rawarg));
    myDomain->myIsInteger(tmp, rawarg);  // must necessarily succeed
    myCodomain->myAssign(rawimage, tmp);
  }



  QuotientRing NewRingFp(const MachineInt& p, GlobalSettings::ResidueRepr repr)
  {
    return QuotientRing(new RingFpImpl(ideal(RingElem(RingZZ(), p)), repr));
  }

  QuotientRing NewRingFp(const BigInt& P)
  {
    return QuotientRing(new RingFpImpl(ideal(RingElem(RingZZ(), P))));
  }

  QuotientRing NewRingFp(const ideal& I)
  {
    if (!IsZZ(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::IdealNotInRing);
    return QuotientRing(new RingFpImpl(I));
  }


  bool IsGoodForRingFp(const MachineInt& p)
  {
    if (IsNegative(p) || !IsSignedLong(p))  return false;
    const long n = AsSignedLong(p);
    return SmallFpImpl::IsGoodCtorArg(n);
  }

  bool IsGoodForRingFp(const BigInt& P)
  {
    if (P <= 0)  return false;
    long p;
    if (!IsConvertible(p, P))  return false;
    return IsGoodForRingFp(p);
  }

  bool IsGoodForRingFp(const ideal& I)
  {
    if (!IsZZ(RingOf(I)))  return false;
    if (IsZero(I))  return false;
    return IsGoodForRingFp(ConvertTo<BigInt>(TidyGens(I)[0]));
  }


  bool IsRingFp(const ring& R)
  {
    return dynamic_cast<const RingFpImpl*>(R.myRawPtr()) != nullptr;
  }

  const SmallFpImpl& ModularArith(const ring& R)
  {
    return dynamic_cast<const RingFpImpl*>(R.myRawPtr())->myModularArith();
  }

} // end of namespace CoCoA
