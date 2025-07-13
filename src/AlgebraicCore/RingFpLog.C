//   Copyright (c)  2002-2011,2014,2021  John Abbott and Anna M. Bigatti

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

#include "CoCoA/RingFpLog.H"

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
#include "CoCoA/SmallFpLogImpl.H"
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

  class RingFpLogImpl: public QuotientRingBase
  {
  private: // data members
    typedef SmallFpLogImpl::value_t value_t;
    const value_t myModulus;
    const SmallFpLogImpl myImpl;
    mutable MemPool myMemMgr;       // MemPool must come *BEFORE* myZeroPtr and myOnePtr
    unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    unique_ptr<RingElem> myFieldGenPtr;   ///< Finite field has primitive elem

  private: // auxiliary functions
    static value_t PrincipalGen(const ideal& I); // used for arg checking in ctor
    static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(RingElemConstRawPtr rawx);

  private:
    RingFpLogImpl(const ideal& I, GlobalSettings::ResidueRepr repr = DefaultResidueRepr()); // called only by NewRingFpLog
    ~RingFpLogImpl();
    friend QuotientRing NewRingFpLog(const MachineInt& p, GlobalSettings::ResidueRepr repr);
    friend QuotientRing NewRingFpLog(const BigInt& P);
    friend QuotientRing NewRingFpLog(const ideal& I);
  public: // disable copy ctor & assignment
    RingFpLogImpl(const RingFpLogImpl&) = delete;
    RingFpLogImpl& operator=(const RingFpLogImpl&) = delete;
  public:

    // functions every ring must implement
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
    void myDelete(RawPtr rawx) const override;                                           // destroys x (incl all resources)

    void mySwap(RawPtr rawx, RawPtr rawy) const override;                                // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;                       // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;                    // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;                        // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                                     // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;

    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;                       // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                                // true iff x is invertible
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = gcd(x,y) in a field

    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;                   // out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;
    void myOutputSelf(std::ostream& out) const override;                                 // out << R
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                             // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;               // OMOut << x

    bool myIsZero(ConstRawPtr rawx) const override;                                      // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                       // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                                  // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                        // always true
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                       // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                         // false iff x overflows
    bool myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const override;// lhs += y*z, result says whether lhs == 0.
    bool myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawy, ConstRawPtr rawz) const override;// lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;                    // x == y

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


  private:
    class InducedHomImpl: public RingHomBase
    {
    public:
      InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom);
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return false; }
    };

  };



  inline RingFpLogImpl::value_t& RingFpLogImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const RingFpLogImpl::value_t& RingFpLogImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }


  // Returns generator of I as a value_t; returns 0 if value is too large to fit.
  RingFpLogImpl::value_t RingFpLogImpl::PrincipalGen(const ideal& I)
  {
    if (IsZero(I)) return 0;
    const BigInt GenI = ConvertTo<BigInt>(TidyGens(I)[0]);
    value_t p;
    if (!IsConvertible(p, GenI))  // check that the value of the principal generator will fit
      return 0;
    return p;
  }


  RingFpLogImpl::RingFpLogImpl(const ideal& I, GlobalSettings::ResidueRepr repr):
      QuotientRingBase(RingZZ(), I),  // confirms that I is an ideal of Z
      myModulus(PrincipalGen(I)),
      myImpl(myModulus, repr),   // also checks that myModulus is a small prime
      myMemMgr(SmallFpLogImpl::ourDatumSize, "RingFpLogImpl.myMemMgr")
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myFieldGenPtr.reset(new RingElem(ring(this), PrimitiveRoot(myImpl.myModulus())));
    myRefCountZero();
  }


  RingFpLogImpl::~RingFpLogImpl()
  {}


  BigInt RingFpLogImpl::myCharacteristic() const
  {
    return BigInt(myModulus);
  }


  long RingFpLogImpl::myLogCardinality() const
  {
    return 1;
  }


  bool RingFpLogImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingFpLogImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingFpLogImpl::IamField() const
  {
    return true;
  }


  bool RingFpLogImpl::IamFiniteField() const
  {
    return true;
  }


  bool RingFpLogImpl::IamExact() const
  {
    return true;
  }


  RingElem RingFpLogImpl::myFieldGen() const
  {
    return *myFieldGenPtr;
  }


  ConstRefRingElem RingFpLogImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingFpLogImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingFpLogImpl::myNew() const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = 0;
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpLogImpl::myNew(const MachineInt& n) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myReduce(n);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpLogImpl::myNew(const BigInt& N) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myReduce(N);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFpLogImpl::myNew(ConstRawPtr rawy) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = import(rawy);
    return RingElemRawPtr(ans);
  }


  void RingFpLogImpl::myDelete(RawPtr rawx) const
  {
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingFpLogImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    std::swap(import(rawx), import(rawy));
  }


  void RingFpLogImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingFpLogImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = myImpl.myReduce(n);
  }


  void RingFpLogImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = myImpl.myReduce(N);
  }


  void RingFpLogImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs) = 0;
  }


  void RingFpLogImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
  }


  void RingFpLogImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = myImpl.myNegate(import(rawx));
  }


  void RingFpLogImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myAdd(import(rawx), import(rawy));
  }


  void RingFpLogImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.mySub(import(rawx), import(rawy));
  }


  void RingFpLogImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myMul(import(rawx), import(rawy));
  }


  void RingFpLogImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    import(rawlhs) = myImpl.myDiv(import(rawx), import(rawy));
  }


  bool RingFpLogImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) return false;
    myDiv(rawlhs, rawx, rawy);
    return true;
  }


  bool RingFpLogImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingFpLogImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void RingFpLogImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    import(rawlhs) = myImpl.myPower(import(rawx), n);
  }

  void RingFpLogImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(N > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    // Use Fermat's Little Theorem to reduce exponent...
    import(rawlhs) = myImpl.myPower(import(rawx), N%(myModulus-1));
  }


  void RingFpLogImpl::myOutput(ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));
    out << myImpl.myExport(import(rawx));
  }


  bool RingFpLogImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return myImpl.myExport(import(rawx)) >= 0;
  }


  bool RingFpLogImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return myImpl.myExport(import(rawx)) < 0;
  }


  void RingFpLogImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    //    out << "FFp(" << myModulus << ")";
    out << "RingWithID(" << myID << ",\"ZZ/(" << myModulus << ")\")";
  }


  void RingFpLogImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("setname2", "GFp");
    OMOut << myModulus;
    OMOut->mySendApplyEnd();
  }


  void RingFpLogImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut << myImpl.myExport(import(rawx));
  }


  bool RingFpLogImpl::myIsZero(ConstRawPtr rawx) const
  {
    return (import(rawx) == 0);
  }


  bool RingFpLogImpl::myIsOne(ConstRawPtr rawx) const
  {
    return (import(rawx) == 1);
  }


  bool RingFpLogImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return (import(rawx) == myImpl.myReduce(-1));
  }


  bool RingFpLogImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    N = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpLogImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    Q = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpLogImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    d = myImpl.myExport(import(rawx));
    return true;
  }


  bool RingFpLogImpl::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  {
    return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  }


  bool RingFpLogImpl::myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  { // same as above: just to avoid calling RingBase::myIsZeroAddMul with 4 args
    return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  }


  bool RingFpLogImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return (import(rawx) == import(rawy));
  }




  ideal RingFpLogImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom RingFpLogImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    // No need to check compatibility -- it was checked when theta and phi were built
    return RingHom(new InducedHomImpl(QuotientRing(this), phi(theta(myQuotientingHomCtor()))));
  }


  bool RingFpLogImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }


  RingElem RingFpLogImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    return RingElem(myReprRing, myImpl.myExport(import(rawx)));
  }


  void RingFpLogImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    BigInt tmp;
    CoCoA_ASSERT(myReprRing->myIsInteger(tmp, rawarg));
    myReprRing->myIsInteger(tmp, rawarg);
    import(rawimage) = myImpl.myReduce(tmp);
  }


  RingHom RingFpLogImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    // Compatibility has already been checked (see InducedHom in QuotientRing.C)
    CoCoA_ASSERT(IsZero(InducingHom(myModulus)));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  //---------------------------------------------------------------------------
  // Functions to do with ring homomorphisms


  RingFpLogImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom):
      RingHomBase(domain, codomain(InducingHom))
  { /* Compatibility already checked in InducedHom in QuotientRing.C */  }


  void RingFpLogImpl::InducedHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    BigInt tmp;  //??? wasteful new/delete
    CoCoA_ASSERT(myDomain->myIsInteger(tmp, rawarg));
    myDomain->myIsInteger(tmp, rawarg);  // must necessarily succeed
    myCodomain->myAssign(rawimage, tmp);
  }



  QuotientRing NewRingFpLog(const MachineInt& p, GlobalSettings::ResidueRepr repr)
  {
    return QuotientRing(new RingFpLogImpl(ideal(RingElem(RingZZ(), p)), repr));
  }

  QuotientRing NewRingFpLog(const BigInt& P)
  {
    return QuotientRing(new RingFpLogImpl(ideal(RingElem(RingZZ(), P))));
  }

  QuotientRing NewRingFpLog(const ideal& I)
  {
    if (!IsZZ(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::IdealNotInRing);
    return QuotientRing(new RingFpLogImpl(I));
  }


  bool IsGoodForRingFpLog(const MachineInt& p)
  {
    if (IsNegative(p) || !IsSignedLong(p))  return false;
    const long n = AsSignedLong(p);
    return SmallFpLogImpl::IsGoodCtorArg(n);
  }

  bool IsGoodForRingFpLog(const BigInt& P)
  {
    if (P <= 0)  return false;
    long p;
    if (!IsConvertible(p, P))  return false;
    return IsGoodForRingFpLog(p);
  }

  bool IsGoodForRingFpLog(const ideal& I)
  {
    if (!IsZZ(RingOf(I)))  return false;
    if (IsZero(I))  return false;
    return IsGoodForRingFpLog(ConvertTo<BigInt>(TidyGens(I)[0]));
  }


} // end of namespace CoCoA
