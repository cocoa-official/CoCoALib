//   Copyright (c)  2015,2021  John Abbott and Anna M. Bigatti

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

#include "CoCoA/RingFqVec.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SmallFqUtils.H"
#include "CoCoA/SmallFqVecImpl.H"
#include "CoCoA/SparsePolyIter.H" // for myReduction
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
#include "CoCoA/error.H"

#include <iostream>
using std::ostream;

//#include <vector>
using std::vector;

#include <memory>
using std::unique_ptr;

namespace CoCoA
{
  class RingFqVecImpl: public QuotientRingBase
  {
  private: // data members
    typedef SmallFpImpl::value FpElem;
    const long myDeg;
    const FFqImpl_vec myImpl;
//    const FpElem myModulus;
//    const SmallFpImpl myImpl;
    mutable MemPool myMemMgr;       // MemPool must come *BEFORE* myZeroPtr and myOnePtr
    unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    unique_ptr<RingElem> myFieldGenPtr;   ///< Every finite field stores its own cyclic generator.

  private: // auxiliary functions
    static FpElem PrincipalGen(const ideal& I); // used for arg checking in ctor
    static FpElem* import(RingElemRawPtr rawx);
    static const FpElem* import(RingElemConstRawPtr rawx);

  private:
    explicit RingFqVecImpl(long p, int d); // called only by NewRingFqVec
    RingFqVecImpl(const ideal& I);
    ~RingFqVecImpl();
    friend ring NewRingFqVec(const MachineInt& p, const MachineInt& d); // create field of size p^d
    friend ring NewRingFqVec(const BigInt& p, const MachineInt& d); // create field of size p^d
    friend ring NewRingFqVec(const MachineInt& p, const BigInt& d); // create field of size p^d
    friend ring NewRingFqVec(const BigInt& p, const BigInt& d); // create field of size p^d
  public: // disable copy ctor & assignment
    RingFqVecImpl(const RingFqVecImpl&) = delete;
    RingFqVecImpl& operator=(const RingFqVecImpl&) = delete;
  public:

    // functions every ring must implement
    BigInt myCharacteristic() const override;
    long myLogCardinality() const override;
    bool IamCommutative() const override;
    bool3 IamIntegralDomain3(bool) const override;
    bool IamField() const override;
    bool IamFiniteField() const override;
    bool IamExact() const override;

    RingElem myFieldGen() const override { return *myFieldGenPtr; }
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
    void myOutputSelf(std::ostream& out) const override;                            // out << R
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                        // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;          // OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                                 // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                  // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                             // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                   // always true
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                  // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                    // false iff x overflows
    // bool myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const override;// lhs += y*z, result says whether lhs == 0.
    // bool myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawy, ConstRawPtr rawz) const override;// lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;              // x == y
    RingElem mySymbolValue(const symbol& /*sym*/) const override;
    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))

    bool myImageLiesInSubfield(const RingHom& phi) const override;

  protected:
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;   // lhs = x^n, n>1, x not -1,0,1
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
    private:
      RingHom myInducingHom;
    };

  };


  inline RingFqVecImpl::FpElem* RingFqVecImpl::import(RingElemRawPtr rawx)
  {
    return static_cast<FpElem*>(rawx.myRawPtr());
  }

  inline const RingFqVecImpl::FpElem* RingFqVecImpl::import(RingElemConstRawPtr rawx)
  {
    return static_cast<const FpElem*>(rawx.myRawPtr());
  }


  // // Returns generator of I as a FpElem; returns 0 if value is too large to fit.
  // RingFqVecImpl::FpElem RingFqVecImpl::PrincipalGen(const ideal& I)
  // {
  //   if (IsZero(I)) return 0;
  //   const BigInt GenI = ConvertTo<BigInt>(TidyGens(I)[0]);
  //   FpElem p;
  //   if (!IsConvertible(p, GenI))  // check that the value of the principal generator will fit
  //     return 0;
  //   return p;
  // }



  RingFqVecImpl::RingFqVecImpl(const ideal& I):
      QuotientRingBase(RingOf(I), I),
      myDeg(deg(gens(I)[0])),
      myImpl(I),
      myMemMgr(myDeg*SmallFpImpl::ourDatumSize, "RingFqVecImpl.myMemMgr"),
      myZeroPtr(),
      myOnePtr(),
      myFieldGenPtr()
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myFieldGenPtr.reset(new RingElem(ring(this)));
    myImpl.myGen(import(raw(*myFieldGenPtr)));
////    import(raw(*myFieldGenPtr)) = myImpl.myGen();
    myRefCountZero();
  }


  RingFqVecImpl::~RingFqVecImpl()
  {}


  BigInt RingFqVecImpl::myCharacteristic() const
  {
    return BigInt(myImpl.myModulus());
  }


  long RingFqVecImpl::myLogCardinality() const
  {
    return myImpl.myExtnDeg();
  }


  bool RingFqVecImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingFqVecImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingFqVecImpl::IamField() const
  {
    return true;
  }


  bool RingFqVecImpl::IamFiniteField() const
  {
    return true;
  }


  bool RingFqVecImpl::IamExact() const
  {
    return true;
  }


  ConstRefRingElem RingFqVecImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingFqVecImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingFqVecImpl::myNew() const
  {
    FpElem* ans = static_cast<FpElem*>(myMemMgr.alloc());
    for (int i=0;i<myDeg;++i) ans[i]=zero(SmallFp);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFqVecImpl::myNew(const MachineInt& n) const
  {
    FpElem* ans = static_cast<FpElem*>(myMemMgr.alloc());
    ans[0] = myImpl.myFpArith().myReduce(n);
    for (int i=1;i<myDeg;++i) ans[i]=zero(SmallFp);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFqVecImpl::myNew(const BigInt& N) const
  {
    FpElem* ans = static_cast<FpElem*>(myMemMgr.alloc());
    ans[0] = myImpl.myFpArith().myReduce(N);
    for (int i=1;i<myDeg;++i) ans[i]=zero(SmallFp);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFqVecImpl::myNew(ConstRawPtr rawy) const
  {
    FpElem* ans = static_cast<FpElem*>(myMemMgr.alloc());
    for (int i=0;i<myDeg;++i) ans[i]=import(rawy)[i];//    *ans = import(rawy);
    return RingElemRawPtr(ans);
  }


  void RingFqVecImpl::myDelete(RawPtr rawx) const
  {
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingFqVecImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    for (int i=0;i<myDeg;++i)
      std::swap(import(rawx)[i], import(rawy)[i]);
  }


  void RingFqVecImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    for (int i=0;i<myDeg;++i)
    import(rawlhs)[i] = import(rawx)[i];
  }


  void RingFqVecImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs)[0] = myImpl.myFpArith().myReduce(n);
    for (int i=1;i<myDeg;++i)
      import(rawlhs)[i] = zero(SmallFp);
  }


  void RingFqVecImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs)[0] = myImpl.myFpArith().myReduce(N);
    for (int i=1;i<myDeg;++i)
      import(rawlhs)[i] = zero(SmallFp);
  }


  void RingFqVecImpl::myAssignZero(RawPtr rawlhs) const
  {
    for (int i=0;i<myDeg;++i)
      import(rawlhs)[i] = zero(SmallFp);
  }


  void RingFqVecImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingFqVecImpl::myRecvTwinFloat");
  }


  void RingFqVecImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    for (int i=0;i<myDeg;++i)
      import(rawlhs)[i] = myImpl.myFpArith().myNegate(import(rawx)[i]);
  }


  void RingFqVecImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myImpl.myAdd(import(rawlhs), import(rawx), import(rawy));
  }


  void RingFqVecImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myImpl.mySub(import(rawlhs), import(rawx), import(rawy));
  }


  void RingFqVecImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myImpl.myMul(import(rawlhs), import(rawx), import(rawy));
  }


  void RingFqVecImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    myImpl.myDiv(import(rawlhs), import(rawx), import(rawy));
  }


  bool RingFqVecImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (import(rawy) == 0) return false;
    myImpl.myDiv(import(rawlhs), import(rawx), import(rawy));
    return true;
  }


  bool RingFqVecImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingFqVecImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void RingFqVecImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    myImpl.myPower(import(rawlhs), import(rawx), n);
  }

  void RingFqVecImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(N > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    // Use Fermat's Little Theorem to reduce exponent...
    const long card = SmallPower(myImpl.myModulus(), myImpl.myExtnDeg());
    myImpl.myPower(import(rawlhs), import(rawx), N%(card-1));
  }


  void RingFqVecImpl::myOutput(ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));

    const SmallFpImpl& ModP(myImpl.myFpArith());
    const SmallFpImpl::value* v = import(rawx);
    out << "[" << ModP.myExport(v[0]);
    for (int i=1; i < myDeg; ++i)
      out << ", " << ModP.myExport(v[i]);
    out <<"]";
  }


  bool RingFqVecImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return false; ///??? BUG BUG BUG ???probably OK if image of an int?
//    return myImpl.myExport(import(rawx)) >= 0;
  }


  bool RingFqVecImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return false; ///??? BUG BUG BUG   STOPGAP impl, see comment in isprintatom
//    return myImpl.myExport(import(rawx)) < 0;
  }


  void RingFqVecImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    //    out << "FFp(" << myModulus << ")";
    out << "RingWithID(" << myID << ",\"FFqVec(p=" << myImpl.myModulus() << ", deg=" << myDeg << ")\")";
  }


  void RingFqVecImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    // OMOut->mySendApplyStart();
    // OMOut << OpenMathSymbol("setname2", "GFp");
    // OMOut << myModulus;
    // OMOut->mySendApplyEnd();
  }


  void RingFqVecImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingFqVecImpl::myOutput");
    OMOut << 0; //???myImpl.myExport(import(rawx));
  }


  bool RingFqVecImpl::myIsZero(ConstRawPtr rawx) const
  {
    for (int i=0;i<myDeg;++i)
      if (!IsZero(import(rawx)[i])) return false;
    return true;
  }


  bool RingFqVecImpl::myIsOne(ConstRawPtr rawx) const
  {
    for (int i=1;i<myDeg;++i)
      if (!IsZero(import(rawx)[i])) return false;
    return IsOne(import(rawx)[0]);
  }


  bool RingFqVecImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    for (int i=1;i<myDeg;++i)
      if (!IsZero(import(rawx)[i])) return false;
    return import(rawx)[0] == myImpl.myFpArith().myReduce(-1);
  }


  bool RingFqVecImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    for (int i=1;i<myDeg;++i)
      if (!IsZero(import(rawx)[i])) return false;
    N = myImpl.myFpArith().myExport(import(rawx)[0]);
    return true;
  }


  bool RingFqVecImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    for (int i=1;i<myDeg;++i)
      if (!IsZero(import(rawx)[i])) return false;
    Q = myImpl.myFpArith().myExport(import(rawx)[0]);
    return true;
  }


  bool RingFqVecImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    for (int i=1;i<myDeg;++i)
      if (!IsZero(import(rawx)[i])) return false;
    d = myImpl.myFpArith().myExport(import(rawx)[0]);
    return true;
  }


  // bool RingFqVecImpl::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  // {
  //   return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  // }


  // bool RingFqVecImpl::myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  // { // same as above: just to avoid calling RingBase::myIsZeroAddMul with 4 args
  //   return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  // }


  bool RingFqVecImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    for (int i=0;i<myDeg;++i)
      if (import(rawx)[i] != import(rawy)[i]) return false;
    return true;
  }


  RingElem RingFqVecImpl::mySymbolValue(const symbol& /*sym*/) const {CoCoA_THROW_ERROR("This ring has no symbols", "RingFqVecImpl::mySymbolValue"); return myZero();}


  ideal RingFqVecImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom RingFqVecImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    // No need to check compatibility -- it was checked when theta and phi were built
//???    return RingHom(new InducedHomImpl(QuotientRing(this), phi(theta(myQuotientingHomCtor()))));
    return myCompose(theta,phi);// BUG BUG BUG just to keep compiler quiet
  }


  bool RingFqVecImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }



  RingElem RingFqVecImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    const RingElemAlias x = indet(myBaseRingValue, 0);
    const vector<long> c = myImpl.myExport(import(rawx));
    RingElem ans = zero(myBaseRingValue);
    for (int i=0; i < myDeg; ++i)
      ans += c[i]*power(x,i);
    return ans;
  }


  void RingFqVecImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
// Reduce univariate poly mod generator of ideal
// convert to vector of coeffs
    RingElemAlias x(myBaseRingValue, rawarg);
    RingElem rem = NR(x, gens(myDefiningIdeal));
    myAssignZero(rawimage);
    for (SparsePolyIter it = BeginIter(rem); !IsEnded(it); ++it)
    {
      BigInt c; IsInteger(c, coeff(it));
      import(rawimage)[deg(PP(it))] = myImpl.myFpArith().myReduce(c);
    }
  }


  RingHom RingFqVecImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    // Compatibility has already been checked (see InducedHom in QuotientRing.C)
//???    CoCoA_ASSERT(IsZero(InducingHom(myModulus)));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  RingFqVecImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom):
      RingHomBase(domain, codomain(InducingHom)),
      myInducingHom(InducingHom)
  {
    CoCoA_ASSERT(IsInKer(DefiningIdeal(domain), InducingHom)); // actually already checked in pseudo-ctor InducedHom (in QuotientRing.C:875)
  }

  void RingFqVecImpl::InducedHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    myCodomain->myAssign(rawimage, raw(myInducingHom(QuotientRing(myDomain)->myCanonicalRepr(rawarg))));
  }


    ring NewRingFqVec(const MachineInt& p, const MachineInt& d)
    {
      if (IsNegative(p) || !IsPrime(p)) CoCoA_THROW_ERROR(ERR::BadSmallFpChar, "NewRingFq");
      if (IsNegative(d) || AsSignedLong(d) < 2) CoCoA_THROW_ERROR(ERR::BadArg, "NewRingFq");
      const long P = AsSignedLong(p);
      const long deg = AsSignedLong(d);
      ring Fpx = NewPolyRing(NewRingFp(P), NewSymbols(1));
  std::vector<SmallFpImpl::value> v = FindGroupGenerator(P,deg);
  RingElem MinPoly = VectorToPoly(Fpx, v);
  // RingElem MinPoly(Fpx);
  // RingElem x = indet(Fpx,0);
  // for (int i=0; i <= deg; ++i) MinPoly += v[i] * power(x,i);
  ideal I(MinPoly);
  return ring(new RingFqVecImpl(I));
    }

    ring NewRingFqVec(const BigInt& p, const MachineInt& d)
    {
      return NewRingFqVec(ConvertTo<long>(p), d);
    }
    ring NewRingFqVec(const MachineInt& p, const BigInt& d)
    {
      return NewRingFqVec(p, ConvertTo<long>(d));
    }
    ring NewRingFqVec(const BigInt& p, const BigInt& d)
    {
      return NewRingFqVec(ConvertTo<long>(p), ConvertTo<long>(d));
    }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingFqVec.C,v 1.20 2023/03/08 19:43:57 abbott Exp $
// $Log: RingFqVec.C,v $
// Revision 1.20  2023/03/08 19:43:57  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.19  2022/10/25 17:29:17  abbott
// Summary: Added FieldGen (finite field) and mem fn myFieldGen (default defn gives error)
//
// Revision 1.18  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.17  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.16  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.15  2021/02/17 17:50:36  abbott
// Summary: Now asserts ostream is in decimal mode (redmine 1547)
//
// Revision 1.14  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.13  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.12  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.11  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.10  2019/03/04 10:58:38  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.9  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.8  2018/05/18 16:42:11  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.7  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.6  2018/05/17 15:42:47  bigatti
// -- added include SparsePolyIter
//
// Revision 1.5  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.4  2018/02/27 10:56:39  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.3  2017/09/06 11:56:29  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.2  2016/01/27 14:00:26  abbott
// Summary: Added new pseudo-ctors which accept BigInts
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//
