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

#include "CoCoA/RingFqLog.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SmallFqLogImpl.H"
#include "CoCoA/SmallFqUtils.H"
#include "CoCoA/SparsePolyIter.H" // for myReduction
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H

#include <iostream>
using std::ostream;
#include <memory>
using std::unique_ptr;


namespace CoCoA
{

  class RingFqLogImpl: public QuotientRingBase
  {
  private: // data members
//???    typedef SmallFpImpl::value value_t;
    typedef FFqLogImpl::repr_t value_t;
    const long myDeg;
    const FFqLogImpl myImpl;
//    const value_t myModulus;
//    const SmallFpImpl myImpl;
    mutable MemPool myMemMgr;       // MemPool must come *BEFORE* myZeroPtr and myOnePtr
    unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    unique_ptr<RingElem> myFieldGenPtr;   ///< Every ring stores its own cyclic generator.

  private: // auxiliary functions
///???    static value_t PrincipalGen(const ideal& I); // used for arg checking in ctor
    static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(RingElemConstRawPtr rawx);

  private:
//?????    RingFqLogImpl(long p, int d); // called only by NewRingFqLog
    explicit RingFqLogImpl(const ideal& I);
    ~RingFqLogImpl();
    friend ring NewRingFqLog(const MachineInt& p, const MachineInt& d); // create field of size p^d
    friend ring NewRingFqLog(const BigInt& p, const MachineInt& d); // create field of size p^d
    friend ring NewRingFqLog(const MachineInt& p, const BigInt& d); // create field of size p^d
    friend ring NewRingFqLog(const BigInt& p, const BigInt& d); // create field of size p^d
  public: // disable copy ctor & assignment
    RingFqLogImpl(const RingFqLogImpl&) = delete;
    RingFqLogImpl& operator=(const RingFqLogImpl&) = delete;
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


  inline RingFqLogImpl::value_t& RingFqLogImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const RingFqLogImpl::value_t& RingFqLogImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }


  // // Returns generator of I as a value_t; returns 0 if value is too large to fit.
  // RingFqLogImpl::value_t RingFqLogImpl::PrincipalGen(const ideal& I)
  // {
  //   if (IsZero(I)) return 0;
  //   const BigInt GenI = ConvertTo<BigInt>(TidyGens(I)[0]);
  //   value_t p;
  //   if (!IsConvertible(p, GenI))  // check that the value of the principal generator will fit
  //     return 0;
  //   return p;
  // }


  RingFqLogImpl::RingFqLogImpl(const ideal& I):
      QuotientRingBase(RingOf(I), I),
      myDeg(deg(gens(I)[0])),
      myImpl(I),
      myMemMgr(myDeg*SmallFpImpl::ourDatumSize, "RingFqLogImpl.myMemMgr"),
      myZeroPtr(),
      myOnePtr(),
      myFieldGenPtr()
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myFieldGenPtr.reset(new RingElem(ring(this)));
    import(raw(*myFieldGenPtr)) = myImpl.myGen();
    myRefCountZero();
  }


  RingFqLogImpl::~RingFqLogImpl()
  {}


  BigInt RingFqLogImpl::myCharacteristic() const
  {
    return BigInt(myImpl.myModulus);
  }


  long RingFqLogImpl::myLogCardinality() const
  {
    return myImpl.myExtnDeg(); ///???BUG??? myImpl.myExtnDeg() * LogCardinality(BaseField);???
  }


  bool RingFqLogImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingFqLogImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingFqLogImpl::IamField() const
  {
    return true;
  }


  bool RingFqLogImpl::IamFiniteField() const
  {
    return true;
  }


  bool RingFqLogImpl::IamExact() const
  {
    return true;
  }


  ConstRefRingElem RingFqLogImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingFqLogImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingFqLogImpl::myNew() const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myZero();
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFqLogImpl::myNew(const MachineInt& n) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myReduce(LeastNNegRemainder(n,myImpl.myModulus));
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFqLogImpl::myNew(const BigInt& N) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = myImpl.myReduce(LeastNNegRemainder(N,myImpl.myModulus));
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingFqLogImpl::myNew(ConstRawPtr rawy) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    *ans = import(rawy);
    return RingElemRawPtr(ans);
  }


  void RingFqLogImpl::myDelete(RawPtr rawx) const
  {
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingFqLogImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    std::swap(import(rawx), import(rawy));
  }


  void RingFqLogImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingFqLogImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = myImpl.myReduce(LeastNNegRemainder(n,myImpl.myModulus));
  }


  void RingFqLogImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = myImpl.myReduce(LeastNNegRemainder(N,myImpl.myModulus));
  }


  void RingFqLogImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs) = myImpl.myZero();
  }


  void RingFqLogImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingFqLogImpl::myRecvTwinFloat");
  }


  void RingFqLogImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = myImpl.myNegate(import(rawx));
  }


  void RingFqLogImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myAdd(import(rawx), import(rawy));
  }


  void RingFqLogImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.mySub(import(rawx), import(rawy));
  }


  void RingFqLogImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    import(rawlhs) = myImpl.myMul(import(rawx), import(rawy));
  }


  void RingFqLogImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    import(rawlhs) = myImpl.myDiv(import(rawx), import(rawy));
  }


  bool RingFqLogImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) return false;
    import(rawlhs) = myImpl.myDiv(import(rawx), import(rawy));
    return true;
  }


  bool RingFqLogImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingFqLogImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void RingFqLogImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const  // assumes n > 1
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    import(rawlhs) = myImpl.myPower(import(rawx), n);
  }

  void RingFqLogImpl::myPowerBigExp(RawPtr rawlhs, ConstRawPtr rawx, const BigInt& N) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(N > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    // Use Fermat's Little Theorem to reduce exponent...
    import(rawlhs) = myImpl.myPower(import(rawx), N%(myImpl.myCard-1));  //??? BUG better to do it inside FFqLogImpl?
  }


  void RingFqLogImpl::myOutput(ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    CoCoA_ASSERT(IsDecimal(out));

    // BUG BUG BUG  lazy STOPGAP IMPL!!!
    const long p = myImpl.myModulus;
    const value_t x = import(rawx);
    long v = myImpl.myExport(x);
    out << '[';
    for (int i=0; i < myDeg-1; ++i)
    {
      out << v%p << ' ';
      v /= p;
    }
    CoCoA_ASSERT(v < p);
    out << v << ']';
  }


  bool RingFqLogImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return false; ///??? NYI probably OK if image of an int?
//    return myImpl.myExport(import(rawx)) >= 0;
  }


  bool RingFqLogImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return false; ///??? BUG BUG BUG  NYI  STOPGAP impl
//    return myImpl.myExport(import(rawx)) < 0;
  }


  void RingFqLogImpl::myOutputSelf(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    //    out << "FFp(" << myModulus << ")";
    out << "RingWithID(" << myID << ",\"FFqLog(p=" << myImpl.myModulus << ", deg=" << myDeg << ")\")";
  }


  void RingFqLogImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    // OMOut->mySendApplyStart();
    // OMOut << OpenMathSymbol("setname2", "GFp");
    // OMOut << myModulus;
    // OMOut->mySendApplyEnd();
  }


  void RingFqLogImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingFqLogImpl::myOutput");
    OMOut << 0; //???myImpl.myExport(import(rawx));
  }


  bool RingFqLogImpl::myIsZero(ConstRawPtr rawx) const
  {
    return myImpl.myIsZero(import(rawx));
  }


  bool RingFqLogImpl::myIsOne(ConstRawPtr rawx) const
  {
    return myImpl.myIsOne(import(rawx));
  }


  bool RingFqLogImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return myImpl.myIsMinusOne(import(rawx));
  }


  bool RingFqLogImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    const long v = myImpl.myExport(import(rawx));
    if (v >= myImpl.myModulus) return false;
    N = v;
    return true;
  }


  bool RingFqLogImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    const long v = myImpl.myExport(import(rawx));
    if (v >= myImpl.myModulus) return false;
    Q = v;
    return true;
  }


  bool RingFqLogImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    const long v = myImpl.myExport(import(rawx));
    if (v >= myImpl.myModulus) return false;
    d = v;
    return true;
  }


  // bool RingFqLogImpl::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  // {
  //   return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  // }


  // bool RingFqLogImpl::myIsZeroAddMul(RawPtr rawlhs, RawPtr /*rawtmp*/, ConstRawPtr rawfact1, ConstRawPtr rawfact2) const
  // { // same as above: just to avoid calling RingBase::myIsZeroAddMul with 4 args
  //   return myImpl.myIsZeroAddMul(import(rawlhs), import(rawfact1), import(rawfact2));
  // }


  bool RingFqLogImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return (import(rawx) == import(rawy));
  }


  RingElem RingFqLogImpl::mySymbolValue(const symbol& /*sym*/) const {CoCoA_THROW_ERROR("This ring has no symbols", "RingFqLogImpl::mySymbolValue"); return myZero();}


  ideal RingFqLogImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom RingFqLogImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    // No need to check compatibility -- it was checked when theta and phi were built
    CoCoA_THROW_ERROR(ERR::NYI, "RingFqLogImpl::myCompose");
//???    return RingHom(new InducedHomImpl(QuotientRing(this), phi(theta(myQuotientingHomCtor()))));
    return theta;// BUG BUG BUG just to keep compiler quiet
  }


  bool RingFqLogImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }



  RingElem RingFqLogImpl::myCanonicalRepr(ConstRawPtr rawx) const
  {
    const RingElemAlias x = indet(myBaseRingValue, 0);
    const long p = myImpl.myModulus;
    long c = myImpl.myExport(import(rawx));
    RingElem ans = zero(myBaseRingValue);
    for (int i=0; i < myDeg; ++i)
    {
      long a = c%p;
      ans += a*power(x,i);
      c /= p;
    }
    return ans;
  }


  void RingFqLogImpl::myReduction(RawPtr rawimage, ConstRawPtr rawarg) const
  {
// Reduce univariate poly mod generator of ideal
// convert to vector of coeffs
    RingElemAlias x(myBaseRingValue, rawarg);
    RingElem rem = NR(x, gens(myDefiningIdeal));
//    myAssignZero(rawimage);
    value_t ans = myImpl.myZero();
    for (SparsePolyIter it = BeginIter(rem); !IsEnded(it); ++it)
    {
//      BigInt c; IsInteger(c, coeff(it));
      const long c = ConvertTo<long>(coeff(it));  /// BUG BUG ??? overflow???
      ans = myImpl.myAdd(ans, myImpl.myMul(myImpl.myReduce(c),
                                           myImpl.myPower(myImpl.myGen(),deg(PP(it)))));
    }
    import(rawimage) = ans;
  }


  RingHom RingFqLogImpl::myInducedHomCtor(const RingHom& InducingHom) const
  {
    // Compatibility has already been checked (see InducedHom in QuotientRing.C)
//???    CoCoA_ASSERT(IsZero(InducingHom(myModulus)));
    return RingHom(new InducedHomImpl(QuotientRing(this), InducingHom));
  }


  RingFqLogImpl::InducedHomImpl::InducedHomImpl(const QuotientRing& domain, const RingHom& InducingHom):
      RingHomBase(domain, codomain(InducingHom)),
      myInducingHom(InducingHom)
  {
    CoCoA_ASSERT(IsInKer(DefiningIdeal(domain), InducingHom)); // actually already checked in pseudo-ctor InducedHom (in QuotientRing.C:875)
  }

  void RingFqLogImpl::InducedHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  {
    myCodomain->myAssign(rawimage, raw(myInducingHom(QuotientRing(myDomain)->myCanonicalRepr(rawarg))));
  }


  ///////////////////////////////////////////////////////
    ring NewRingFqLog(const MachineInt& p, const MachineInt& d)
    {
      if (IsNegative(p) || !IsPrime(p)) CoCoA_THROW_ERROR(ERR::BadSmallFpChar, "NewRingFqLog");
      if (IsNegative(d) || AsSignedLong(d) < 2) CoCoA_THROW_ERROR(ERR::BadArg, "NewRingFqLog");
      const long P = AsSignedLong(p);
      const long deg = AsSignedLong(d);
      ring Fpx = NewPolyRing(NewRingFp(P),NewSymbols(1));
  std::vector<SmallFpImpl::value> v = FindGroupGenerator(P,deg);
  RingElem MinPoly = VectorToPoly(Fpx, v);
  // RingElem MinPoly(Fpx);
  // RingElem x = indet(Fpx,0);
  // for (int i=0; i <= deg; ++i) MinPoly += v[i] * power(x,i);
  ideal I(MinPoly);
  return ring(new RingFqLogImpl(I));
    }

    ring NewRingFqLog(const BigInt& p, const MachineInt& d)
    {
      return NewRingFqLog(ConvertTo<long>(p), d);
    }
    ring NewRingFqLog(const MachineInt& p, const BigInt& d)
    {
      return NewRingFqLog(p, ConvertTo<long>(d));
    }
    ring NewRingFqLog(const BigInt& p, const BigInt& d)
    {
      return NewRingFqLog(ConvertTo<long>(p), ConvertTo<long>(d));
    }



} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingFqLog.C,v 1.23 2023/03/08 19:43:57 abbott Exp $
// $Log: RingFqLog.C,v $
// Revision 1.23  2023/03/08 19:43:57  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.22  2022/10/25 17:29:17  abbott
// Summary: Added FieldGen (finite field) and mem fn myFieldGen (default defn gives error)
//
// Revision 1.21  2022/08/05 20:45:50  abbott
// Summary: Added two NYI comments
//
// Revision 1.20  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.19  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.18  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.17  2021/02/17 17:50:36  abbott
// Summary: Now asserts ostream is in decimal mode (redmine 1547)
//
// Revision 1.16  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.15  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.14  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.13  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.12  2019/03/18 11:18:08  abbott
// Summary: Added ERR::NYI to myCompose (removed code which triggered compiler warning)
//
// Revision 1.11  2019/03/04 10:58:03  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.10  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.9  2018/05/18 16:42:11  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.8  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.7  2018/05/17 15:42:47  bigatti
// -- added include SparsePolyIter
//
// Revision 1.6  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.5  2018/02/27 10:56:55  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.4  2017/09/06 11:56:29  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.3  2016/01/27 14:48:44  abbott
// Summary: Removed debugging printout
//
// Revision 1.2  2016/01/27 14:00:26  abbott
// Summary: Added new pseudo-ctors which accept BigInts
//
// Revision 1.1  2015/12/18 15:25:07  abbott
// Summary: Added impls of non-prime finite fields
//
//
