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
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingFpLogImpl::myRecvTwinFloat");
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
    if (!IsZZ(RingOf(I))) CoCoA_THROW_ERROR(ERR::IdealNotInRing, "NewRingFpLog(I)");
    return QuotientRing(new RingFpLogImpl(I));
  }


  bool IsGoodForRingFpLog(const MachineInt& p)
  {
    if (IsNegative(p) || !IsSignedLong(p)) return false;
    const long n = AsSignedLong(p);
    return SmallFpLogImpl::IsGoodCtorArg(n);
  }

  bool IsGoodForRingFpLog(const BigInt& P)
  {
    if (P <= 0) return false;
    long p;
    if (!IsConvertible(p, P)) return false;
    return IsGoodForRingFpLog(p);
  }

  bool IsGoodForRingFpLog(const ideal& I)
  {
    if (!IsZZ(RingOf(I))) return false;
    if (IsZero(I)) return false;
    return IsGoodForRingFpLog(ConvertTo<BigInt>(TidyGens(I)[0]));
  }


} // end of namespace CoCoA


// RCS header/log
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingFpLog.C,v 1.55 2023/03/08 19:43:57 abbott Exp $
// $Log: RingFpLog.C,v $
// Revision 1.55  2023/03/08 19:43:57  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.54  2022/10/25 17:29:17  abbott
// Summary: Added FieldGen (finite field) and mem fn myFieldGen (default defn gives error)
//
// Revision 1.53  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.52  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.51  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.50  2021/03/04 21:03:46  abbott
// Summary: enum revision and renaming (redmine 894)
//
// Revision 1.49  2021/02/17 17:50:36  abbott
// Summary: Now asserts ostream is in decimal mode (redmine 1547)
//
// Revision 1.48  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.47  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.46  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.45  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.44  2019/09/25 13:33:34  abbott
// Summary: Changed "compact" printed form of finite fields from FFp(...)  to ZZ/(...)
//
// Revision 1.43  2019/03/04 10:57:28  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.42  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.41  2017/09/06 11:56:29  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.40  2014/07/30 14:09:05  abbott
// Summary: Changed name AmbientRing --> RingOf
// Author: JAA
//
// Revision 1.39  2014/07/28 16:04:56  abbott
// Summary: Renamed myQuotientingHom to myQuotientingHomCtor
// Author: JAA
//
// Revision 1.38  2014/07/28 15:50:06  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.37  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.36  2014/07/02 16:35:04  bigatti
// -- new way of printing ring with ID
//
// Revision 1.35  2014/06/17 10:13:26  abbott
// Summary: Added (void)(phi) to avoid compiler warning about unused param
// Author: JAA
//
// Revision 1.34  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.33  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.32  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.31  2014/01/28 10:02:48  abbott
// Replaced some calls to IsInteger by calls to ConvertTo<BigInt>.
//
// Revision 1.30  2013/03/25 17:04:19  abbott
// Major clean-up of interface to SmallFpImpl/SmallFpLogImpl/SmallFpDoubleImpl
// (underlying impl remains much the same).  Removed lots of cruft.
// Consequential changes to RingFp* classes; small change to SparsePolyRing.
//
// Revision 1.29  2012/09/07 15:21:13  abbott
// First stage of revision of SmallFpImpl interface (and SmallFpLog, SmallFpDouble).
//
// Revision 1.28  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.27  2012/04/27 15:04:10  abbott
// Added mem fns IamFiniteField & myLogCardinality.
//
// Revision 1.26  2012/02/10 10:29:07  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.25  2012/02/08 15:13:30  bigatti
// -- changed: Z,Q -> ZZ,QQ
//
// Revision 1.24  2012/01/30 23:17:51  abbott
// Realigned errors
//
// Revision 1.23  2012/01/25 13:34:27  bigatti
// -- added myIsZeroAddMul with 4 args
// -- some tidying
//
// Revision 1.22  2011/11/09 14:11:58  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.21  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.20  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.19  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.18  2011/05/24 14:54:29  abbott
// Consquential changes from removing several ctors for principal ideals.
//
// Revision 1.17  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.16  2011/05/20 09:45:04  abbott
// Harmonized RingFp, RingFpLog, RingFpDouble -- they are now almost ready to be merged into a single class!
//
// Revision 1.15  2011/05/19 14:38:27  abbott
// Updated small prime finite field impls to allow user to specify
// separately for each whether to use symmetric or non-negative
// residues for export operations (myExport and printing).
//
// Revision 1.14  2011/03/22 20:00:37  abbott
// Added IsGoodForXXX fns to test whether a given arg is suitable as
// characteristic for the given type of small prime finite field.
//
// Revision 1.13  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.12  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.11  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.10  2009/09/25 14:01:11  abbott
// Cleaned up include directives.
//
// Revision 1.9  2009/09/24 16:22:44  abbott
// Tidied up include directives.  Removed some unnecessary "std::" prefixes.
//
// Revision 1.8  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.7  2009/06/05 12:08:27  abbott
// Changed return type of operator%(ZZ,MachineInteger); it is now unsigned long
// instead of ZZ.
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInteger in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.3  2007/05/21 12:44:25  abbott
// Changed impl of powering routine as a consequence of changed signature
// to the modulus operator for ZZ modulo machine integer.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.13  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.12  2007/03/08 10:23:29  bigatti
// -- CanonHom --> CanonicalHom
//
// Revision 1.11  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.10  2007/03/03 14:07:23  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.9  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.8  2007/01/15 16:11:41  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.7  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.6  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.5  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.4  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/10/06 10:15:52  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.7  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.6  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.5  2006/04/21 15:01:36  cocoa
// Changed default implementation of RingBase::myGcd -- it now gives a SERIOUS
// error.  All fields must now handle a call to gcd explicitly: they can use
// the new myGcdInField function.  It's now cleaner than it was.
//
// Revision 1.4  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.3  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.2  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.6  2005/10/14 15:25:07  cocoa
// Major tidying and cleaning to small prime finite fields.
// Several consequential changes.  Improved their documentation.
//
// Added Makefile and script to include/CoCoA/ directory to
// keep library.H up to date.
//
// Revision 1.5  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.4  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.3  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
//
// Revision 1.2  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/02/11 16:45:24  cocoa
// Removed the useless and misleading functions myInit and myKill
// from the SmallFp*Impl classes; various consequential changes.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.16  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.15  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.14  2004/11/09 15:47:29  cocoa
// -- changed myOutput: FF --> FFp
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
// Revision 1.11  2004/11/02 15:09:45  cocoa
// -- added name to MemPool
//
// Revision 1.10  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.9  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
// Revision 1.8  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.7  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.6  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.5  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.4  2004/01/30 14:07:10  cocoa
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
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:15:45  cocoa
// New coding convention for rings.
//
// Revision 1.14  2003/06/23 16:55:43  abbott
// Minor cleaning prior to public release.
// Name change and alignment with RingBase.
//
// Revision 1.13  2003/04/23 10:16:46  abbott
// Numerous changes, almost all consequential from changes made to ring_FpLog.H
//  * ring_FpLog is now derived from QuotientRing (so has a different ctor,
//    and must implement the functions specific to QuotientRings)
//  * implementation of ring_FpLog::hom added
//  * uses the generic FieldIdeal for ideals
//  * now uses reset member function of auto_ptr instead of assignment from raw ptr
//  * ring_FpLog::power is now simpler (since exponent must be non-negative)
//
// Revision 1.12  2002/11/14 18:22:00  abbott
// Consequential changes following renaming in ring.H.
// Slightly improved some error messages.
// Changed division by zero check into an assertion.
//
// Revision 1.11  2002/07/05 15:24:12  abbott
// Added definition of member function IsDivisible.
//
// Revision 1.10  2002/06/28 12:10:06  abbott
// Now using C++ vectors and auto_ptr to achieve exception safety.
// Better argument checking in the constructor.
//
// Revision 1.9  2002/06/22 17:08:18  abbott
// Changed name of "equal" member function to "IsEqual" (as per new ring.H).
//
// Revision 1.8  2002/06/04 13:04:02  abbott
// Changed type of elements in the log/exp tables; now the type is FpTableElem
// (which is a typedef in ring_FpLog.H)
//
// Revision 1.7  2002/05/30 14:31:29  abbott
// Added definitions for member functions IsField and zero (as required by the new ring.H).
// Removed bogus function ngens.
//
// Revision 1.6  2002/05/15 15:43:47  abbott
// Fixed bug in power for negative exponents.
//
// Revision 1.5  2002/05/15 15:25:31  abbott
// Added definitions for new member functions of ring_FpLog:
//   characteristic, negate, and power.
// Other minor changes to stay in line with new ring.H
//
// Revision 1.4  2002/03/19 13:59:47  abbott
// Changed class name to ring_FpLog.  Added two "using" directives.
//
// Revision 1.3  2002/02/08 11:20:14  bigatti
// - changed syntax to IsZeroAddMul
//
// Revision 1.2  2002/01/31 11:01:37  abbott
// Corrected the constructor in the case of characteristic 2.
//
// Revision 1.1  2002/01/30 17:26:27  abbott
// Initial revision
//
