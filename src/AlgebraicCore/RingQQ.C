//   Copyright (c)  2003-2011,2014,2021  John Abbott and Anna M. Bigatti

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

#include "CoCoA/RingQQ.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/utils-gmp.H"

#include <memory>
using std::unique_ptr;
#include <cmath>
using std::ldexp;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;


namespace CoCoA
{

  class RingQQImpl: public FractionFieldBase
  {
  private:
    typedef mpq_t value_t; // mpq_t is the actual type of the values in a RingQQImpl
    static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(const RingElemConstRawPtr rawx);

  private: // data members
    mutable MemPool myMemMgr;           // MemPool must come before myZeroPtr and myOnePtr
    unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.

  private:
    RingQQImpl(const ring& ZZ);
    RingQQImpl(const RingQQImpl&) = delete;           // copy ctor disabled
    RingQQImpl operator=(const RingQQImpl&) = delete; // assignment disabled
    ~RingQQImpl();
    friend FractionField MakeUniqueInstanceOfRingQQ(const ring&); // the only function allowed to call the constructor
    friend bool RingQQStillInUse(const FractionField& QQ);

  public:
    typedef RingElemRawPtr RawPtr;
    typedef RingElemConstRawPtr ConstRawPtr;

    // functions that every ring must implement
    ConstRefRingElem myZero() const override;
    ConstRefRingElem myOne() const override;
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& n) const override;
    RingElemRawPtr myNew(const BigRat& Q) const override;
    RingElemRawPtr myNew(ConstRawPtr rawt) const override;
    void myDelete(RawPtr rawx) const override;                      // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;                        // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;               // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;        // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;                // lhs = N
    void myAssign(RawPtr rawlhs, const BigRat& Q) const override;                    // lhs = Q
    void myAssignZero(RawPtr rawlhs) const override;                             // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;               // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                                // true iff x is invertible
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = gcd(x,y) in a field
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;// lhs = x^n, n>1, x not -1,0,1
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;                   // out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;                                 // x^n may be printed without parentheses
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;
    void myOutputSelf(std::ostream& out) const override       { out << "QQ"; }
    void myOutputSelfShort(std::ostream& out) const override  { out << "QQ"; }
    void myOutputSelfLong(std::ostream& out) const override;
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                       // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;         // OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                                // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                 // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                            // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                  // true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                 // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                   // false iff x overflows
    // Use default definition myIsZeroAddMul (see ring.C)
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;             // x == y
    int myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const override;                  // result is <0, =0, >0 according as x<y, x=y, x>y
    int myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const override;               // equiv to myCmp(abs(x),abs(y))
    int mySign(ConstRawPtr rawx) const override;                                   // -1,0,+1 according as x <0,=0,>0
    BigInt myFloor(ConstRawPtr rawx) const override;                               // floor(x)
    BigInt myCeil(ConstRawPtr rawx) const override;                                // ceil(x)
    BigInt myNearestInt(ConstRawPtr rawx) const override;                          // NearestInt(x)

    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))

    bool myImageLiesInSubfield(const RingHom& phi) const override;

    // functions every FractionField must implement
    ConstRawPtr myRawNum(ConstRawPtr rawq) const override;  ///< result belongs to BaseRing!!
    ConstRawPtr myRawDen(ConstRawPtr rawq) const override;  ///< result belongs to BaseRing!!
    RingHom myEmbeddingHomCtor() const override;
    RingHom myInducedHomCtor(const RingHom&) const override;


  private:
    class InducedHomImpl: public RingHomBase
    {
    public:
      InducedHomImpl(const FractionField& domain, const RingHom& phi);
      void myApply(RawPtr rawimage, ConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return !IsZero(characteristic(myCodomain)) || !IsField(myCodomain); }
    private:
      RingHom myInducingHomValue;
    };

  };


  // These two inline fns are used only in this file.
  inline RingQQImpl::value_t& RingQQImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const RingQQImpl::value_t& RingQQImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }



  RingQQImpl::RingQQImpl(const ring& ZZ):
      FractionFieldBase(ZZ),
      myMemMgr(sizeof(value_t), "RingQQImpl.myMemMgr")
  {
    CoCoA_ASSERT(IsZZ(ZZ));
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myRefCountZero();
  }


  RingQQImpl::~RingQQImpl()
  {}


  ConstRefRingElem RingQQImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingQQImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingQQImpl::myNew() const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpq_init(*ans);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingQQImpl::myNew(const MachineInt& n) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpq_init(*ans);
    if (IsNegative(n))
      mpq_set_si(*ans, AsSignedLong(n), 1);
    else
      mpq_set_ui(*ans, AsUnsignedLong(n), 1);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingQQImpl::myNew(const BigInt& N) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpq_init(*ans);
    mpq_set_z(*ans, mpzref(N));
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingQQImpl::myNew(const BigRat& Q) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpq_init(*ans);
    mpq_set(*ans, mpqref(Q));
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingQQImpl::myNew(ConstRawPtr rawcopy) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpq_init(*ans);
    mpq_set(*ans, import(rawcopy));
    return RingElemRawPtr(ans);
  }


  void RingQQImpl::myDelete(RawPtr rawx) const
  {
    mpq_clear(import(rawx));
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingQQImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    mpq_swap(import(rawx), import(rawy));
  }


  void RingQQImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    mpq_set(import(rawlhs), import(rawx));
  }


  void RingQQImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    if (IsNegative(n))
      mpq_set_si(import(rawlhs), AsSignedLong(n), 1);
    else
      mpq_set_ui(import(rawlhs), AsUnsignedLong(n), 1);
  }


  void RingQQImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    mpq_set_z(import(rawlhs), mpzref(N));
  }


  void RingQQImpl::myAssign(RawPtr rawlhs, const BigRat& Q) const
  {
    mpq_set(import(rawlhs), mpqref(Q));
  }


  void RingQQImpl::myAssignZero(RawPtr rawlhs) const
  {
    mpq_set_ui(import(rawlhs), 0, 1);
  }


  void RingQQImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingQQImpl::myRecvTwinFloat");
  }

  void RingQQImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    mpq_neg(import(rawlhs), import(rawx));
  }


  void RingQQImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpq_add(import(rawlhs), import(rawx), import(rawy));
  }


  void RingQQImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpq_sub(import(rawlhs), import(rawx), import(rawy));
  }


  void RingQQImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpq_mul(import(rawlhs), import(rawx), import(rawy));
  }


  void RingQQImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    mpq_div(import(rawlhs), import(rawx), import(rawy));
  }


  bool RingQQImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) return false;
    mpq_div(import(rawlhs), import(rawx), import(rawy));
    return true;
  }


  bool RingQQImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingQQImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myGcdInField(rawlhs, rawx, rawy);
  }


  void RingQQImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    mpz_pow_ui(mpq_numref(import(rawlhs)), mpq_numref(import(rawx)), n);
    mpz_pow_ui(mpq_denref(import(rawlhs)), mpq_denref(import(rawx)), n);
  }


  void RingQQImpl::myOutput(std::ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    constexpr int base = 10;
    const size_t BufferSize = 3 + mpz_sizeinbase(mpq_numref(import(rawx)), base)
                                + mpz_sizeinbase(mpq_denref(import(rawx)), base);
    vector<char> buffer(BufferSize);
    mpq_get_str(&buffer[0], base, import(rawx));
    out << &buffer[0];
  }


  bool RingQQImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return (mpz_cmp_si(mpq_denref(import(rawx)), 1)==0 && mpq_sgn(import(rawx)) >= 0);
  }


  bool RingQQImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return (mpq_sgn(import(rawx)) < 0);
  }


  void RingQQImpl::myOutputSelfLong(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "RingWithID(" << myID << ", \"QQ\")";
  }


  void RingQQImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut << OpenMathSymbol("setname1", "QQ");
  }


  void RingQQImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "NormalizedRational");
    OMOut << BigIntFromMPZ(mpq_numref(import(rawx)))
          << BigIntFromMPZ(mpq_denref(import(rawx)));
    OMOut->mySendApplyEnd();
  }


  bool RingQQImpl::myIsZero(ConstRawPtr rawx) const
  {
    return mpq_sgn(import(rawx)) == 0;
  }


  bool RingQQImpl::myIsOne(ConstRawPtr rawx) const
  {
    return mpq_cmp_ui(import(rawx), 1, 1) == 0;
  }


  bool RingQQImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return mpq_cmp_si(import(rawx), -1, 1) == 0;
  }


  bool RingQQImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(mpz_sgn(mpq_denref(import(rawx))) > 0);
    if (myIsZero(rawx)) { N = 0; return true; }
    if (mpz_cmp_ui(mpq_denref(import(rawx)), 1) != 0) return false;
    mpz_set(mpzref(N), mpq_numref(import(rawx)));
    return true;
  }


  bool RingQQImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    mpq_set(mpqref(Q), import(rawx));
    return true;
  }


  bool RingQQImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    long exp = 0; // pointless initialization to keep compiler quiet
    d = mpq_get_d_2exp(&exp, import(rawx)); // BUG, ignore possible overflow in exp
    if (numeric_limits<double>::radix != 2) CoCoA_THROW_ERROR(ERR::NYI, "RingQQImpl::myIsDouble");
    if (exp < numeric_limits<double>::min_exponent) { d=0; return true; }  // ???false also for underflow???
    if (exp > numeric_limits<double>::max_exponent) return false;
    d = ldexp(d,exp);
    return true;
  }


  bool RingQQImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return mpq_cmp(import(rawx), import(rawy)) == 0;
  }


  int RingQQImpl::myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return sign(mpq_cmp(import(rawx), import(rawy)));
  }


  int RingQQImpl::myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return sign(mpq_cmpabs(import(rawx), import(rawy)));
  }


  int RingQQImpl::mySign(ConstRawPtr rawx) const
  {
    return mpq_sgn(import(rawx));
  }


  BigInt RingQQImpl::myFloor(ConstRawPtr rawx) const
  {
    BigInt N;
    mpz_fdiv_q(mpzref(N), mpq_numref(import(rawx)), mpq_denref(import(rawx)));
    return N;
  }

  BigInt RingQQImpl::myCeil(ConstRawPtr rawx) const
  {
    BigInt N;
    mpz_cdiv_q(mpzref(N), mpq_numref(import(rawx)), mpq_denref(import(rawx)));
    return N;
  }

  BigInt RingQQImpl::myNearestInt(ConstRawPtr rawx) const
  {
    BigInt N;
    mpq_round(mpzref(N), import(rawx));
    return N;
  }


  ideal RingQQImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return NewFieldIdeal(ring(this), gens);
  }


  RingHom RingQQImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    return myInducedHomCtor(phi(theta(myEmbeddingHomCtor())));
  }

  bool RingQQImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true;
  }


  //--------------------------------------------------
  // Below are functions special to a FractionField

  // const ring& RingQQImpl::myBaseRing() const
  // {
  //   return myZZ;
  // }


  RingElemConstRawPtr RingQQImpl::myRawNum(ConstRawPtr rawq) const
  {
    return RingElemConstRawPtr(mpq_numref(import(rawq)));
  }

  RingElemConstRawPtr RingQQImpl::myRawDen(ConstRawPtr rawq) const
  {
    return RingElemConstRawPtr(mpq_denref(import(rawq)));
  }


  RingHom RingQQImpl::myEmbeddingHomCtor() const
  {
    return RingHom(ZZEmbeddingHom(myBaseRing(), FractionField(this)));
  }


  RingHom RingQQImpl::myInducedHomCtor(const RingHom& phi) const
  {
    return RingHom(new InducedHomImpl(FractionField(this), phi));
  }



  RingQQImpl::InducedHomImpl::InducedHomImpl(const FractionField& domain, const RingHom& phi):
      RingHomBase(domain, codomain(phi)),
      myInducingHomValue(phi)
  { /*??? MUST CHECK IT MAKES SENSE -- e.g.  given ZZ->ZZ[x] ker=0, but cannot do QQ->ZZ[x]*/ }


  void RingQQImpl::InducedHomImpl::myApply(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    RingHom phi = myInducingHomValue;
    const FractionField& QQ = myDomain;
    CoCoA_ASSERT(domain(phi) == BaseRing(QQ));
    const ring& ZZ = BaseRing(QQ);
    RingElemAlias N(ZZ, QQ->myRawNum(rawarg));
    RingElemAlias D(ZZ, QQ->myRawDen(rawarg));

    const RingElem ImN = phi(N);
    const RingElem ImD = phi(D);

    if (!myCodomain->myIsDivisible(rawimage, raw(ImN), raw(ImD)))
      CoCoA_THROW_ERROR(ERR::BadPartialRingHomArg, "RingQQImpl::InducedHomImpl::myApply");
  }



  // This fn is declared in GlobalManager.C and called only by ctor of CoCoA::GlobalManager.
  FractionField MakeUniqueInstanceOfRingQQ(const ring& GlobalZZ)
  {
    return FractionField(new RingQQImpl(GlobalZZ));
  }


  const FractionField& RingQQ()
  {
    if (!IsInitialized())  CoCoA_THROW_ERROR1(ERR::NoGlobalMgr);
    return GlobalManager::ourGlobalDataPtr->myRingQQ();
  }


  bool IsQQ(const ring& R)
  {
    return dynamic_cast<const RingQQImpl*>(R.myRawPtr()) != 0;
  }


  RingHom QQEmbeddingHom(const ring& codomain)
  {
    return InducedHom(RingQQ(), ZZEmbeddingHom(codomain));
  }


  // This fn is used only in the dtor for GlobalManager.
  bool RingQQStillInUse(const FractionField& QQ)
  {
    const RingQQImpl* QQptr = dynamic_cast<const RingQQImpl*>(QQ.myRawPtr());
#ifdef CoCoA_DEBUG
    if (QQptr->myRefCount() > 1)
      std::cerr << "ERROR!!!  RingQQ refcount = " << QQptr->myRefCount() << " but should be 1." << std::endl;
#endif
    return QQptr->myRefCount() > 1; // copy in GlobalManager
  }


} // end of namespace CoCoA
