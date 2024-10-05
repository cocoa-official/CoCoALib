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
    if (GlobalManager::ourGlobalDataPtr == 0)
      CoCoA_THROW_ERROR(ERR::NoGlobalMgr, "RingQQ()");
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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingQQ.C,v 1.38 2024/03/08 19:59:59 abbott Exp $
// $Log: RingQQ.C,v $
// Revision 1.38  2024/03/08 19:59:59  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.37  2022/08/05 12:41:02  abbott
// Summary: Minor tidying
//
// Revision 1.36  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.35  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.34  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.33  2021/01/31 10:02:20  abbott
// Summary: Removed old (& slightly buggy) impl of mpq_get_d_2exp; now in utils-gmp
//
// Revision 1.32  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.31  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.30  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.29  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.28  2019/03/04 10:59:04  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.27  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.26  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.25  2017/09/06 11:56:29  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.24  2017/05/09 13:51:30  bigatti
// -- changed BadRingHomArg2 --> BadPartialRingHomArg
//
// Revision 1.23  2016/03/25 20:42:37  abbott
// Summary: Renamed utils_gmp to utils-gmp
//
// Revision 1.22  2016/03/25 20:06:49  abbott
// Summary: Introduced new mem fns myCeil and myNearestInt; impls in ZZ, QQ, & TwinFloat. Renamed NearestInteger to NearestInt.
//
// Revision 1.21  2015/11/23 18:22:26  abbott
// Summary: Changed arg name GlobalZ into GlobalZZ
//
// Revision 1.20  2015/07/29 11:04:55  bigatti
// -- added space after comma in printing rings
//
// Revision 1.19  2014/07/28 16:05:18  abbott
// Summary: Renamed myEmbeddingHom to myEmbeddingHomCtor
// Author: JAA
//
// Revision 1.18  2014/07/28 15:50:59  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.17  2014/07/11 15:37:00  bigatti
// -- added myOutputSelfShort and myOutputSelfLong
//
// Revision 1.16  2014/07/09 11:37:15  abbott
// Summary: Removed several fns which are now in FractionFieldBase
// Author: JAA
//
// Revision 1.15  2014/07/08 15:23:16  abbott
// Summary: Updated to fit in with new FractionFieldBase
// Author: JAA
//
// Revision 1.14  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.13  2014/07/08 08:36:43  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.12  2014/06/14 19:45:07  abbott
// Summary: Added new fn CmpAbs for RingElem (via memfn myCmpAbs)
// Author: JAA
//
// Revision 1.11  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.10  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.9  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.8  2013/05/21 11:32:09  abbott
// Replaced mem fns FractionFieldBase::myGetNum (and Den) by FractionFieldBase::myRawNum (and Den).
//
// Revision 1.7  2013/02/21 14:14:42  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.6  2012/10/24 13:36:54  abbott
// Corrected alignment of a comment.
//
// Revision 1.5  2012/05/28 10:35:32  abbott
// Changed default defn of IsTrueGCDDomain (makes RingQQImpl a bit simpler).
//
// Revision 1.4  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.3  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.2  2012/04/27 15:03:39  abbott
// Added mem fn IamFiniteField
//
// Revision 1.1  2012/02/10 12:16:13  bigatti
// -- was RingQ.C
//
// Revision 1.23  2012/02/10 10:30:03  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.22  2012/02/08 13:45:18  bigatti
// -- changed Z,Q --> ZZ,QQ
//
// Revision 1.21  2011/11/09 14:11:58  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.20  2011/09/06 15:28:58  abbott
// Changed impl of  "myCmp" so that its return value is in {-1,0,+1}
//
// Revision 1.19  2011/08/24 10:28:49  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.18  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.17  2011/08/03 09:05:44  abbott
// Minor change to silence a mistaken compiler warning.
//
// Revision 1.16  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.15  2011/05/19 14:46:07  abbott
// Added defn of myIsDouble.
//
// Revision 1.14  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.13  2010/11/02 15:32:42  bigatti
// -- fixed myIsPrintAtom
//
// Revision 1.12  2010/10/27 20:58:45  abbott
// Major reorganization of GlobalManager and GMPAllocator.
//
// Revision 1.11  2010/09/30 14:42:17  abbott
// Changed copyright date.
//
// Revision 1.10  2010/09/30 14:41:32  abbott
// Added new fn to check whether RingQ has ref count greater than 1: this will
// happen only if RingQ is referred to by objects other than the GlobalManager.
// The fn is called only from the dtor of GlobalManager.
//
// RingQImpl ctor now increments its ref count at the start: this is necessary
// for exception cleanliness (see doc for details).
//
// Minor changes to MakeUniqueInstanceOfRingQ which must now take RingZ as arg
// (since no global RingZ exists when the fn is called during construction of
// the GlobalManager).
//
// Revision 1.9  2009/10/29 18:33:02  abbott
// Changed order of include directives (now alphabetical by file name).
// Removed some unnecessary std:: prefixes.
//
// Revision 1.8  2009/10/26 15:40:23  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.7  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2008/11/20 10:49:04  abbott
// Added definition of myFloor member fn.
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.2  2007/03/28 10:06:13  abbott
// Now gives error when you use RingZ() or RingQ() without creating GlobalManager.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.13  2007/03/05 21:25:57  cocoa
// Forgot to check these in a few minutes ago.
//
// Revision 1.12  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.11  2007/02/26 15:00:54  bigatti
// -- added placeholders for new syntax based on unique Z implementation
//
// Revision 1.10  2007/01/17 12:32:39  cocoa
// Changed a few more "raw" variable names so that the code compiles fine
// also when CoCoA_DEBUG is set.
//
// Revision 1.9  2007/01/15 16:08:26  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.8  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.7  2006/11/27 14:26:24  cocoa
// -- reorganised #include files
//
// Revision 1.6  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.5  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.4  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.3  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.2  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
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
// Revision 1.3  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.2  2005/06/22 14:42:16  cocoa
// Renamed MemPool data member to myMemMgr
// (seems more sensible than myMemory).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.14  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.13  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.12  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.11  2004/11/05 15:34:33  cocoa
// Consequential change following from the renaming of
// FieldIdealImpl and the introduction of the new pseudo-ctor.
//
// Revision 1.10  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.9  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.8  2004/07/14 16:40:42  cocoa
// Separated RingFpLog from its implementation which now resides in
// a new class: SmallFpLogImpl.  This is analogous to the change made
// to RingFp yesterday.
//
// Some tidying and other sundry minor changes.
//
// Revision 1.7  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.6  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.5  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.4  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.3  2004/01/30 14:07:10  cocoa
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
// Revision 1.2  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.1  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
//
