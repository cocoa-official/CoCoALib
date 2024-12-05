//   Copyright (c)  2005-2011,2014,2021  John Abbott and Anna M. Bigatti

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

#include "CoCoA/RingZZ.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector; // for RingZZImpl::output and RingZZImpl::IdealImpl


namespace CoCoA
{

  // This is the concrete class which does all the work.
  class RingZZImpl: public RingBase
  {
  private:
    typedef mpz_t value_t; // mpz_t is the actual type of the values in a RingZZImpl
     static value_t& import(RingElemRawPtr rawx);
    static const value_t& import(RingElemConstRawPtr rawx);
  private: // data members
    mutable MemPool myMemMgr;           // MemPool must come BEFORE myZeroPtr and myOnePtr
    std::unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    std::unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.

  private:
    RingZZImpl();                        ///< Called only by NewRingZZ
    RingZZImpl(const RingZZImpl&) = delete;           // copy ctor disabled
    RingZZImpl operator=(const RingZZImpl&) = delete; // assignment disabled
    ~RingZZImpl();
    friend ring MakeUniqueInstanceOfRingZZ(); // the only function allowed to call the constructor
    friend bool RingZZStillInUse(const ring& ZZ);

    static const mpz_t& AsMPZ(ConstRefRingElem r); // handy in a few member fns

  public: // functions that every ring must implement
    const ring& myBaseRing() const override;
    BigInt myCharacteristic() const override;
    bool IamCommutative() const override;
    bool3 IamIntegralDomain3(bool) const override;
    bool IamOrderedDomain() const override;
    bool IamField() const override;
    bool IamFiniteField() const override;
    bool IamExact() const override;
    ConstRefRingElem myZero() const override;
    ConstRefRingElem myOne() const override;
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& N) const override;
    RingElemRawPtr myNew(ConstRawPtr rawcopy) const override;
    void myDelete(RawPtr rawx) const override;                      // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;                          // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;                 // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;              // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;                  // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                               // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;                 // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                          // true iff x is invertible
    bool myIsIrred(ConstRawPtr rawx) const override;                               ///< true iff x is irreducible (if defined)
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // lhs = gcd(x,y) if TrueGCDDomain;
    void myExgcd(RawPtr rawlhs, RawPtr rawxcofac, RawPtr rawycofac, ConstRawPtr rawx, ConstRawPtr rawy) const;            // lhs = gcd(x,y) = xcofac*x + ycofac*y  if TrueGCDDomain;
    //    void myNormalizeFrac: use default defn
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;  // lhs = x^n, n>1, x not -1,0,1
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;             // out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;
    void myOutputSelf(std::ostream& out) const override      { out << "ZZ"; }
    void myOutputSelfShort(std::ostream& out) const override { out << "ZZ"; }
    void myOutputSelfLong(ostream& out) const override;
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                       // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;         // OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                                // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                                 // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                            // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                  // true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;                     // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                   // false iff x overflows
    //    virtual bool myIsZeroAddMul: use default definition
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;             // x == y
    int myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const override;                  // result is <0, =0, >0 according as x<y, x=y, x>y
    int myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const override;               // equiv to myCmp(abs(x),abs(y))
    int mySign(ConstRawPtr rawx) const override;                                   // -1,0,+1 according as x <0,=0,>0
    BigInt myFloor(ConstRawPtr rawx) const override;                               // identity fn
    BigInt myCeil(ConstRawPtr rawx) const override;                                // identity fn
    BigInt myNearestInt(ConstRawPtr rawx) const override;                          // identity fn

    RingElem mySymbolValue(const symbol& /*sym*/) const override  {CoCoA_THROW_ERROR("This ring has no symbols", "RingZZImpl::mySymbolValue"); return myZero();}
    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))

    bool myImageLiesInSubfield(const RingHom& phi) const override;

    // Function special to RingZZBase
    RingHom myZZEmbeddingHomCtor(const ring& codomain) const;

  private: // homomorphism class
    class HomImpl: public RingHomBase
    {
    public:
      HomImpl(const ring& ZZ, const ring& codomain);
      void myApply(RawPtr rawimage, ConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return false; }
    };

  private: // ideal class

    class IdealImpl: public IdealBase
    {
    public:
      IdealImpl(const ring& ZZ, const std::vector<RingElem>& gens);
      // Default copy ctor works fine.
      // Assignment disabled (see below)
      virtual ~IdealImpl();
    public: // disable assignment
      IdealImpl& operator=(const IdealImpl&) = delete;
    public:
      IdealImpl* myClone() const override;
//???    virtual void swap(ideal& other);

    public: // functions every ideal must have
      const ring& myRing() const override;
      bool IamZero() const override;
      bool IamOne() const override;
      bool IhaveElem(RingElemConstRawPtr rawx) const override;
      void myReduceMod(RingElemRawPtr rawx) const override;
      void myAdd(const ideal&) override;
      void myMul(const ideal&) override;
      void myIntersect(const ideal&) override;
      void myColon(const ideal&) override;
      bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const override; // true iff quotient exists & is unique, sets lhs = num/den modulo I iff result is true

      const std::vector<RingElem>& myGens() const override; // gens as specified by user
      const std::vector<RingElem>& myTidyGens(const CpuTimeLimit&) const override; // tidier set of gens

    protected:
      void myTestIsMaximal() const override;
      void myTestIsPrimary() const override;
      void myTestIsPrime() const override;
      void myTestIsRadical() const override;

    private: // data members
      const ring myR;
      std::vector<RingElem> myGensValue;
      std::vector<RingElem> myTidyGensValue;
      static const IdealImpl* GetPtr(const ideal& I);
    };
  };



  inline RingZZImpl::value_t& RingZZImpl::import(RingElemRawPtr rawx)
  {
    return *static_cast<value_t*>(rawx.myRawPtr());
  }

  inline const RingZZImpl::value_t& RingZZImpl::import(RingElemConstRawPtr rawx)
  {
    return *static_cast<const value_t*>(rawx.myRawPtr());
  }


  inline const mpz_t& RingZZImpl::AsMPZ(ConstRefRingElem x)
  {
    CoCoA_ASSERT(IsZZ(owner(x)));
    return RingZZImpl::import(raw(x));
  }



  RingZZImpl::RingZZImpl():
      myMemMgr(sizeof(value_t), "RingZZImpl.myMemMgr"),
      myZeroPtr(),
      myOnePtr()
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myRefCountZero();
  }


  RingZZImpl::~RingZZImpl()
  {}


  const ring& RingZZImpl::myBaseRing() const
  {
    CoCoA_THROW_ERROR(ERR::BadRing, "BaseRing");
    return RingZZ();
  }


  BigInt RingZZImpl::myCharacteristic() const
  {
    return BigInt(0);
  }


  bool RingZZImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingZZImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingZZImpl::IamOrderedDomain() const
  {
    return true;
  }


  bool RingZZImpl::IamField() const
  {
    return false;
  }


  bool RingZZImpl::IamFiniteField() const
  {
    return false;
  }


  bool RingZZImpl::IamExact() const
  {
    return true;
  }


  ConstRefRingElem RingZZImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingZZImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingZZImpl::myNew() const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpz_init(*ans);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingZZImpl::myNew(const MachineInt& n) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    if (IsNegative(n))
      mpz_init_set_si(*ans, AsSignedLong(n));
    else
      mpz_init_set_ui(*ans, AsUnsignedLong(n));
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingZZImpl::myNew(const BigInt& N) const
  {
    value_t* ans = static_cast<value_t*>(myMemMgr.alloc());
    mpz_init_set(*ans, mpzref(N));
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingZZImpl::myNew(ConstRawPtr rawcopy) const
  {
    value_t *ans = static_cast<value_t*>(myMemMgr.alloc());
    mpz_init_set(*ans, import(rawcopy));
    return RingElemRawPtr(ans);
  }


  void RingZZImpl::myDelete(RawPtr rawx) const
  {
    mpz_clear(import(rawx)); // EXPLICIT dtor call here
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingZZImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    mpz_swap(import(rawx), import(rawy));
  }


  void RingZZImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    mpz_set(import(rawlhs), import(rawx));
  }


  void RingZZImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    if (IsNegative(n))
      mpz_set_si(import(rawlhs), AsSignedLong(n));
    else
      mpz_set_ui(import(rawlhs), AsUnsignedLong(n));
  }


  void RingZZImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    mpz_set(import(rawlhs), mpzref(N));
  }


  void RingZZImpl::myAssignZero(RawPtr rawlhs) const
  {
    mpz_set_ui(import(rawlhs), 0);
  }


  void RingZZImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingZZImpl::myRecvTwinFloat");
  }

  void RingZZImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    mpz_neg(import(rawlhs), import(rawx));
  }


  void RingZZImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpz_add(import(rawlhs), import(rawx), import(rawy));
  }


  void RingZZImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpz_sub(import(rawlhs), import(rawx), import(rawy));
  }


  void RingZZImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpz_mul(import(rawlhs), import(rawx), import(rawy));
  }


  void RingZZImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    mpz_divexact(import(rawlhs), import(rawx), import(rawy));
  }


  bool RingZZImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) return false;
    BigInt remainder;
    mpz_tdiv_qr(import(rawlhs), mpzref(remainder), import(rawx), import(rawy));
    return IsZero(remainder);
  }


  bool RingZZImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return myIsOne(rawx) || myIsMinusOne(rawx);
  }


  bool RingZZImpl::myIsIrred(ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!myIsZero(rawx));
//    if (myIsZero(rawx)) CoCoA_THROW_ERROR(ERR::ZeroRingElem, "IsIrred");
    return IsPrime(abs(BigIntFromMPZ(import(rawx))));
  }


  void RingZZImpl::myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpz_gcd(import(rawlhs), import(rawx), import(rawy));
  }


  void RingZZImpl::myExgcd(RawPtr rawlhs, RawPtr rawxcofac, RawPtr rawycofac, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    mpz_gcdext(import(rawlhs), import(rawxcofac), import(rawycofac), import(rawx), import(rawy));
  }


  void RingZZImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    mpz_pow_ui(import(rawlhs), import(rawx), n);
  }


  void RingZZImpl::myOutput(ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    vector<char> buffer(mpz_sizeinbase(import(rawx),10)+2);
    mpz_get_str(&buffer[0], 10, import(rawx));
    out << &buffer[0];
  }


  bool RingZZImpl::myIsPrintAtom(ConstRawPtr rawx) const
  {
    return mySign(rawx)>=0;
  }


  BigInt RingZZImpl::myFloor(ConstRawPtr rawx) const
  {
    BigInt N;
    myIsInteger(N, rawx);
    return N;
  }

  BigInt RingZZImpl::myCeil(ConstRawPtr rawx) const
  {
    BigInt N;
    myIsInteger(N, rawx);
    return N;
  }

  BigInt RingZZImpl::myNearestInt(ConstRawPtr rawx) const
  {
    BigInt N;
    myIsInteger(N, rawx);
    return N;
  }


  bool RingZZImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return mySign(rawx)<0;
  }


  void RingZZImpl::myOutputSelfLong(ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "RingWithID(" << myID << ", \"ZZ\")";
  }


  void RingZZImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut << OpenMathSymbol("setname1", "Z");
  }


  void RingZZImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    OMOut << BigIntFromMPZ(import(rawx));  // have to get value of x as a BigInt
  }


  bool RingZZImpl::myIsZero(ConstRawPtr rawx) const
  {
    return mpz_sgn(import(rawx)) == 0;
  }


  bool RingZZImpl::myIsOne(ConstRawPtr rawx) const
  {
    return mpz_cmp_ui(import(rawx), 1) == 0;
  }


  bool RingZZImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return mpz_cmp_si(import(rawx), -1) == 0;
  }


  bool RingZZImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    mpz_set(mpzref(N), import(rawx));
    return true;
  }


  bool RingZZImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    mpq_set_z(mpqref(Q), import(rawx));
    return true;
  }

  bool RingZZImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    long exp;
    d = mpz_get_d_2exp(&exp, import(rawx)); // BUG, ignore possible overflow in exp
    if (numeric_limits<double>::radix != 2) CoCoA_THROW_ERROR(ERR::NYI, "RingZZImpl::myIsDouble");
    if (exp < numeric_limits<double>::min_exponent) { d=0; return true; }  // ???false also for underflow???
    if (exp > numeric_limits<double>::max_exponent) return false;
    d = ldexp(d,exp);
    return true;
  }


  bool RingZZImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return mpz_cmp(import(rawx), import(rawy)) == 0;
  }


  int RingZZImpl::myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return sign(mpz_cmp(import(rawx), import(rawy)));
  }


  int RingZZImpl::myCmpAbs(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return sign(mpz_cmpabs(import(rawx), import(rawy)));
  }


  int RingZZImpl::mySign(ConstRawPtr rawx) const
  {
    return mpz_sgn(import(rawx)); // Needs modification for GMP older than 4.1
  }


  ideal RingZZImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new RingZZImpl::IdealImpl(ring(this), gens));
  }


  //---------------------------------------------------------------------------
  // Functions to do with ring homomorphisms

  RingHom RingZZImpl::myZZEmbeddingHomCtor(const ring& codomain) const
  {
    if (IsZZ(codomain)) return IdentityHom(ring(this));
    return RingHom(new RingZZImpl::HomImpl(ring(this), codomain));
  }


  RingHom RingZZImpl::myCompose(const RingHom& phi, const RingHom& /*theta*/) const
  {
    return myZZEmbeddingHomCtor(codomain(phi));
  }


  bool RingZZImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return false;
  }


  //---------------------------------------------------------------------------
  // Homomorphisms

  RingZZImpl::HomImpl::HomImpl(const ring& ZZ, const ring& codomain):
      RingHomBase(ZZ, codomain)
  {}


  void RingZZImpl::HomImpl::myApply(RawPtr rawimage, ConstRawPtr rawarg) const
  {
    // Unfortunately we have to copy the value of arg.
    BigInt N;  // wasteful new/delete???
    myDomain->myIsInteger(N, rawarg);  // must necessarily succeed
    myCodomain->myAssign(rawimage, N);
  }


  //---------------------------------------------------------------------------
  // Ideals

  inline const RingZZImpl::IdealImpl* RingZZImpl::IdealImpl::GetPtr(const ideal& I)
  {
    return dynamic_cast<const RingZZImpl::IdealImpl*>(I.myIdealPtr());
  }


  RingZZImpl::IdealImpl::IdealImpl(const ring& ZZ, const std::vector<RingElem>& gens):
      myR(ZZ),
      myGensValue(gens),
      myTidyGensValue()
  {
    RingElem g(ZZ);
    for (long i=0; i < len(gens); ++i) // break out if g=1???
      g = gcd(g, gens[i]);
    if (IsZero(g)) return;
    myTidyGensValue.push_back(g);
  }


  RingZZImpl::IdealImpl::~IdealImpl()
  {}


  RingZZImpl::IdealImpl* RingZZImpl::IdealImpl::myClone() const
  {
    return new IdealImpl(*this);
  }


//???  void RingZZImpl::IdealImpl::swap(CoCoA::ideal& J) {}


  const ring& RingZZImpl::IdealImpl::myRing() const
  {
    return myR;
  }


  bool RingZZImpl::IdealImpl::IamZero() const
  {
    return myTidyGensValue.empty();
  }


  bool RingZZImpl::IdealImpl::IamOne() const
  {
    return !myTidyGensValue.empty() && IsOne(myTidyGensValue[0]);
  }


  void RingZZImpl::IdealImpl::myReduceMod(RawPtr rawx) const
  {
    if (IamZero()) return;
    mpz_fdiv_r(import(rawx), import(rawx), AsMPZ(myTidyGensValue[0]));
  }


  bool RingZZImpl::IdealImpl::IhaveElem(const RingElemConstRawPtr rawx) const
  {
    if (myR->myIsZero(rawx)) return true;
    if (IamZero()) return false;
    RingElem junk(myR);
    return myR->myIsDivisible(raw(junk), rawx, raw(myTidyGensValue[0]));
  }


  void RingZZImpl::IdealImpl::myAdd(const ideal& J)
  {
    const IdealImpl* J1 = GetPtr(J);
    // ANNA: clever insert (skip 0s, deal with 1s)
    myGensValue.insert(myGensValue.end(), gens(J).begin(), gens(J).end());
    if (IsZero(J)) return;
    if (IamZero())
    {
      myTidyGensValue.push_back(J1->myTidyGensValue[0]);
      IamPrime3Flag = IamMaximal3Flag = J1->IamMaximal3Flag;
      return;
    }
    myTidyGensValue[0] = gcd(myTidyGensValue[0], J1->myTidyGensValue[0]);
    IamPrime3Flag = IamMaximal3Flag = uncertain3;
  }


  void RingZZImpl::IdealImpl::myMul(const ideal& J)
  {
    if (IamZero()) return;
    //    CoCoA_THROW_ERROR(ERR::NYI, "RingZZImpl::IdealImpl::mul");
    if (IsZero(J)) { myGensValue.clear(); myTidyGensValue.clear(); IamPrime3Flag = true3; IamMaximal3Flag = false3; return; }
    const IdealImpl* J1 = GetPtr(J);
    myTidyGensValue[0] = myTidyGensValue[0]*J1->myTidyGensValue[0];
    myGensValue = myTidyGensValue;
    IamPrime3Flag = IamMaximal3Flag = uncertain3;
  }


  void RingZZImpl::IdealImpl::myIntersect(const ideal& J)
  {
    if (IamZero() || IsOne(J)) return;
    if (IsZero(J)) { myGensValue.clear(); myTidyGensValue.clear(); IamPrime3Flag = true3; IamMaximal3Flag = false3; return; }
    const IdealImpl* J1 = GetPtr(J);
    myTidyGensValue[0] = lcm(myTidyGensValue[0], J1->myTidyGensValue[0]);
    myGensValue = myTidyGensValue;
    IamPrime3Flag = IamMaximal3Flag = uncertain3;
  }


  void RingZZImpl::IdealImpl::myColon(const ideal& J)
  {
    if (IsZero(J))
    {
      if (IamZero()) myTidyGensValue.push_back(one(myR));
      else myTidyGensValue[0] = 1;
      myGensValue = myTidyGensValue;
      IamPrime3Flag = IamMaximal3Flag = false3;
      return;
    }
    if (IamZero()) return;
    const IdealImpl* J1 = GetPtr(J);
    myTidyGensValue[0] = myTidyGensValue[0]/gcd(myTidyGensValue[0], J1->myTidyGensValue[0]);
    myGensValue = myTidyGensValue;
    IamPrime3Flag = IamMaximal3Flag = uncertain3;
  }


  bool RingZZImpl::IdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
    const ring& ZZ = myRing();
    if (IamOne()) { ZZ->myAssignZero(rawlhs); return true; }
    if (IamZero())
    {
      return ZZ->myIsDivisible(rawlhs, rawnum, rawden);
    }

    const BigInt m = ConvertTo<BigInt>(myTidyGensValue[0]);
    BigInt g, RecipDen, junk;
    mpz_gcdext(mpzref(g), mpzref(RecipDen), mpzref(junk), import(rawden), mpzref(m));
    if (!IsOne(g)) { return false; } // (gcd != 1) ==>  quot not uniquely defd in this ring
//2013-06-15    BigInt quo, rem;
//2013-06-15    mpz_fdiv_qr(mpzref(quo), mpzref(rem), import(rawnum), mpzref(g));
//2013-06-15    if (!IsZero(rem)) { ZZ->myAssignZero(rawlhs); return; }
//???    if (!CoCoA::IsZero(rem)) CoCoA_THROW_ERROR(ERR::DivByZero, "RingZZImpl::IdealImpl::inverse");
//2013-06-15    mpz_mul(import(rawlhs), mpzref(quo), mpzref(RecipDen));
    mpz_mul(mpzref(junk), import(rawnum), mpzref(RecipDen));
    mpz_mod(import(rawlhs), mpzref(junk), mpzref(m));
    return true;
  }


  const std::vector<RingElem>& RingZZImpl::IdealImpl::myGens() const
  {
    return myGensValue;
  }


  const std::vector<CoCoA::RingElem>& RingZZImpl::IdealImpl::myTidyGens(const CpuTimeLimit&) const
  {
    return myTidyGensValue;
  }


  void RingZZImpl::IdealImpl::myTestIsMaximal() const
  {
    myAssignMaximalFlag(!IamZero() && !IamOne() && mpz_probab_prime_p(AsMPZ(myTidyGensValue[0]), 25)); // trust to luck in randomized primality test
  }


  void RingZZImpl::IdealImpl::myTestIsPrimary() const
  { CoCoA_THROW_ERROR(ERR::NYI, "myTestIsPrimary ZZ"); }


  void RingZZImpl::IdealImpl::myTestIsPrime() const
  { myAssignPrimeFlag(IamZero() || IamMaximal()); }


  void RingZZImpl::IdealImpl::myTestIsRadical() const
  { CoCoA_THROW_ERROR(ERR::NYI, "myTestIsRadical ZZ"); }


  //---------------------------------------------------------------------------

  // This fn should be called only by ctor for CoCoA::GlobalManager (see GlobalManager.C)
  ring MakeUniqueInstanceOfRingZZ()
  { return ring(new RingZZImpl()); }


  const ring& RingZZ()
  {
    if (GlobalManager::ourGlobalDataPtr == nullptr)
      CoCoA_THROW_ERROR(ERR::NoGlobalMgr, "RingZZ()");
    return GlobalManager::ourGlobalDataPtr->myRingZZ();
  }


  bool IsZZ(const ring& R)
  { return dynamic_cast<const RingZZImpl*>(R.myRawPtr()) != nullptr; }


  RingHom ZZEmbeddingHom(const ring& codomain)
  {
    // Ugly down cast of raw ptr in RingZZ() to type RingZZImpl*.
    const RingZZImpl* const ZZPtr = static_cast<const RingZZImpl*>(RingZZ().myRawPtr());
    return ZZPtr->myZZEmbeddingHomCtor(codomain);
  }

  RingHom ZZEmbeddingHom(const ring& ZZ, const ring& codomain)
  {
    CoCoA_ASSERT(IsZZ(ZZ));
    // Ugly down cast of raw ptr in RingZZ() to type RingZZImpl*.
    const RingZZImpl* const ZZPtr = static_cast<const RingZZImpl*>(ZZ.myRawPtr());
    return ZZPtr->myZZEmbeddingHomCtor(codomain);
  }


  // This fn is used only in the dtor for GlobalManager.
  bool RingZZStillInUse(const ring& ZZ)
  {
    const RingZZImpl* ZZptr = dynamic_cast<const RingZZImpl*>(ZZ.myRawPtr());
#ifdef CoCoA_DEBUG
    if (ZZptr->myRefCount() > 3)
      std::cerr << "ERROR!!!  RingZZ refcount = " << ZZptr->myRefCount() << " but should be 3." << std::endl;
#endif
    return ZZptr->myRefCount() > 3; // copy in GlobalManager & as base ring of RingQQ & in RingQQ embedding hom
  }


} // end of namespace CoCoA
