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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingZZ.C,v 1.42 2024/07/23 19:54:13 bigatti Exp $
// $Log: RingZZ.C,v $
// Revision 1.42  2024/07/23 19:54:13  bigatti
// Summary: new error syntax: CoCoA_THROW_ERROR1/2
//
// Revision 1.41  2023/03/08 19:43:57  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.40  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.39  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.38  2021/10/29 19:47:59  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.37  2021/01/07 15:16:53  abbott
// Summary: Corrected copyright
//
// Revision 1.36  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.35  2020/02/12 09:01:47  bigatti
// -- changed myTestIsMaximal etc to return void (and consequences)
//
// Revision 1.34  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.33  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.32  2019/10/15 11:54:08  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.31  2019/03/04 11:11:54  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.30  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.29  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.28  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.27  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.26  2018/03/20 11:38:07  bigatti
// -- changed iAm***Test --> myTestIs***;  and it returns bool
//
// Revision 1.25  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.24  2018/02/27 10:56:21  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.23  2017/11/08 14:05:24  abbott
// Summary: Corrected myIsIrred
//
// Revision 1.22  2017/09/06 11:56:30  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.21  2016/09/08 14:12:32  bigatti
// -- mySetMaximalFlag --> myAssignMaximalFlag
// -- mySetPrimeFlag --> myAssignPrimeFlag
// -- updated the related code
// -- (still "old design", but better aligned)
//
// Revision 1.20  2016/03/25 20:06:49  abbott
// Summary: Introduced new mem fns myCeil and myNearestInt; impls in ZZ, QQ, & TwinFloat. Renamed NearestInteger to NearestInt.
//
// Revision 1.19  2015/07/29 11:04:54  bigatti
// -- added space after comma in printing rings
//
// Revision 1.18  2014/07/30 14:09:43  abbott
// Summary: Changed myAmbientRing --> myRing
// Author: JAA
//
// Revision 1.17  2014/07/11 15:49:55  bigatti
// -- added  myOutputSelfLong
//
// Revision 1.16  2014/07/08 13:14:41  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.15  2014/06/14 19:45:07  abbott
// Summary: Added new fn CmpAbs for RingElem (via memfn myCmpAbs)
// Author: JAA
//
// Revision 1.14  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.13  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.12  2014/03/27 17:18:11  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.11  2014/03/27 16:41:57  bigatti
// -- added myMul
// -- fixed intersect
//
// Revision 1.10  2014/01/28 09:45:32  abbott
// Now uses ConvertTo instead of IsInteger.
//
// Revision 1.9  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.8  2013/06/17 08:55:45  abbott
// Improved RingZZImpl::IdealImpl::myDivMod (but there's a design issue!!)
//
// Revision 1.7  2013/02/21 14:14:41  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.6  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.5  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.4  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.3  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.2  2012/04/27 15:03:16  abbott
// Added mem fn IamFiniteField;  corrected several Z into ZZ.
//
// Revision 1.1  2012/02/10 12:15:53  bigatti
// -- was RingZ.C
//
// Revision 1.28  2012/02/10 10:30:03  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.27  2012/02/08 13:43:50  bigatti
// -- changed Z --> ZZ
//
// Revision 1.26  2012/01/26 17:13:34  abbott
// Cleaned up impl of  RingZImpl::IdealImpl::myAdd
//
// Revision 1.25  2012/01/26 15:02:33  abbott
// Replaced call to copy(..., back_inserter...) by a call to range insert
// [see Scott Meyers "Effective STL" section 5, page 26]
//
// Revision 1.24  2011/11/09 14:29:37  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.23  2011/09/06 15:28:58  abbott
// Changed impl of  "myCmp" so that its return value is in {-1,0,+1}
//
// Revision 1.22  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.21  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.20  2011/06/23 16:04:46  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.19  2011/05/24 14:54:29  abbott
// Consquential changes from removing several ctors for principal ideals.
//
// Revision 1.18  2011/05/19 14:46:47  abbott
// Added defn of myIsDouble.
//
// Revision 1.17  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.16  2010/10/27 20:58:45  abbott
// Major reorganization of GlobalManager and GMPAllocator.
//
// Revision 1.15  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.14  2010/10/01 15:45:17  bigatti
// -- added mySymbolValue
//
// Revision 1.13  2010/09/30 14:35:34  abbott
// Added new pseudo-ctor for ZEmbeddingHom which requires the domain as first arg:
// this form is essential for the ctor for RingQImpl!
//
// Added (hidden) new fn which checks whether RingZ has ref count greater than 3:
// this fn tells us whether any object other than the GlobalManager and RingQ
// refers to RingZ.  The fn is called only from the dtor of GlobalManager when
// checking that no CoCoALib values still survive.
//
// Revision 1.12  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.11  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.10  2009/10/29 18:41:42  abbott
// Made a change to  RingZImpl::IdealImpl::myDivMod  but I don't recall why now
// (the new version does avoid some wasteful copying).
//
// Revision 1.9  2009/10/26 15:42:55  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.8  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.7  2009/06/05 12:14:55  abbott
// Major change:
//   created new files NumTheory.H/C  which contain basic number theory operations
//   removed several basic number theory operations from ZZ.H/C
//   removed gcd from MachineInt.H/C
//   changed names of some basic fns:
//      IsPPrime -> IsProbPrime
//      invmod -> InvMod    (changed signature too)
//      powermod -> PowerMod  (changed signature too)
//   added new fns
//      NextProbPrime & PrevProbPrime
//   consequent changes to other code and tests and examples
//
// Revision 1.6  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.5  2008/11/20 10:49:17  abbott
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
// Revision 1.12  2007/03/07 14:08:11  bigatti
// -- minor: commented argument names for -Wextra
//
// Revision 1.11  2007/03/05 21:25:57  cocoa
// Forgot to check these in a few minutes ago.
//
// Revision 1.10  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.9  2007/01/17 12:32:39  cocoa
// Changed a few more "raw" variable names so that the code compiles fine
// also when CoCoA_DEBUG is set.
//
// Revision 1.8  2007/01/15 16:04:48  cocoa
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
// Revision 1.5  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.4  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
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
// Revision 1.8  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.7  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.6  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.5  2006/03/21 09:43:13  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
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
// Revision 1.4  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
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
// Revision 1.16  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.15  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.14  2004/11/09 15:53:01  cocoa
// -- fixed IamPrime  (now check also for equality to ideal(0))
//
// Revision 1.13  2004/11/04 18:47:42  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.12  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.11  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.10  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.9  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.8  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.7  2004/01/30 14:07:10  cocoa
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
// Revision 1.6  2004/01/28 15:49:06  cocoa
// RingZ now explicitly uses mpz_t to represent its values
// (previously it used ZZ objects, with unhappy ownership questions).
//
// Revision 1.5  2003/11/14 13:06:05  cocoa
// -- New function "myIsPrintAtom" for printing polynomials and fractions
//
// Revision 1.4  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.3  2003/10/09 14:55:19  cocoa
// - minor debugging after merge
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.14  2003/06/23 16:41:19  abbott
// Prior to public release:
//   changed class name to follow coding conventions
//   updated characteristic following RingBase
//   added IsIntegralDomain following RingBase
//
// Revision 1.13  2003/05/14 17:15:42  abbott
// Several consequential changes stemming from revisions to ring.H;
// most notably added support for homomorphisms and ideals.
//
// Revision 1.12  2002/11/15 11:26:06  abbott
// Updated following the renaming in ring.H.
// Commented out the exception unsafe IsZeroAddMul, and opted to
// use the default definition supplied by AbstractRing.
// Fixed a memory leak in IsDivisible, and made it exception safe.
//
// Revision 1.11  2002/07/05 15:21:45  abbott
// Revised ctor and dtor since zero is now an auto_ptr.
// Added definition of member function IsDivisible.
// Member function div now calls mpz_divexact.
//
// Revision 1.10  2002/06/22 17:10:57  abbott
// Changed name of "equal" member function to "IsEqual" (as per new ring.H).
//
// Revision 1.9  2002/06/03 13:36:10  abbott
// Added function definitions required to satisfy the OrderedRing interface.
//
// Revision 1.8  2002/05/30 13:50:29  abbott
// Added IsGCDDomain, IsField, and zero functions as required by the new ring.H.
//
// Revision 1.7  2002/03/21 15:13:42  bigatti
// - added: IsZeroAddMul
// - corrected: pool name
//
// Revision 1.6  2002/02/15 11:57:16  bigatti
// - added: gcd
// - corrected: equal
//
// Revision 1.5  2001/12/07 18:22:30  abbott
// Changed names in accordance with the new coding conventions.
//
// Revision 1.4  2001/11/23 20:59:09  abbott
// Added assignment from a long.
//
// Revision 1.3  2001/11/19 20:10:03  abbott
// Added using commands for the identifiers required from the standard
// library (required by g++ version 3.0.2, and more properly correct).
//
// Revision 1.2  2001/11/06 16:44:41  abbott
// Two main changes:
//   *  alignment with the use of a union for a ring::raw_elem
//      (instead of the previous void*);
//
//   *  the elements of a ring_Z are now pointers to mpz_t
//      (instead of pointer to objects of class ZZ)
//
// Revision 1.1  2001/10/31 20:37:15  abbott
// Initial revision
//
