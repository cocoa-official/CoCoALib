//   Copyright (c)  2004-2010,2014,2016,2021  John Abbott and Anna M. Bigatti

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

#include "CoCoA/RingTwinFloat.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/FieldIdeal.H"
#include "CoCoA/FloatApprox.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/ideal.H"
#include "CoCoA/verbose.H"

#include <iostream>
using std::ostream;
#include <memory>
using std::unique_ptr;
#include <vector>
using std::vector; // only in RingTwinFloatImpl::output
#include <algorithm>
using std::max;    // only in <anon>::BitAccuracy  and in RingTwinFloatImpl::myCmp
using std::min;    // only in RingTwinFloatImpl::myIsEqualNZIgnoreSign
#include <limits>
using std::numeric_limits;
#include <cmath>
//using std::log;   // only in myDebugPrint
//using std::fabs;  // only in BitAccuracy
//using std::ceil;  // only in RingTwinFloatImpl::myOutput  &  RingTwinFloatImpl::myDebugPrint
//using std::floor; // only in RingTwinFloatImpl::myOutput

namespace CoCoA
{

  RingElem RingTwinFloatBase::mySymbolValue(const symbol& /*sym*/) const
  {
    CoCoA_THROW_ERROR("This ring has no symbols", "RingTwinFloatBase::mySymbolValue");
    return myZero();
  }


  const RingTwinFloatBase* RingTwinFloatPtr(const ring& R)
  {
    return dynamic_cast<const RingTwinFloatBase*>(R.myRawPtr());
  }

  // const RingTwinFloatBase* RingTwinFloatPtr(const ring& R, const char* const FnName)
  // {
  //   const RingTwinFloatBase* ptr = RingTwinFloatPtr(R);
  //   if (ptr == nullptr) CoCoA_THROW_ERROR(ERR::NotRingTwinFloat, FnName);
  //   return ptr;
  // }

  const RingTwinFloatBase* RingTwinFloatPtr(const ring& R, const ErrorContext& ErrCtx)
  {
    const RingTwinFloatBase* ptr = RingTwinFloatPtr(R);
    if (ptr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::NotRingTwinFloat, ErrCtx);
    return ptr;
  }



  namespace // anonymous namespace for file local functions
  {

    inline long EnforceLimits(long lwb, long value, long upb)
    {
      if (value < lwb) return lwb;
      if (value > upb) return upb;
      return value;
    }


    typedef mpf_t* MultipleFloat_t;
    typedef mpf_t const* ConstMultipleFloat_t;

    inline MultipleFloat_t import(RingElemRawPtr rawx)
    {
      return static_cast<MultipleFloat_t>(rawx.myRawPtr());
    }

    inline ConstMultipleFloat_t import(RingElemConstRawPtr rawx)
    {
      return static_cast<ConstMultipleFloat_t>(rawx.myRawPtr());
    }


    // ASSUMES that a & b are non-zero and have the same sign.
    // Returns largest positive integer k (up to MaxLog) s.t.
    // 2^(-k) >= rho  where rho = |a-b|/|a| is the relative difference.  
    // If no such k exists, it returns 0.
    // Equivalently k = min(max(0,-1-ceil(log2(rho))), MaxLog)
    // NB this impl is off by 1 if rho is exactly a power of 1/2.
    long LogRelDiff(const mpf_t a, const mpf_t b, long MaxLog)
    {
      CoCoA_ASSERT(MaxLog > 0);
      CoCoA_ASSERT(mpf_sgn(a) != 0 && mpf_sgn(b) != 0);
      CoCoA_ASSERT(mpf_sgn(a) == mpf_sgn(b));

      mpf_t RelDiff; mpf_init2(RelDiff, 32); // a very low precision suffices
      mpf_reldiff(RelDiff, a, b);
      mpf_abs(RelDiff, RelDiff);
      // Get rid of two (rare?) awkward cases:
      if (mpf_sgn(RelDiff) == 0) { mpf_clear(RelDiff); return MaxLog; }
      if (mpf_cmp_ui(RelDiff, 1) >= 0) { mpf_clear(RelDiff); return 0; }

      long LogRelDiff;
      mpf_get_d_2exp(&LogRelDiff, RelDiff);
      // LogRelDiff is now 1+floor(log2(RelDiff))
      mpf_clear(RelDiff);
      if (LogRelDiff == 0) return 0;
      return min(-LogRelDiff, MaxLog);
    }


  } // end of anonymous namespace


  /*-----------------------------------------------------------------*/
  /** \include RingTwinFloat.txt  */
  /*-----------------------------------------------------------------*/
  class RingTwinFloatImpl: public RingTwinFloatBase
  {
  private: // data members
    mutable MemPool myMemMgr;         // MemPool must come before myZeroPtr, myOnePtr, and myMinusOnePtr
    unique_ptr<RingElem> myZeroPtr;     ///< Every ring stores its own zero.
    unique_ptr<RingElem> myOnePtr;      ///< Every ring stores its own one.
    unique_ptr<RingElem> myMinusOnePtr; ///< Useful for myIsMinusOne
    const long myAccuracyBits;        ///< the precision requested (or 10, whichever is the greater)
    const long myBufferBits;          ///< number of buffer bits (see article)
    const long myNoiseBits;           ///< initially the last myNoiseBits are "noise"
    const long myWorkingBits;         ///< the full precision actually used
    enum { NoOverlap, AllowOverlap } myEqTestMode; // used only in myIsEqualNZIgnoreSign -- see redmine 859
    mutable gmp_randstate_t myRandomSource;

  private:
    RingTwinFloatImpl(long AccuracyBits, long BufferBits, long NoiseBits);
    RingTwinFloatImpl(const RingTwinFloatImpl&) = delete;           ///< disable copy ctor
    RingTwinFloatImpl operator=(const RingTwinFloatImpl&) = delete; ///< disable assignment
    ~RingTwinFloatImpl();
    // There are the only two functions allowed to call the constructor.
    friend RingTwinFloat NewRingTwinFloat(const MachineInt& AccuracyBits);
    friend RingTwinFloat NewRingTwinFloat(const MachineInt& AccuracyBits, const MachineInt& BufferBits, const MachineInt& NoiseBits);
    static void ThrowInsuffPrec() { throw RingTwinFloat::InsufficientPrecision(); };

  public: // functions every ring must implement
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
    RingElemRawPtr myNew(const BigRat& Q) const override;
    RingElemRawPtr myNew(ConstRawPtr rawt) const override;
    void myDelete(RawPtr rawx) const override;                                   // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;                        // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;               // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;            // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;                // lhs = N
    void myAssign(RawPtr rawlhs, const BigRat& Q) const override;                // lhs = Q
    void myAssignZero(RawPtr rawlhs) const override;                             // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    MantExp2 myExport(ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;               // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x-y
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x*y
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                        // true iff x is invertible
    // No GCD virtual void myGcd(...) const override;
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;// lhs = x^n, n>1, x not -1,0,1
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;           // out << x
    void myOutputSelf(std::ostream& out) const override;                         // out << R
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                     // OMOut << R
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;       // OMOut << x
    void myDebugPrint(std::ostream& out, ConstRawPtr rawx) const;
    bool myIsZero(ConstRawPtr rawx) const override;                              // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                               // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                          // x == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;                // true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;               // true iff x is rational
    bool myIsDouble(double& d, ConstRawPtr rawx) const override;                 // false iff x overflows
    //    bool myIsZeroAddMul: use default definition
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;           // x == y
    int myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const override;                // result is <0, =0, >0 according as x<y, x=y, x>y
    int mySign(ConstRawPtr rawx) const override;                                 // -1,0,+1 according as x <0,=0,>0
    BigInt myFloor(ConstRawPtr rawx) const override;                             // floor(x)
    BigInt myCeil(ConstRawPtr rawx) const override;                              // ceil(x)
    BigInt myNearestInt(ConstRawPtr rawx) const override;                        // NearestInt(x)

    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override;  // phi(theta(...))
    RingHom myHomCtor(const ring& codomain) const; //  not override???
    bool myImageLiesInSubfield(const RingHom& phi) const override;

    long myPrecisionBits() const override;
  private: // impl details
    enum RelErrRating {TooSmall, OK, TooBig};
    RelErrRating myCheckRelErr(ConstRawPtr rawx) const;
    void myCheckValidity(RawPtr rawx) const;
    void myOuterSemiwidth(mpf_t outersw, const ConstMultipleFloat_t& x) const;
    void myInnerSemiwidth(mpf_t innersw, const ConstMultipleFloat_t& x) const;
    void myPerturb(MultipleFloat_t x) const;
    bool3 myIsEqualNZIgnoreSign(ConstRawPtr rawx, ConstRawPtr rawy) const; // test whether abs(x) == abs(y), may return "uncertain"
    bool myIsRational(BigInt& N, BigInt& D, ConstRawPtr rawx) const;    // does the work for public myIsRational

    friend RingHom NewApproxHom(const ring& TwinFloat, const ring& R);
    class HomExactImpl: public RingHomBase
    {
    public:
      HomExactImpl(const RingTwinFloat& domain, const ring& codomain);
      void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const override;
      bool IamPartial() const override  { return true; }
//???      bool IamPartial() const { return !IsZero(characteristic(myCodomain)) || !IsField(myCodomain); }
    };

    class HomApproxImpl: public RingHomBase
    {
    public:
      HomApproxImpl(const RingTwinFloat& domain, const ring& codomain);
      void myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const override;
      bool IamPartial() const override  { return true; }
    };
  };


  void RingTwinFloatImpl::myOuterSemiwidth(mpf_t outersw, const ConstMultipleFloat_t& x) const
  {
    mpf_sub(outersw, x[0], x[1]);
    mpf_abs(outersw, outersw);
  }

  void RingTwinFloatImpl::myInnerSemiwidth(mpf_t innersw, const ConstMultipleFloat_t& x) const
  {
    mpf_sub(innersw, x[0], x[1]);
    mpf_abs(innersw, innersw);
    mpf_div_2exp(innersw, innersw, myNoiseBits/2);
  }


  RingTwinFloatImpl::RingTwinFloatImpl(long AccuracyBits, long BufferBits, long NoiseBits):
      myMemMgr(2*sizeof(mpf_t), "RingTwinFloatImpl.myMemMgr"),
      myAccuracyBits(EnforceLimits(8,AccuracyBits,8388608)),     // force value between 8 and 8388608
      myBufferBits(EnforceLimits(8,BufferBits,8388608)),         // force value between 8 and 8388608
      myNoiseBits(2*EnforceLimits(32/2,(NoiseBits+1)/2,1024/2)), // force myNoiseBits to be even & between 32 and 1024
      myWorkingBits(myAccuracyBits+myBufferBits+myNoiseBits),
      myEqTestMode(AllowOverlap)
  {
    gmp_randinit_default(myRandomSource);
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myMinusOnePtr.reset(new RingElem(ring(this), -1));
    myRefCountZero();
  }


  RingTwinFloatImpl::~RingTwinFloatImpl()
  {
    gmp_randclear(myRandomSource);
  }


  const ring& RingTwinFloatImpl::myBaseRing() const
  {
    return RingQQ(); ///???????????????????????????
  }


  BigInt RingTwinFloatImpl::myCharacteristic() const
  {
    return BigInt(0);
  }


  bool RingTwinFloatImpl::IamCommutative() const
  {
    return true;
  }


  bool3 RingTwinFloatImpl::IamIntegralDomain3(bool) const
  {
    return true3;
  }


  bool RingTwinFloatImpl::IamOrderedDomain() const
  {
    return true;;
  }


  bool RingTwinFloatImpl::IamField() const
  {
    return true;
  }


  bool RingTwinFloatImpl::IamFiniteField() const
  {
    return false;
  }


  bool RingTwinFloatImpl::IamExact() const
  {
    return false;
  }


  ConstRefRingElem RingTwinFloatImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingTwinFloatImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingTwinFloatImpl::myNew() const
  {
    MultipleFloat_t ans = static_cast<MultipleFloat_t>(myMemMgr.alloc());
    mpf_init2(ans[0], myWorkingBits);
    mpf_init2(ans[1], myWorkingBits);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingTwinFloatImpl::myNew(const MachineInt& n) const
  {
    if (IsZero(n)) { return myNew(); }
    MultipleFloat_t ans = static_cast<MultipleFloat_t>(myMemMgr.alloc());
    mpf_init2(ans[0], myWorkingBits);
    mpf_init2(ans[1], myWorkingBits);

    if (IsNegative(n))
      mpf_set_si(ans[0], AsSignedLong(n));
    else
      mpf_set_ui(ans[0], AsUnsignedLong(n));
    myPerturb(ans);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingTwinFloatImpl::myNew(const BigInt& N) const
  {
    if (IsZero(N)) { return myNew(); }
    MultipleFloat_t ans = static_cast<MultipleFloat_t>(myMemMgr.alloc());
    mpf_init2(ans[0], myWorkingBits);
    mpf_init2(ans[1], myWorkingBits);

    mpf_set_z(ans[0], mpzref(N));
    myPerturb(ans);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingTwinFloatImpl::myNew(const BigRat& Q) const
  {
    if (IsZero(Q)) { return myNew(); }
    if (IsOneDen(Q)) { return myNew(num(Q)); }
    MultipleFloat_t ans = static_cast<MultipleFloat_t>(myMemMgr.alloc());
    mpf_init2(ans[0], myWorkingBits);
    mpf_init2(ans[1], myWorkingBits);

    mpf_set_q(ans[0], mpqref(Q));
    myPerturb(ans);
    return RingElemRawPtr(ans);
  }


  RingElemRawPtr RingTwinFloatImpl::myNew(ConstRawPtr rawcopy) const
  {
    MultipleFloat_t ans = static_cast<MultipleFloat_t>(myMemMgr.alloc());
    ConstMultipleFloat_t rhs(import(rawcopy));
    // NB cannot use mpf_init_set here as it clobbers the precision
    mpf_init2(ans[0], myWorkingBits);
    mpf_init2(ans[1], myWorkingBits);

    mpf_set(ans[0], rhs[0]);
    mpf_set(ans[1], rhs[1]);

    // Should I perturb the copy???
    return RingElemRawPtr(ans);
  }


  void RingTwinFloatImpl::myDelete(RawPtr rawx) const
  {
    MultipleFloat_t val(import(rawx));
    mpf_clear(val[1]);
    mpf_clear(val[0]);
    myMemMgr.free(rawx.myRawPtr());
  }


  void RingTwinFloatImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    MultipleFloat_t lhs(import(rawx));
    MultipleFloat_t rhs(import(rawy));
    mpf_swap(lhs[0], rhs[0]);
    mpf_swap(lhs[1], rhs[1]);
  }


  void RingTwinFloatImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    MultipleFloat_t x(import(rawlhs));
    ConstMultipleFloat_t y(import(rawx));
    mpf_set(x[0], y[0]);
    mpf_set(x[1], y[1]);
  }


  void RingTwinFloatImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    if (IsZero(n)) { myAssignZero(rawlhs); return; }
    MultipleFloat_t x(import(rawlhs));

    if (IsNegative(n))
      mpf_set_si(x[0], AsSignedLong(n));
    else
      mpf_set_ui(x[0], AsUnsignedLong(n));

    myPerturb(x);
  }


  void RingTwinFloatImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    if (IsZero(N)) { myAssignZero(rawlhs); return; }
    MultipleFloat_t x(import(rawlhs));
    mpf_set_z(x[0], mpzref(N));

    myPerturb(x);
  }


  void RingTwinFloatImpl::myAssign(RawPtr rawlhs, const BigRat& Q) const
  {
    if (IsZero(Q)) { myAssignZero(rawlhs); return; }
    if (IsOneDen(Q)) { myAssign(rawlhs, num(Q)); return; }

    MultipleFloat_t x(import(rawlhs));
    mpf_set_q(x[0], mpqref(Q));

    myPerturb(x);
  }


  void RingTwinFloatImpl::myAssignZero(RawPtr rawlhs) const
  {
    MultipleFloat_t x(import(rawlhs));
    mpf_set_ui(x[0], 0);
    mpf_set_ui(x[1], 0);
  }


  void RingTwinFloatImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    if (myCheckRelErr(rawx) == TooBig)  // CAREFUL HERE because rawx belongs to a DIFFERENT RingTwinFloat!!!
      ThrowInsuffPrec();
    myAssign(rawlhs, rawx);  // let GMP do the precision conversion
    myCheckValidity(rawlhs); // might do a myPerturb
  }


  void RingTwinFloatImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    MultipleFloat_t x(import(rawlhs));
    ConstMultipleFloat_t y(import(rawx));
    mpf_neg(x[0], y[0]);
    mpf_neg(x[1], y[1]);
  }


  void RingTwinFloatImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    // Eliminate the trivial cases of x=0 or y=0; simplifies later code.
    if (myIsZero(rawx)) { myAssign(rawlhs, rawy); return; }
    if (myIsZero(rawy)) { myAssign(rawlhs, rawx); return; }

    ConstMultipleFloat_t b(import(rawx));
    ConstMultipleFloat_t c(import(rawy));

    // Heuristic check for "perfect cancellation"
    if (mpf_sgn(b[0]) == -mpf_sgn(c[0]))
    {
      const bool3 cancel = myIsEqualNZIgnoreSign(rawx, rawy);
      if (IsUncertain3(cancel))
        ThrowInsuffPrec();
      if (IsTrue3(cancel))
      {
        myAssignZero(rawlhs);
        return;
      }
    }
    // Now we actually compute the sum.
    MultipleFloat_t a(import(rawlhs));
    mpf_add(a[0], b[0], c[0]);
    mpf_add(a[1], b[1], c[1]);

    myCheckValidity(rawlhs); // could throw InsufficientPrecision
  }


  void RingTwinFloatImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    // Eliminate the trivial cases of x=0 or y=0; simplifies later code.
    if (myIsZero(rawx)) { myNegate(rawlhs, rawy); return; }
    if (myIsZero(rawy)) { myAssign(rawlhs, rawx); return; }

    ConstMultipleFloat_t b(import(rawx));
    ConstMultipleFloat_t c(import(rawy));

    // Heuristic check for perfect cancellation
    if (mpf_sgn(b[0]) == mpf_sgn(c[0]))
    {
      const bool3 cancel = myIsEqualNZIgnoreSign(rawx, rawy);
      if (IsUncertain3(cancel))
        ThrowInsuffPrec();
      if (IsTrue3(cancel))
      {
        myAssignZero(rawlhs);
        return;
      }
    }

    // Now we actually compute the difference.
    MultipleFloat_t a(import(rawlhs));
    mpf_sub(a[0], b[0], c[0]);
    mpf_sub(a[1], b[1], c[1]);

    myCheckValidity(rawlhs); // could throw InsufficientPrecision
  }


  void RingTwinFloatImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawx) || myIsZero(rawy)) { myAssignZero(rawlhs); return; }
//???     if (myIsOne(rawx)) { myAssign(rawlhs, rawy); return; }
//???     if (myIsOne(rawy)) { myAssign(rawlhs, rawx); return; }
    MultipleFloat_t a(import(rawlhs));
    ConstMultipleFloat_t b(import(rawx));
    ConstMultipleFloat_t c(import(rawy));
    mpf_mul(a[0], b[0], c[0]);
    mpf_mul(a[1], b[1], c[1]);

    myCheckValidity(rawlhs); // could throw InsufficientPrecision
  }


  void RingTwinFloatImpl::myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    CoCoA_ASSERT(!myIsZero(rawy));
    if (myIsZero(rawx)) { myAssignZero(rawlhs); return; }
//???    if (myIsOne(rawy)) { myAssign(rawlhs, rawx); return; }
    MultipleFloat_t a(import(rawlhs));
    ConstMultipleFloat_t b(import(rawx));
    ConstMultipleFloat_t c(import(rawy));
    mpf_div(a[0], b[0], c[0]);
    mpf_div(a[1], b[1], c[1]);

    myCheckValidity(rawlhs); // could throw InsufficientPrecision
  }


  bool RingTwinFloatImpl::myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    if (myIsZero(rawy)) return false;
    myDiv(rawlhs, rawx, rawy); // could throw InsufficientPrecision
    return true;
  }


  bool RingTwinFloatImpl::myIsInvertible(ConstRawPtr rawx) const
  {
    return !myIsZero(rawx);
  }


  void RingTwinFloatImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
    myBinaryPower(rawlhs, rawx, n);
  }


  bool RingTwinFloatImpl::myIsPrintedWithMinus(ConstRawPtr rawx) const
  {
    return mySign(rawx) < 0;
  }


  // The test for myIsInteger is important (see redmine 853).
  // Easy case if myIsInteger throws or returns true.
  // If myIsInteger returns false then result is floor of primary component.
  BigInt RingTwinFloatImpl::myFloor(ConstRawPtr rawx) const
  {
    BigInt N;
    if (myIsInteger(N,rawx)) // IMPORTANT! This will throw if "outer interval" contains more than 1 integer
      return N;

    // We get here only if the outer interval contains no integer.
    ConstMultipleFloat_t X(import(rawx));
    // Next two lines put floor of primary component into N.
    mpz_set_f(mpzref(N), X[0]); // truncates towards 0...
    if (mpf_sgn(X[0]) < 0) //... so must adjust if negative.
      --N;
    return N;
  }

  // The test for myIsInteger is important (see redmine 853).
  // Easy case if myIsInteger throws or returns true.
  // If myIsInteger returns false then result is ceil of primary component.
  BigInt RingTwinFloatImpl::myCeil(ConstRawPtr rawx) const
  {
    BigInt N;
    if (myIsInteger(N,rawx)) // IMPORTANT! This will throw if "outer interval" contains more than 1 integer
      return N;

    // We get here only if the outer interval contains no integer.
    ConstMultipleFloat_t X(import(rawx));
    // Next two lines put ceil of primary component into N.
    mpz_set_f(mpzref(N), X[0]); // truncates towards 0...
    if (mpf_sgn(X[0]) > 0) // ... so must adjust if positive
      ++N;
    return N;
  }

  
  // Separate positive and negative cases to implement easily
  // the strategy that halves round away from zero.
  // Current impl creates a wasteful temporary -- does it really matter?
  BigInt RingTwinFloatImpl::myNearestInt(ConstRawPtr rawx) const
  {
    if (mySign(rawx) >= 0)
      return floor(BigRat(1,2) + RingElemAlias(ring(this), rawx));
    else
      return ceil(BigRat(-1,2) + RingElemAlias(ring(this), rawx));
  }


  // Normal printing function: see myDebugPrint for special debugging print out
  void RingTwinFloatImpl::myOutput(std::ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    using std::ceil;
    // Special case if x happens to be an integer or rational
    try
    {
      BigRat q;
      if (myIsRational(q, rawx)) // could throw
      {
        out << q;
        return;
      }
    }
    catch (const RingTwinFloat::InsufficientPrecision&) {} // discard InsufficientPrecision if it happens
    ConstMultipleFloat_t val(import(rawx));
    CoCoA_ASSERT(mpf_sgn(val[0]) != 0);

    const long EstimatedBitPrecision = -1+LogRelDiff(val[0], val[1], myWorkingBits);
    const long SigFig = static_cast<long>(ceil(EstimatedBitPrecision*0.30103)); // 0.30103 approx log(2)/log(10)

    // Print format is 0.123456789*10^123
    // Print format is -0.123456789*10^(-123)
    // In the line below 10 = strlen("-0.") + strlen("*10^(-") + strlen(")")
    const long NumChars = 10+ SigFig + numeric_limits<long>::digits10;
    const int base = 10;
    vector<char> buffer(NumChars);
    mp_exp_t exp;
    const char* const mantissa = mpf_get_str(&buffer[0], &exp, base, SigFig, val[0]);
    if (mpf_sgn(val[0]) < 0)
      out << "-0." << &mantissa[1]; // NB &mantissa[1] is a trick to skip the "-" in mantissa[0]
    else
      out << "0." << mantissa;
    if (exp > 0) out << "*10^" << exp;
    if (exp < 0) out << "*10^(" << exp << ")";
  }


  void RingTwinFloatImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

//     out << "RingTwinFloat(AccuracyBits=" << myAccuracyBits
//         << ", BufferBits=" << myBufferBits
//         << ", NoiseBits=" << myNoiseBits
//         << ")";
    out << "RingWithID(" << myID 
        << ",\"RingTwinFloat(" << myAccuracyBits
        << ")\")";
  }


  void RingTwinFloatImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("cocoa", "RingTwinFloat");
    OMOut << myAccuracyBits
          << myBufferBits
          << myNoiseBits;
    OMOut->mySendApplyEnd();
  }


  void RingTwinFloatImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "Sending Twin Floats via OpenMath not yet properly defined");
    using std::ceil;
    ConstMultipleFloat_t val(import(rawx));
//?????    if (mpf_sgn(val[0]) == 0) return out << "0";

    const long EstimatedBitPrecision = -1+LogRelDiff(val[0], val[1], myWorkingBits);
    const long SigFig = static_cast<long>(ceil(EstimatedBitPrecision*0.30103)); // 0.30103 approx log(2)/log(10)

    // Print format is 0.123456789*10^123
    // Print format is -0.123456789*10^(-123)
    // In the line below 10 = strlen("-0.") + strlen("*10^(-") + strlen(")")
    const long NumChars = 10+ SigFig + numeric_limits<long>::digits10;
    const int base = 10;
    vector<char> buffer(NumChars);
    mp_exp_t exp;
    mpz_t mantissa;
    mpz_init(mantissa);
    mpf_get_str(&buffer[0], &exp, base, SigFig, val[0]);
    mpz_set_str(mantissa, &buffer[0], 10);
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("bigfloat1","bigfloat");
    OMOut << BigIntFromMPZ(mantissa)
          << long(base)
          << static_cast<long>(exp)-SigFig;  // BUG: slight risk of overflow
    OMOut->mySendApplyEnd();

////   BUG???  Should the secondary component be sent too???
  }


  bool RingTwinFloatImpl::myIsZero(ConstRawPtr rawx) const
  {
    return mpf_sgn(import(rawx)[0]) == 0;
  }


  bool RingTwinFloatImpl::myIsOne(ConstRawPtr rawx) const
  {
    return myIsEqual(rawx, raw(*myOnePtr));
  }


  bool RingTwinFloatImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return myIsEqual(rawx, raw(*myMinusOnePtr));
  }


  bool RingTwinFloatImpl::myIsInteger(BigInt& N, ConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) { N = 0; return true; }
    // General case: check there is exactly 1 integer in outer interval,
    // and that this integer lies in the inner interval.
    ConstMultipleFloat_t x = import(rawx);
    const long prec = myWorkingBits;
    bool3 FinalResult; // default undecided --> ThrowInsuffPrec
    mpf_t semiwidth, upb, lwb, candidate;
    mpf_init2(semiwidth,prec);
    mpf_init2(upb,prec);
    mpf_init2(lwb,prec);
    mpf_init2(candidate,prec);
    myOuterSemiwidth(semiwidth, x);
    mpf_add(upb, x[0], semiwidth);
    mpf_sub(lwb, x[0], semiwidth);
    mpf_floor(candidate, upb);
    if (mpf_cmp(candidate, lwb) <= 0) // no integer in outer interval
    {
      FinalResult = false3;
      goto clean_and_return;
    }
    mpf_add_ui(lwb, lwb, 1);
    if (mpf_cmp(candidate, lwb) > 0) // more than one integer in outer interval
      goto clean_and_return; // will throw InsuffPrec

    // Now check candidate is in inner interval
    // We could simply check that abs(candidate-x[0]) < InnerSemiwidth...
    myInnerSemiwidth(semiwidth, x);
    mpf_add(upb, x[0], semiwidth);
    mpf_sub(lwb, x[0], semiwidth);

    if (mpf_cmp(upb, candidate) >= 0 && mpf_cmp(lwb, candidate) <= 0)
    {
      // Candidate is in inner interval, so we accept it!
      mpz_set_f(mpzref(N), candidate);
      FinalResult = true3;
    }
    clean_and_return:
    mpf_clear(candidate);
    mpf_clear(lwb);
    mpf_clear(upb);
    mpf_clear(semiwidth);
    if (IsUncertain3(FinalResult))
      ThrowInsuffPrec();
    return IsTrue3(FinalResult);
  }


  bool RingTwinFloatImpl::myIsRational(BigInt& N, BigInt& D, ConstRawPtr rawx) const
  {
    // First check whether value is integer (or too close to tell for certain).
    if (myIsInteger(N, rawx)) { D = 1; return true; }  // this may throw (& we want it to)

    ConstMultipleFloat_t x = import(rawx);

// First block below calls SimplestBigRatBetween; the code should be
// simpler and clearer than the code which reimplements the algorithm.
// However tests on my aging MacBook Pro show that this simpler code
// is about 7 times slower than the reimplementation (in the second block).
// See redmine issue #860  (before deleting or altering the code below!)
// See also redmine issue 897 (why is SimplestBigRatBetween so slow?)
#if 0
    long exp;
    mpf_get_d_2exp(&exp, x[0]);
    if (exp < -myWorkingBits) ThrowInsuffPrec();
    mpf_t tmp; mpf_init2(tmp, myWorkingBits+1);
    mpf_t semiwidth; mpf_init2(semiwidth, myWorkingBits);
    BigRat upb, lwb;
    myOuterSemiwidth(semiwidth, x);
    mpf_abs(tmp, x[0]);
    mpf_add(tmp, tmp, semiwidth);
    mpq_set_f(mpqref(upb), tmp);
    mpf_abs(tmp, x[0]);
    mpf_sub(tmp, tmp, semiwidth);
    mpq_set_f(mpqref(lwb), tmp);
    BigRat candidate = SimplestBigRatBetween(lwb, upb);
    if (mpf_sgn(x[0]) < 0)
      mpq_neg(mpqref(candidate), mpqref(candidate));
    mpf_set_q(tmp, mpqref(candidate));
    mpf_sub(tmp, tmp, x[0]);
    mpf_abs(tmp, tmp);
    myInnerSemiwidth(semiwidth, x);
    if (mpf_cmp(tmp, semiwidth) > 0) ThrowInsuffPrec();
    N = num(candidate);
    D = den(candidate);
    return true;
#else
    
    // Now we know that x is not an integer (in ptic x is non-zero), and
    // not huge (since myIsInteger would throw InsuffPrec).
    const long prec = myWorkingBits + 32; // (lazy) extra 32 bits to absorb rounding error.

    // Store the sign of x; the main algm works with positive values, will restore sign at the end.
    const int sign = mpf_sgn(x[0]);

    // This is the main part; put in hi and lo upper and lower bounds for rational value.
    mpf_t lo, hi, tmp; // tmp gets used for various temporary results.
    mpf_init2(lo, prec); mpf_init2(hi, prec); mpf_init2(tmp, prec);

    // Next few lines set lo and hi to lwb and upb of outer(x) made positive.
    myOuterSemiwidth(tmp, x);
    mpf_abs(lo, x[0]);
    mpf_sub(lo, lo, tmp);
    mpf_abs(hi, x[0]);
    mpf_add(hi, hi, tmp);
//       std::cout<<"flt_lo=";mpf_out_str(stdout,10,0,lo);std::cout<<std::endl;
//       std::cout<<"flt_hi=";mpf_out_str(stdout,10,0,hi);std::cout<<std::endl;

    // Next we emulate SimplestRatBetween.
    // Now generate continued fraction expansions of lwb and upb;
    // keep going while the expansions agree.
    mpz_t N0, N1, D0, D1;
    mpz_init_set_ui(N0, 1);    mpz_init_set_ui(N1, 0);    mpz_init_set_ui(D0, 0);    mpz_init_set_ui(D1, 1);
    mpz_t int_lo, int_hi, ztmp;
    mpz_init(int_lo); mpz_init(int_hi); mpz_init(ztmp);

    mpz_set_f(int_lo, lo);  // int_lo = floor(lo)
    mpz_set_f(int_hi, hi);  // int_hi = floor(hi)
    while (mpz_cmp(int_lo, int_hi) == 0)
    {
//       std::cout<<"quot_lo=";mpz_out_str(stdout,10,int_lo);std::cout<<std::endl;
//       std::cout<<"quot_hi=";mpz_out_str(stdout,10,int_hi);std::cout<<std::endl;
//       std::cout<<"flt_lo=";mpf_out_str(stdout,10,0,lo);std::cout<<std::endl;
//       std::cout<<"flt_hi=";mpf_out_str(stdout,10,0,hi);std::cout<<std::endl;
      mpz_mul(ztmp, int_lo, N0);
      mpz_add(N1, N1, ztmp);
      mpz_mul(ztmp, int_lo, D0);
      mpz_add(D1, D1, ztmp);
      mpz_swap(N0, N1);
      mpz_swap(D0, D1);
      if (mpf_integer_p(lo)) break;
      mpf_set_z(tmp, int_lo);
      mpf_sub(lo, lo, tmp);
      mpf_ui_div(lo, 1, lo);
      mpf_sub(hi, hi, tmp);
      mpf_ui_div(hi, 1, hi);
      mpf_swap(lo, hi);
      mpz_set_f(int_lo, lo);
      mpz_set_f(int_hi, hi);
    }
    if (mpz_cmp(int_lo, int_hi) != 0)
    {
      mpf_ceil(lo, lo);
      mpz_set_f(int_lo, lo);
      mpz_mul(ztmp, int_lo, N0);
      mpz_add(N1, N1, ztmp);
      mpz_mul(ztmp, int_lo, D0);
      mpz_add(D1, D1, ztmp);
      mpz_swap(N0, N1);
      mpz_swap(D0, D1);
    }

    // Check that rational is "simple" compared to precision used.
    bool AnswerIsGood = true;
    if (sign < 0) mpz_neg(N0, N0); // make negative if input was negative

    // Next lines compare  |x[0] - N0/D0|  with  myInnerSemiwidth
    // i.e. decide whether  N0/D0  lies inside the inner interval.
    mpf_t tmp2; mpf_init2(tmp2, prec);
    mpf_set_z(tmp, N0);
    mpf_set_z(tmp2, D0);
    mpf_div(tmp, tmp, tmp2);
    mpf_sub(tmp, tmp, x[0]);
    mpf_abs(tmp, tmp);
    myInnerSemiwidth(tmp2, x);
    if (mpf_cmp(tmp, tmp2) > 0)
      AnswerIsGood = false;
    mpf_clear(tmp2);

    if (AnswerIsGood)
    {
      mpz_swap(mpzref(N), N0);
      mpz_swap(mpzref(D), D0);
    }
    mpz_clear(ztmp); mpz_clear(int_hi); mpz_clear(int_lo);
    mpz_clear(D1);      mpz_clear(D0);      mpz_clear(N1);      mpz_clear(N0);
    mpf_clear(tmp); mpf_clear(hi), mpf_clear(lo);
    if (!AnswerIsGood)
      ThrowInsuffPrec();
    return true;
#endif
  }


  bool RingTwinFloatImpl::myIsRational(BigRat& Q, ConstRawPtr rawx) const
  {
    BigInt N,D;
    const bool OK = myIsRational(N,D,rawx);
    if (OK)
    {
      mpz_set(mpq_numref(mpqref(Q)), mpzref(N));
      mpz_set(mpq_denref(mpqref(Q)), mpzref(D));
      // No need to call mpq_canonicalize since N & D are coprime.
    }
    return OK;
  }


  bool RingTwinFloatImpl::myIsDouble(double& d, ConstRawPtr rawx) const
  {
    ConstMultipleFloat_t x = import(rawx);
    long exp;
    d = mpf_get_d_2exp(&exp, x[0]);
    if (numeric_limits<double>::radix != 2) CoCoA_THROW_ERROR(ERR::NYI, "RingTwinFloatImpl::myIsDouble");
    if (exp < numeric_limits<double>::min_exponent) { d=0; return true; }  // ???false also for underflow???
    if (exp > numeric_limits<double>::max_exponent) return false;
    d = ldexp(d,exp);
    return true;
  }


  bool RingTwinFloatImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return myCmp(rawx, rawy) == 0; // myCmp may throw InsufficientPrecision
  }


  int RingTwinFloatImpl::myCmp(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    ConstMultipleFloat_t a(import(rawx));
    ConstMultipleFloat_t b(import(rawy));
    // Deal with obvious cases based just on sign information.
    const int signx = mpf_sgn(a[0]);
    const int signy = mpf_sgn(b[0]);
    if (signy == 0) return signx;
    if (signx == 0) return -signy;
    if (signx != signy) return signx;

    // Now we know that a[0] and b[0] have the same sign.
    const bool3 equal = myIsEqualNZIgnoreSign(rawx, rawy);
    if (IsUncertain3(equal))
      ThrowInsuffPrec();
    if (IsTrue3(equal))
      return 0;
    // Now we know that a[0] and b[0] are definitely unequal.
    if (mpf_cmp(a[0], b[0]) > 0) return 1;
    return -1;
  }


  int RingTwinFloatImpl::mySign(ConstRawPtr rawx) const
  {
    return mpf_sgn(import(rawx)[0]);
  }


  //---------------------------------------------------------------------------

  RingHom RingTwinFloatImpl::myCompose(const RingHom& phi, const RingHom& theta) const
  {
    if (IsExact(codomain(phi)) && IsExact(codomain(theta)))
      return NewApproxHom(domain(theta), codomain(phi));
//!!!PHILOSOPHICAL QUESTIONS TO ANSWER!!!
//NYI!!!    return RingHomComposite(phi,theta);
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingTwinFloatImpl::myCompose -- how did you get here?");
    return phi; // just to keep compiler quiet
  }


  bool RingTwinFloatImpl::myImageLiesInSubfield(const RingHom& phi) const
  {
    (void)(phi); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(codomain(phi) == ring(this));
    return true; // regard RingTwinFloat as a field.
  }


  ideal RingTwinFloatImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    // Is it wise to allow one to create an ideal in RingTwinFloat???
    return NewFieldIdeal(ring(this), gens);
  }


  // Check that the pair of non-zero bigfloat values in rawx is a valid
  // representation for a twin float i.e. that the two components agree
  // to at least myAccuracyBits bits of accuracy.  Will call myPerturb on rawx if
  // the two components agree too closely (more than myWorkingBits-myNoiseBits/2).
  // If rawx is NOT A VALID REPR, then THROWS InsufficientPrecision.
  RingTwinFloatImpl::RelErrRating RingTwinFloatImpl::myCheckRelErr(ConstRawPtr rawx) const
  {
    VerboseLog VERBOSE("TwinFloat::myCheckRelErr");
    if (VerbosityLevel() >= 200)
    {
      VERBOSE(200) <<"myCheckValidity:  INPUT\n";    myDebugPrint(VERBOSE(200),rawx);
    }
    ConstMultipleFloat_t X(import(rawx));
    CoCoA_ASSERT(mpf_sgn(X[0]) != 0);
    CoCoA_ASSERT(mpf_sgn(X[1]) != 0);
    if (mpf_sgn(X[1]) != mpf_sgn(X[0]))
      return TooBig;

    const long accuracy = LogRelDiff(X[0], X[1], myWorkingBits);
    VERBOSE(200) << "acc=" << accuracy << std::endl;
    if (accuracy < myAccuracyBits)
      return TooBig;

    // Add more noise, if rel diff is too small.
    if (accuracy >= myWorkingBits-myNoiseBits/2)
      return TooSmall;
    return OK;
  }

  void RingTwinFloatImpl::myCheckValidity(RawPtr rawx) const
  {
    const RelErrRating outcome = myCheckRelErr(rawx);
    if (outcome == TooBig)
      ThrowInsuffPrec();
    if (outcome == TooSmall)    // Add more noise, if rel diff is too small.
      myPerturb(import(rawx));
  }


  // Produce a valid secondary component of the twin-float x by setting it to
  // primary component plus a uniform random perturbation of relative size at most 1/2^(myAcc+myBuff).
  // We guarantee that the relative difference is AT LEAST 1/2^(myAcc+myBuff+myNoiseBits/2)
  void RingTwinFloatImpl::myPerturb(MultipleFloat_t x) const
  {
    // Generate uniform random "noise" of length myNoiseBits in interval [-1,1],
    // but also insist that noise has magnitude at least 2^(-1/2*myNoiseBits).
    mpf_t noise;
    mpf_init2(noise, myNoiseBits);
    while (true)
    {
      do
        mpf_urandomb(noise, myRandomSource, myNoiseBits);
      while (mpf_sgn(noise) == 0); // JAA: I doubt this makes any practical difference
      // At this point noise is uniform in the interval (0,1)
      mpf_mul_2exp(noise, noise, 1);
      mpf_sub_ui(noise, noise, 1);
      // Now noise is uniform in the interval (-1,1)
      long ExpNoise;
      mpf_get_d_2exp(&ExpNoise, noise);
      const long NumZeroes = -ExpNoise;
      if (NumZeroes < myNoiseBits/2) break;
    }

    // Now set the secondary component to be primary component plus a relative
    // difference of noise/2^(myAccuracyBits+myBufferBits).
    mpf_div_2exp(noise, noise, myAccuracyBits+myBufferBits);
    mpf_mul(noise, noise, x[0]);
    mpf_add(x[1], x[0], noise);
    mpf_clear(noise);
  }


  // Func below is JUST FOR DEBUGGING!!!
  std::ostream& operator<<(std::ostream& out, mpf_t x)
  {
    if (!out) return out;  // short-cut for bad ostreams
//     long exp;
//     const double mant = mpf_get_d_2exp(&exp, x);
//     return out << mant <<"*2^(" << exp << ")";

    using std::ceil;
    const long bits = mpf_get_prec(x);
    const long SigFig = static_cast<long>(ceil(bits*0.30103)); // 0.30103 approx log(2)/log(10)
    // Print format is 0.123456789*10^123
    // Print format is -0.123456789*10^(-123)
    // In the line below 10 = strlen("-0.") + strlen("*10^(-") + strlen(")")
    const long NumChars = 10+ SigFig + numeric_limits<long>::digits10;
    const int base = 10;
    vector<char> buffer(NumChars);
    mp_exp_t exp;
    const char* mantissa = mpf_get_str(&buffer[0], &exp, base, SigFig, x);
    if (mpf_sgn(x) < 0)
      out << "-0." << &mantissa[1]; // NB &mantissa[1] is a trick to skip the "-" in mantissa[0]
    else
      out << "0." << mantissa;
    if (exp > 0) out << "*10^" << exp;
    if (exp < 0) out << "*10^(" << exp << ")";
    return out;
  }


  // ASSUMES x and y are non-zero.
  // IGNORES signs of x and y; ignores the trivial case of exponents differing by more than 1.
  // MAY THROW InsufficientPrecision!
  bool3 RingTwinFloatImpl::myIsEqualNZIgnoreSign(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    /***********************************/
    /* WARNING: uses flag myEqTestMode */
    /* See also redmine issue 859      */
    /***********************************/
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsZero(rawy));

    ConstMultipleFloat_t x(import(rawx)); // an alias
    ConstMultipleFloat_t y(import(rawy)); // an alias

    bool3 outcome; // initially "uncertain3"

    mpf_t AbsX; mpf_init2(AbsX, myWorkingBits);
    mpf_t AbsY; mpf_init2(AbsY, myWorkingBits);
    mpf_abs(AbsX, x[0]);
    mpf_abs(AbsY, y[0]);
    const int CmpXY = mpf_cmp(AbsX, AbsY);
    if (CmpXY == 0) { outcome = true3; goto TidyUp2; }

    mpf_t XOuterSW; mpf_init2(XOuterSW, myWorkingBits);
    mpf_t XInnerSW; mpf_init2(XInnerSW, myWorkingBits);
    myOuterSemiwidth(XOuterSW, x);
    myInnerSemiwidth(XInnerSW, x);

    mpf_t YOuterSW; mpf_init2(YOuterSW, myWorkingBits);
    mpf_t YInnerSW; mpf_init2(YInnerSW, myWorkingBits);
    myOuterSemiwidth(YOuterSW, y);
    myInnerSemiwidth(YInnerSW, y);

    // By swapping ensure that AbsX > AbsY
    if (CmpXY < 0)
    {
      mpf_swap(AbsX, AbsY);
      mpf_swap(XOuterSW, YOuterSW);
      mpf_swap(XInnerSW, YInnerSW);
    }

    mpf_t XYdist; mpf_init2(XYdist, myWorkingBits);
    mpf_sub(XYdist, AbsX, AbsY); // distance between AbsX and AbsY

    mpf_t tmp; mpf_init2(tmp, myWorkingBits);

    // If AbsX-AbsY > outersw(x)+outersw(y) the surely the outer intervals
    // are disjoint, so X and Y are surely unequal.
    mpf_add(tmp, XOuterSW, YOuterSW);
    if (mpf_cmp(XYdist, tmp) >= 0) { outcome = false3; goto TidyUp; }

    // The outer intervals do meet, so check the inner intervals.
    // Just check equiv condition AbsX-AbsY < innersw(x)+innersw(y)
    mpf_add(tmp, XInnerSW, YInnerSW);
    if (mpf_cmp(XYdist, tmp) < 0) { outcome = true3; goto TidyUp; }

    // Now, if requested, check whether permitted amount of overlap occurs
    if (myEqTestMode == AllowOverlap)
    {
      if (mpf_cmp(XYdist, XOuterSW) > 0 && mpf_cmp(XYdist, YOuterSW) > 0) outcome = false3;
    }
      
    TidyUp:
    mpf_clear(tmp);
    mpf_clear(XYdist);
    mpf_clear(YInnerSW);
    mpf_clear(YOuterSW);
    mpf_clear(XInnerSW);
    mpf_clear(XOuterSW);
    TidyUp2:
    mpf_clear(AbsY);
    mpf_clear(AbsX);

    return outcome;
  }


  RingTwinFloat::InsufficientPrecision::InsufficientPrecision():
      ErrorInfo(ERR::InsuffPrec, "RingTwinFloat arithmetic")
  {}


  RingTwinFloatImpl::HomExactImpl::HomExactImpl(const RingTwinFloat& domain, const ring& codomain):
      RingHomBase(domain, codomain)
  {
    CoCoA_ASSERT(IsExact(codomain));
  }


  void RingTwinFloatImpl::HomExactImpl::myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const
  {
    BigRat tmp;
    if (!myDomain->myIsRational(tmp, arg)) CoCoA_THROW_ERROR(ERR::BadArg, "Applying RingTwinFloat exact hom");
    myCodomain->myAssign(image, tmp);
  }


  RingTwinFloatImpl::HomApproxImpl::HomApproxImpl(const RingTwinFloat& domain, const ring& codomain):
      RingHomBase(domain, codomain)
  {
    CoCoA_ASSERT(!IsExact(codomain));
    CoCoA_ASSERT(IsRingTwinFloat(codomain));
  }


  void RingTwinFloatImpl::HomApproxImpl::myApply(RingElemRawPtr image, RingElemConstRawPtr arg) const
  {
    if (myDomain->myIsZero(arg)) { myCodomain->myAssignZero(image); return; }
    RingTwinFloatPtr(myCodomain)->myRecvTwinFloat(image,arg);  // CAREFUL HERE image & arg (probably) belong to different RingTwinFloats!!!
  }


  //---------------------------------------------------------------------------

  bool IsPracticallyEqual(ConstRefRingElem x, ConstRefRingElem y)
  {
    if (owner(x) != owner(y))
      CoCoA_THROW_ERROR(ERR::MixedRings, "IsPracticallyEqual(RingElem,RingElem)");
    if (!IsRingTwinFloat(owner(x)))
      CoCoA_THROW_ERROR(ERR::NotRingTwinFloat, "IsPracticallyEqual(RingElem,RingElem)");

    try
    {
      return x == y;
    }
    catch (const RingTwinFloat::InsufficientPrecision&)
    {
      return false;
    }
  }


  RingHom NewApproxHom(const ring& RR, const ring& S)
  {
    if (!IsRingTwinFloat(RR))
      CoCoA_THROW_ERROR(ERR::NotRingTwinFloat, "NewApproxHom (1st arg)");
    if (IsExact(S))
      return RingHom(new RingTwinFloatImpl::HomExactImpl(RR, S));
    if (IsRingTwinFloat(S))
      return RingHom(new RingTwinFloatImpl::HomApproxImpl(RR, S));
    CoCoA_THROW_ERROR(ERR::NYI, "NewApproxHom unhandled case");
    return IdentityHom(S); // just to keep compiler quiet
  }


  void RingTwinFloatImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams

    using std::endl;
    using std::ceil;
    ConstMultipleFloat_t val(import(rawx));
    if (mpf_sgn(val[0]) == 0) { out << "0"; return; } // BUG??? Can this condition ever be true???

    const long SigFig = static_cast<long>(2+ceil(myWorkingBits*0.30103)); // 0.30103 approx log(2)/log(10)
    // Print format is 0.123456789*10^123
    // Print format is -0.123456789*10^(-123)
    // In the line below 10 = strlen("-0.") + strlen("*10^(-") + strlen(")")
    const long NumChars = 10+ SigFig + numeric_limits<long>::digits10;
    const int base = 10;
    vector<char> buffer(NumChars);
    for (int component=0; component <= 1; ++component)
    {
      out << "TwinFloat component " << component << ": ";
      mp_exp_t exp;
      const char* const mantissa = mpf_get_str(&buffer[0], &exp, base, SigFig, val[component]);
      if (mpf_sgn(val[component]) < 0)
        out << "-0." << &mantissa[1]; // NB &mantissa[1] is a trick to skip the "-" in mantissa[0]
      else
        out << "0." << mantissa;
      if (exp > 0) out << "*10^" << exp;
      if (exp < 0) out << "*10^(" << exp << ")";
      out << std::endl;
    }

    mpf_t tmp;
    mpf_init2(tmp,64);
    myInnerSemiwidth(tmp, import(rawx));
    const double InnerSemiwidth = mpf_get_d(tmp);
    myOuterSemiwidth(tmp, import(rawx));
//    mpf_sub(tmp, val[0], val[1]);
//    mpf_abs(tmp, tmp);
    const double OuterSemiwidth = mpf_get_d(tmp);
    out << "InnerSemiwidth: " << InnerSemiwidth << endl;
    out << "OuterSemiwidth: " << OuterSemiwidth << std::endl;
    const long EstimatedBitPrecision = -1+LogRelDiff(val[0], val[1], myWorkingBits);
    out << "EstBitPrec: " << EstimatedBitPrecision << endl;
    mpf_reldiff(tmp, val[0], val[1]);
    mpf_abs(tmp, tmp);
    const double RelDiff = mpf_get_d(tmp);
    out << "Rel diff: " << RelDiff << endl;
    out << "log2(..): " << std::log(RelDiff)/std::log(2.0) << endl;
    mpf_clear(tmp);
  }

  void DebugPrint(std::ostream& out, ConstRefRingElem x)
  {
    if (!IsRingTwinFloat(owner(x)))
      CoCoA_THROW_ERROR("Only for elems of RingTwinFloat", "DebugPrint");
    dynamic_cast<const RingTwinFloatImpl*>(owner(x).myRawPtr())->myDebugPrint(out, raw(x));
  }


  RingTwinFloat::RingTwinFloat(const ring& R):
      ring(RingTwinFloatPtr(R, CoCoA_ERROR_CONTEXT))
  {}


  RingTwinFloat::RingTwinFloat(const RingTwinFloatBase* RingPtr):
      ring(RingPtr)
  {}


  RingTwinFloat NewRingTwinFloat(const MachineInt& AccuracyBits)
  {
    if (IsNegative(AccuracyBits) || ! IsSignedLong(AccuracyBits))
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "NewRingTwinFloat(A)");
    const long A = AsSignedLong(AccuracyBits);
    return RingTwinFloat(new RingTwinFloatImpl(A, A, max(32l, A/4)));
  }


  RingTwinFloat NewRingTwinFloat(const MachineInt& AccuracyBits, const MachineInt& BufferBits, const MachineInt& NoiseBits)
  {
    if (IsNegative(AccuracyBits) || !IsSignedLong(AccuracyBits) ||
        IsNegative(BufferBits) || !IsSignedLong(BufferBits) ||
        IsNegative(NoiseBits) || !IsSignedLong(NoiseBits))
      CoCoA_THROW_ERROR(ERR::BadArg, "NewRingTwinFloat(A,B,N): args must be non negative");
    const long A = AsSignedLong(AccuracyBits);
    const long B = AsSignedLong(BufferBits);
    const long N = AsSignedLong(NoiseBits);
    return RingTwinFloat(new RingTwinFloatImpl(A, B, N));
  }


  long RingTwinFloatImpl::myPrecisionBits() const
  {
    return myAccuracyBits;
  }


  long PrecisionBits(const RingTwinFloat& RR)
  {
    return RR->myPrecisionBits();
  }


  MantExp2 RingTwinFloatImpl::myExport(RingElemConstRawPtr rawx) const
  {
    if (myIsZero(rawx)) return MantExp2(0,0, BigInt(0), 0);
    
    VerboseLog VERBOSE("TwinFloat::myExport");
    ConstMultipleFloat_t X(import(rawx));
    const long MaxPrec = myAccuracyBits+myBufferBits;
    const long CurrPrec = LogRelDiff(X[0],X[1],MaxPrec); // 3rd arg should never be used
    VERBOSE(200) << "CurrPrec=" << CurrPrec << std::endl;
    long exp;
    mpf_get_d_2exp(&exp, X[0]);
    VERBOSE(200) << "exp=" << exp << std::endl;
    mpf_t tmp; mpf_init2(tmp, CurrPrec+32);
    const long shift = CurrPrec+1-exp;  // extra +1 to facilitate rounding
    if (shift >= 0)
      mpf_mul_2exp(tmp, X[0], shift);
    else
      mpf_div_2exp(tmp, X[0], -shift);

    int s=1;
    if (mpf_sgn(tmp) == -1) { s = -1; mpf_neg(tmp, tmp); }
    mpf_add_ui(tmp, tmp, 1u);  // these 3 lines should round the value
    mpf_div_2exp(tmp, tmp, 1); //
    mpf_floor(tmp, tmp);       //
    BigInt mant;
    mpz_set_f(mpzref(mant), tmp);
    VERBOSE(200) << "mant=" << mant << std::endl;
    mpf_clear(tmp);
    if (FloorLog2(mant) == CurrPrec) // check for mantissa overflow
    { mant /= 2; ++exp; }
    return MantExp2(s, exp-1, mant, CurrPrec);
  }

  MantExp2 MantissaAndExponent2(const RingElem& x)
  {
    if (!IsRingTwinFloat(owner(x))) CoCoA_THROW_ERROR(ERR::NotRingTwinFloat, "MantissaAndExponent2");
    if (IsZero(x)) return MantExp2(0,0,BigInt(0),0);
    return RingTwinFloatPtr(owner(x))->myExport(raw(x));
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingTwinFloat.C,v 1.67 2024/07/22 15:43:58 abbott Exp $
// $Log: RingTwinFloat.C,v $
// Revision 1.67  2024/07/22 15:43:58  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.66  2024/07/02 15:36:14  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.65  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.64  2023/03/08 19:43:57  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.63  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.62  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.61  2022/02/05 20:01:33  abbott
// Summary: Minor changes to avoid warnings from clang (redmine 1528)
//
// Revision 1.60  2021/10/29 19:47:59  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.59  2021/01/31 10:03:17  abbott
// Summary: Minor change (after revising bool3), nothing important.
//
// Revision 1.58  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.57  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.56  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.55  2020/02/11 16:12:19  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.54  2020/01/26 14:41:59  abbott
// Summary: Revised includes after splitting NumTheory (redmine 1161)
//
// Revision 1.53  2019/11/14 17:53:04  abbott
// Summary: Removed some commented out cruft; replaced use of clog but VERBOSE
//
// Revision 1.52  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.51  2019/03/04 11:03:40  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.50  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.49  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.48  2018/04/20 18:51:25  abbott
// Summary: Changed ctors for BigInt/BigRat from string or from MPZ/MPQ
//
// Revision 1.47  2017/09/06 11:56:29  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.46  2016/10/12 18:29:43  abbott
// Summary: Fixed memory leak in myIsInteger
//
// Revision 1.45  2016/08/02 12:50:30  abbott
// Summary: More detail in DebugPrint
//
// Revision 1.44  2016/06/29 10:55:57  abbott
// Summary: Minor update after redmine 859
//
// Revision 1.43  2016/04/07 14:42:48  abbott
// Summary: Cleaned up myIsRational
//
// Revision 1.42  2016/03/30 16:27:24  abbott
// Summary: Removed cruft from myNearestInt
//
// Revision 1.41  2016/03/30 16:26:18  abbott
// Summary: Improved comments
//
// Revision 1.40  2016/03/30 13:10:24  abbott
// Summary: Extensive cleaning (esp. myIsEqualNZIgnoreSign, myFloor, myCeil, myNearestInt)
//
// Revision 1.39  2016/03/30 09:47:21  abbott
// Summary: Changed myIsEqualNZIgnoreSign to return a bool3 (and never throw); many consequential changes
//
// Revision 1.38  2016/03/25 20:09:56  abbott
// Summary: Minor cleaning; improved some comments.
//
// Revision 1.37  2016/03/25 20:06:49  abbott
// Summary: Introduced new mem fns myCeil and myNearestInt; impls in ZZ, QQ, & TwinFloat. Renamed NearestInteger to NearestInt.
//
// Revision 1.36  2014/08/16 13:25:17  abbott
// Summary: Changed printing: if IsRational(x) then it prints as a rational
// Author: JAA
//
// Revision 1.35  2014/07/30 14:28:27  abbott
// Summary: Added myExport for MantExp2
// Author: JAA
//
// Revision 1.34  2014/07/09 11:45:08  abbott
// Summary: Removed AsRingTwinFloat
// Author: JAA
//
// Revision 1.33  2014/07/08 13:14:41  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.32  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.31  2014/07/02 16:44:47  bigatti
// -- new way of printing rings with ID
//
// Revision 1.30  2014/06/17 10:13:47  abbott
// Summary: Added (void)(phi) to avoid compiler warning about unused param
// Author: JAA
//
// Revision 1.29  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.28  2014/04/02 10:57:46  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.27  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.26  2013/02/21 14:14:41  abbott
// First attempt at implementing PartialRingHom -- some problems remain!!
//
// Revision 1.25  2012/10/24 13:37:29  abbott
// Replaced ConstRefRingElem by RingElemAlias in local variable type.
//
// Revision 1.24  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.23  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.22  2012/04/27 15:05:48  abbott
// Added mem fn IamFiniteField
//
// Revision 1.21  2011/12/23 14:56:21  bigatti
// -- changed log(2) --> log(2.0)
//
// Revision 1.20  2011/11/09 14:11:58  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.19  2011/08/24 10:29:55  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.18  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.17  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.16  2011/05/19 14:46:26  abbott
// Added defn of myIsDouble.
//
// Revision 1.15  2011/03/14 10:31:16  abbott
// Changed size_t into long (in fn interfaces).
//
// Revision 1.14  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.13  2011/01/19 16:35:26  bigatti
// -- changed ERR::BadArg into ERR::NotNonNegative
//
// Revision 1.12  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.11  2010/02/03 22:32:35  abbott
// Replaced use of std::log2 by std::log (to avoid trouble on cygwin).
//
// Revision 1.10  2010/02/01 22:43:10  abbott
// Changed arg to max so that it compiles on all platforms.
// Changed commented out debugging print cmds so that
// they print on clog.
//
// Revision 1.9  2010/01/20 15:55:27  abbott
// Modified limits for the values of the parameters to the RingTwinFloatImpl ctor.
//
// Revision 1.8  2010/01/19 17:39:32  abbott
// MAJOR OVERHAUL!  Code is now in agreement with the (upcoming) article
// about Twin Floats (see the paper for details).  It is also much more correct.
//
// Revision 1.7  2009/10/26 15:42:29  bigatti
// -- added CopyFromMPZ in ZZ ctor
//
// Revision 1.6  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class QQ, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.5  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.4  2008/11/18 14:32:33  bigatti
// -- added std:: for floor and ceil
//
// Revision 1.3  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/09/24 14:24:41  abbott
// Added parentheses to make code more readable (and to shut up gcc-4.3).
//
// Revision 1.1  2007/06/21 21:29:47  abbott
// Changed name of RingFloat into RingTwinFloat.
//
// Revision 1.3  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.2  2007/03/21 14:52:49  bigatti
// -- added AsRingFloat(R)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.15  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.14  2007/03/07 14:07:39  bigatti
// -- minor: commented argument names for -Wextra
//
// Revision 1.13  2007/01/17 12:32:39  cocoa
// Changed a few more "raw" variable names so that the code compiles fine
// also when CoCoA_DEBUG is set.
//
// Revision 1.12  2007/01/15 14:31:35  bigatti
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.11  2007/01/09 15:52:08  cocoa
// Changed QBGenerator to use std::vector instead of std::list for the result.
// Minor mod to configure script.
//
// Revision 1.10  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.9  2006/12/06 17:36:58  cocoa
// -- commented out 2 lines to keep compiler quiet
//
// Revision 1.8  2006/11/27 14:25:53  cocoa
// -- reorganised #include files
//
// Revision 1.7  2006/11/27 13:06:22  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.6  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
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
// Revision 1.2  2006/07/17 19:40:53  cocoa
// Extensive changes following a slight variation in the semantics.  Added fn IsRational.
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
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.3  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.2  2006/02/20 22:41:19  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
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
// Revision 1.13  2004/12/16 13:34:12  cocoa
// Fixed a memory leak, and noted some areas needing improvement.
//
// Revision 1.12  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.11  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.10  2004/11/05 15:34:33  cocoa
// Consequential change following from the renaming of
// FieldIdealImpl and the introduction of the new pseudo-ctor.
//
// Revision 1.9  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.8  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.7  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.6  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.5  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.4  2004/04/02 16:10:56  cocoa
// -- output: symbol "e" changed into "10^" to be read by cocoa-4
//
// Revision 1.3  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.2  2004/02/03 16:16:20  cocoa
// Removed pointless IamGCDDomain functions from several concrete rings.
// Added IamOrderedDomain functions where appropriate.
// Tidied ctors for the small finite fields.
//
// Revision 1.1  2004/01/28 15:54:09  cocoa
// Sundry additions.
//
