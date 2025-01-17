//   Copyright (c)  2007,2010  Anna Bigatti

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


// Source code for RingDenseUPolyCleanImpl

#include "CoCoA/RingDenseUPolyClean.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseUPolyClean.H"
#include "CoCoA/DenseUPolyRing.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/VectorOps.H" // only for debugging
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"

#include <iostream>
using std::ostream;
// #include <memory>   // from MemPool.H
using std::unique_ptr;
#include <vector>
using std::vector;


namespace CoCoA
{
  class RingDenseUPolyCleanImpl: public DenseUPolyRingBase
  {
  private:
    typedef DenseUPolyClean value_t; // DenseUPolyClean is the actual type of the values in a RingDenseUPolyCleanImpl
    static value_t& import(RingElemRawPtr rawf);
    static const value_t& import(RingElemConstRawPtr rawf);

  public:
    RingDenseUPolyCleanImpl(const ring& R, const symbol& x, long MinCapacity);
    ~RingDenseUPolyCleanImpl() {}


  private: // Data members of RingDenseUPolyCleanImpl
    const ring myCoeffRingValue;  ///< the coefficient ring
    mutable MemPool myDUPPool; ///< memory manager for polynomials
    std::unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    std::unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    std::vector<RingElem> myIndetVector; ///< Vector for compatibility with SparsePolyRing ???
    symbol myIndetSymbolValue; ///< the indet name
    long myMinCapacity;  // the minimum capacity for all coeff vectors

  public:  // functions which every ring must implement
    BigInt myCharacteristic() const override  { return characteristic(myCoeffRingValue); }
    bool IamCommutative() const override  { return true; /*assume CoeffRing comm*/ }
    bool IamIntegralDomain() const /*override*/  { return IsIntegralDomain(myCoeffRingValue); }
    bool IamTrueGCDDomain() const override  { return IsTrueGCDDomain(myCoeffRingValue) || IsField(myCoeffRingValue); }
    bool IamField() const override  { return false; /*??? (myNumIndetsValue==0 && IsField(myCoeffRingValue)) */}
    bool IamFiniteField() const override  { return false; /*??? (myNumIndetsValue==0 && IsFiniteField(myCoeffRingValue)) */}
    ConstRefRingElem myZero() const override  { return *myZeroPtr; }
    ConstRefRingElem myOne() const override  { return *myOnePtr; }
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& N) const override;
    RingElemRawPtr myNew(ConstRawPtr rawcopy) const override;
    void myDelete(RawPtr rawx) const override;                             // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;                  // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;         // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;      // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;          // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                       // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;         // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; // lhs = x-y
    //    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const; // lhs = x^n, n>1, x not -1,0,1
    std::string myImplDetails() const override  {return "RingDenseUPolyClean";}
    //    bool myIsZeroAddMul: use default definition

    // functions which every DenseUPolyRing must implement
    void myAddMulLM(RawPtr rawf, ConstRawPtr rawc, long d, ConstRawPtr rawg) const override;

    // functions which every PolyRing must implement
    const ring& myCoeffRing() const override  { return myCoeffRingValue; }
    const std::vector<RingElem>& myIndets() const override  { return myIndetVector; }

    ///@name Simple functions on polynomials
    //@{
    void myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< WEAK EXCEPTION GUARANTEE
    bool myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< WEAK EXCEPTION GUARANTEE
    void myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const override; ///< lhs = deriv(f, x)
    //@}
    RingHom myCoeffEmbeddingHomCtor() const override;

    //----------------------------------------------------------------------
    // Functions which every DenseUPolyRing must implement:
    //----------------------------------------------------------------------

    ///@name Member functions every concrete DenseUPolyRing implementation must have in addition to those of PolyRingBase.
    //@{
    // ANNA ??? some of these should be in PolyRing
    const symbol& myIndetSymbol() const override  { return myIndetSymbolValue; }
    long myDegPlus1(ConstRawPtr rawf) const override;          ///<  standard degree+1 of f, 0 for zero poly
    long mySize(ConstRawPtr rawf) const override;              ///<  size of f
    void myMulByXExp(RawPtr rawf, unsigned long n) const override;
    void myMulBy1MinusXExp(RawPtr rawf, unsigned long n) const override;
    void myResize(RawPtr rawf, long NewSize) const override;
    void myResetDeg(RawPtr rawf) const override; ///< reset the correct value of deg (assumed to be less than the current value)
    void myAssignZeroCoeff(RawPtr rawf, long d) const override; ///< f_d = 0, no check on size nor degree
    void myAssignNonZeroCoeff(RawPtr rawf, ConstRawPtr rawc, long d) const override; ///< f_d = c, no check on size nor degree
    RingElemAlias myCoeff(ConstRawPtr rawf, long d) const override;
    //@}

    ///@name   Functions for creating/building polynomials
    //@{
    //@}

  private: // homomorphism class
    class HomImpl: public RingHomBase
    {
    public:
      HomImpl(const ring& codomain);
      void myApply(RawPtr rawimage, ConstRawPtr rawarg) const override;
    };

  };


  //----------------------------------------------------------------------

  inline RingDenseUPolyCleanImpl::value_t& RingDenseUPolyCleanImpl::import(RingElemRawPtr rawf)
  {
    return *static_cast<value_t*>(rawf.myRawPtr());
  }
  
  inline const RingDenseUPolyCleanImpl::value_t& RingDenseUPolyCleanImpl::import(RingElemConstRawPtr rawf)
  {
    return *static_cast<const value_t*>(rawf.myRawPtr());
  }
  

  RingDenseUPolyCleanImpl::RingDenseUPolyCleanImpl(const ring& R, const symbol& x, long MinCapacity):
    myCoeffRingValue(R),
    myDUPPool(sizeof(DenseUPolyClean), "RingDenseUPolyCleanImpl::myDUPPool"),
    myIndetSymbolValue(x)
  {
    CoCoA_ASSERT(IsCommutative(myCoeffRingValue));
    myMinCapacity = MinCapacity;
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myIndetVector.resize(1, *myZeroPtr);
    import(raw(myIndetVector[0])).myResize(2);  // deg+1
    import(raw(myIndetVector[0])).myAssignNonZeroCoeff(one(R), 1);
    myRefCountZero(); // otherwise it is 2 + NumIndets and won't be destroyed
  }


  //---- RingDenseUPolyClean specific functions with RingElem

  RingElemRawPtr RingDenseUPolyCleanImpl::myNew() const
  {
    void* ptr = myDUPPool.alloc();
    new(ptr) DenseUPolyClean(myCoeffRingValue, myMinCapacity); // placement new
    return RingElemRawPtr(ptr);
  }


  RingElemRawPtr RingDenseUPolyCleanImpl::myNew(const MachineInt& n) const
  {
    if (IsZero(n)) return myNew();  // not really necessary
    unique_ptr<DenseUPolyClean> ans(new(myDUPPool.alloc()) DenseUPolyClean(myCoeffRingValue, myMinCapacity)); // placement new
    *ans = n;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDenseUPolyCleanImpl::myNew(const BigInt& N) const
  {
    if (N == 0) return myNew();  // not really necessary
    unique_ptr<DenseUPolyClean> ans(new(myDUPPool.alloc()) DenseUPolyClean(myCoeffRingValue, myMinCapacity)); // placement new
    *ans = N;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDenseUPolyCleanImpl::myNew(ConstRawPtr rawcopy) const
  {
    unique_ptr<DenseUPolyClean> ans(new(myDUPPool.alloc()) DenseUPolyClean(import(rawcopy), myMinCapacity)); // placement new
    //    *ans = import(rawcopy);
    return RingElemRawPtr(ans.release());
  }


  void RingDenseUPolyCleanImpl::myDelete(RawPtr rawx) const
  {
    import(rawx).~DenseUPolyClean();
    myDUPPool.free(rawx.myRawPtr());
  }


  void RingDenseUPolyCleanImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    swap(import(rawx), import(rawy));
  }


  void RingDenseUPolyCleanImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingDenseUPolyCleanImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = n;
  }


  void RingDenseUPolyCleanImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = N;
  }


  void RingDenseUPolyCleanImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs).myAssignZero();
  }


  void RingDenseUPolyCleanImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myCoeffRingValue));
    RingElem tmp(myCoeffRingValue);
    myCoeffRingValue->myRecvTwinFloat(raw(tmp), rawx);
    myCoeffEmbeddingHomCtor()->myApply(rawlhs, raw(tmp));
  }


  void RingDenseUPolyCleanImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
//     rawlhs.clear();
//     transform(rawx.begin(), rawx.end(), back_inserter(rawlhs), myNegate());
    import(rawlhs) = import(rawx);
    import(rawlhs).myNegate();
  }


  void RingDenseUPolyCleanImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    add(import(rawlhs), import(rawx), import(rawy));
  }


  void RingDenseUPolyCleanImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    sub(import(rawlhs), import(rawx), import(rawy));
  }


//   void RingDenseUPolyCleanImpl::myDiv(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr rawy) const
//   {
//     CoCoA_ASSERT(!myIsZero(rawy));
//     CoCoA_THROW_ERROR(ERR::NYI, "RingDenseUPolyCleanImpl::myDiv");
//   }


//   bool RingDenseUPolyCleanImpl::myIsDivisible(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr rawy) const
//   {
//     CoCoA_ASSERT(!myIsZero(rawy));
//     CoCoA_THROW_ERROR(ERR::NYI, "RingDenseUPolyCleanImpl::myDiv");
//     return true;
//   }


//   void RingDenseUPolyCleanImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
//   {
//     // Assert that we have a genuinely non-trivial case.
//     CoCoA_ASSERT(n > 1);
//     CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));
//     long NoUse;
//     if (myIsIndet(NoUse, rawx))
//       myIndetPower(rawlhs, 0, n);
//     else
//       CoCoA_THROW_ERROR(ERR::NYI, "RingDenseUPolyCleanImpl::myPowerSmallExp");
//   }


  void RingDenseUPolyCleanImpl::myResize(RawPtr rawf, long NewSize) const
  {
    CoCoA_ASSERT(NewSize < 1000000000);  // to catch silly mistakes
    if (NewSize > myDegPlus1(rawf))  // ??? should this check be here?
    //  (NewSize > mySize())
      import(rawf).myResize(NewSize);
  }


  void RingDenseUPolyCleanImpl::myResetDeg(RawPtr rawf) const
  {
    import(rawf).myResetDeg();
  }


  long RingDenseUPolyCleanImpl::mySize(ConstRawPtr rawf) const
  {
    return import(rawf).mySize();
  }


  long RingDenseUPolyCleanImpl::myDegPlus1(ConstRawPtr rawf) const
  {
    return import(rawf).myDegPlus1();
  }


  void RingDenseUPolyCleanImpl::myAssignZeroCoeff(RawPtr rawf, long d) const
  {
    import(rawf).myAssignZeroCoeff(d);
  }


  void RingDenseUPolyCleanImpl::myAssignNonZeroCoeff(RawPtr rawf, ConstRawPtr rawc, long d) const
  {
    import(rawf).myAssignNonZeroCoeff(RingElemAlias(myCoeffRingValue,rawc), d);
  }


  RingElemAlias RingDenseUPolyCleanImpl::myCoeff(ConstRawPtr rawf, long d) const
  {
    return import(rawf).myCoeff(d);
  }


  //----------------------------------------------------------------------
  // Functions which every PolyRing must implement
  //----------------------------------------------------------------------


  void RingDenseUPolyCleanImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawc, long d, ConstRawPtr rawg) const
  {   // not exception clean  f += c*indet^d*g
    import(rawf).myAddMulLM(RingElemAlias(myCoeffRing(),rawc), d, import(rawg));
  }


  void RingDenseUPolyCleanImpl::myMulByXExp(RawPtr rawf, unsigned long n) const
  {
    import(rawf).myMulByXExp(n);
  }


  void RingDenseUPolyCleanImpl::myMulBy1MinusXExp(RawPtr rawf, unsigned long n) const
  {
    import(rawf).myMulBy1MinusXExp(n);
  }


  bool RingDenseUPolyCleanImpl::myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const
  {
    import(rawf).myDivByCoeff(RingElemAlias(myCoeffRing(),rawc));
    return true;
  }
  
  
  void RingDenseUPolyCleanImpl::myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const
  {
    import(rawf).myMulByCoeff(RingElemAlias(myCoeffRing(),rawc));
  }


  void RingDenseUPolyCleanImpl::myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const
  {
    (void)(rawx); // to avoid compiler warning about unused parameter
    CoCoA_ASSERT(IsIndet(RingElemAlias(ring(this), rawx)));
    deriv(import(rawlhs), import(rawf));
  }


  RingHom RingDenseUPolyCleanImpl::myCoeffEmbeddingHomCtor() const
  {
    return RingHom(new CoeffEmbeddingHomImpl(DenseUPolyRing(this)));
  }


  //----------------------------------------------------------------------
  // Pseudo-ctors for polynomial rings.

  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing, const symbol& IndetName, long MinCapacity)
  {
    if (!IsCommutative(CoeffRing))
      CoCoA_THROW_ERROR(ERR::NotCommutative, "NewPolyRing_DUP(R, x) pseudo ctor");
    if (!IsGoodIndetName(CoeffRing, IndetName))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPolyRing_DUP(R, x) pseudo ctor");
    if (MinCapacity<=0) // or too big?
      CoCoA_THROW_ERROR(ERR::ReqNonNegative, "NewPolyRing_DUP(R, x) pseudo ctor");
    return DenseUPolyRing(new RingDenseUPolyCleanImpl(CoeffRing, IndetName, MinCapacity));
  }

  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing, const symbol& IndetName)
  {
    return NewPolyRing_DUP(CoeffRing, IndetName, 10);
  }

  DenseUPolyRing NewPolyRing_DUP(const ring& CoeffRing)
  {
    return NewPolyRing_DUP(CoeffRing, symbol("x"), 10);
  }
  

} // end of namespace CoCoA
