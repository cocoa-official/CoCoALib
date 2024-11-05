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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingDenseUPolyClean.C,v 1.46 2024/07/02 15:36:14 bigatti Exp $
// $Log: RingDenseUPolyClean.C,v $
// Revision 1.46  2024/07/02 15:36:14  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.45  2023/03/08 19:43:56  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.44  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.43  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.42  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.41  2019/03/04 10:36:13  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.40  2018/05/18 12:15:03  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.39  2018/05/17 15:40:44  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.38  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.37  2014/07/31 14:45:18  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.36  2014/07/28 15:46:56  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.35  2014/07/11 15:36:12  bigatti
// -- removed myOutputSelf (default impl) and added myImplDetails()
//
// Revision 1.34  2014/07/07 16:49:21  bigatti
// -- updated RingWithID
//
// Revision 1.33  2014/07/07 12:47:09  abbott
// Summary: Added (void)(rawx) to avoid a compiler warning
// Author: JAA
//
// Revision 1.32  2014/07/02 16:53:16  bigatti
// -- new way of printing ring with ID
//
// Revision 1.31  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.30  2014/04/30 16:11:11  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.29  2014/01/28 10:58:15  bigatti
// -- just a comment
//
// Revision 1.28  2013/05/22 13:23:29  bigatti
// -- fixed myDeriv, myRecvTwinFloat
//
// Revision 1.27  2013/05/21 16:34:57  abbott
// Added cheap hack to avoid compiler warning when debugging is off.
//
// Revision 1.26  2013/05/20 16:21:19  abbott
// Completed impl of myDeriv (NB work is done by DenseUPolyClean).
//
// Revision 1.25  2012/10/24 12:18:05  abbott
// Changed return type of myCoeff.
// Replaced ConstRefRingElem by RingElemAlias in ctor calls.
//
// Revision 1.24  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.23  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.22  2012/04/27 15:06:09  abbott
// Added mem fn IamFiniteField
//
// Revision 1.21  2011/11/09 14:10:49  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.20  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.19  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.18  2011/04/27 08:21:48  bigatti
// -- added gcd with coefficients in GCDDomain
//
// Revision 1.17  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.16  2010/11/30 11:25:53  bigatti
// -- renamed myIndetName --> myIndetSymbol
//
// Revision 1.15  2010/10/08 14:18:58  bigatti
// -- changed printing style: RingDistrMPolyXXImpl(..) -->  RingDistrMPolyXX(..)
//
// Revision 1.14  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.13  2010/10/01 15:46:31  bigatti
// -- added mySymbolValue
// -- one-liner cleaning
//
// Revision 1.12  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.11  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.10  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.9  2008/06/30 17:15:03  abbott
// Commented out unused parameter names (to avoid annoying compiler warnings).
//
// Revision 1.8  2007/12/21 12:29:08  bigatti
// -- abstract implementation in DenseUPolyRing of myDiv, myIsDivisible, myIsInvertible, myGcd
// -- abstract implementation in DenseUPolyRing of ideal operations
// -- some cleaning
//
// Revision 1.7  2007/12/04 13:44:50  bigatti
// -- fixed problem with myMinCapacity in constructor
// -- commented out myPowerSmallExp in RingDenseUPolyCleanImpl
// -- refined myPowerSmallExp in DenseUPolyRing
//
// Revision 1.6  2007/11/29 17:42:21  bigatti
// -- fixed: myMinCapacity initialization moved up in RingDenseUPolyCleanImpl ctor
//
// Revision 1.5  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/10/19 10:04:23  bigatti
// -- RingDenseUPolyClean now allow to specify the MinCapacity for all
//    coeff vectors (to avoid too many reallocations)
//
// Revision 1.3  2007/10/10 14:02:37  bigatti
// -- added myMulBy1MinusXExp
// -- fixed a few little bugs
//
// Revision 1.2  2007/10/05 16:00:55  bigatti
// -- changed argument of myMulByXExp to std::size_t fon compatibility,
//    but they all have to be reconsidered
//
// Revision 1.1  2007/10/05 15:28:55  bigatti
// -- added abstract class DenseUPolyRing for representing dense
//    univariate polynomials
// -- added concrete class RingDenseUPolyClean, the cleanest
//    implementation
// -- just for testing, still horribly inefficient and incomplete
//
