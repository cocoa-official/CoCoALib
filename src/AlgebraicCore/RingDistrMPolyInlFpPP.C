//   Copyright (c)  2005-2007,2010,2021  John Abbott and Anna M. Bigatti

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


// Source code for RingDistrMPolyInlFpPPImpl

#include "CoCoA/RingDistrMPolyInlFpPP.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DistrMPolyInlFpPP.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoidOv.H" // includes "CoCoA/symbol.H"
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::sort;
using std::copy;
#include <iostream>
using std::ostream;
#include <memory>
using std::unique_ptr;
#include <new>
//for placement new
//#include <vector> // included in RingDistrMPolyInlFpPP.H
using std::vector;


namespace CoCoA
{

  class RingDistrMPolyInlFpPPImpl: public SparsePolyRingBase
  {
  private:
    typedef DistrMPolyInlFpPP value_t; // DistrMPolyInlFpPP is the actual type of the values in a RingDistrMPolyInlFpPP
    static value_t& import(RingElemRawPtr rawf);
    static const value_t& import(RingElemConstRawPtr rawf);

  public:
    RingDistrMPolyInlFpPPImpl(const ring& R, const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    ~RingDistrMPolyInlFpPPImpl();

  private: // Data members of RingDistrMPolyInlFpPPImpl
    const ring myCoeffRingValue;  ///< the coefficient ring (used for input and output)
    const DistrMPolyInlFpPP::InlineFpImpl myInlineCoeffImplValue;  ///< to be used instead of ring for inlining operations on coefficients
    const PPMonoid myPPMValue; ///< the monoid of the power-products
    const OrdvArith::reference myOrdvArith;
    mutable MemPool myDMPPool; ///< memory manager for polynomials
    long myNumIndetsValue; ///< number of indeteminates
    mutable MemPool mySummandPool; ///< memory manager for summands; MemPool MUST COME BEFORE myZeroPtr, myOnePtr, and myIndetVector!
    std::unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    std::unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    std::vector<RingElem> myIndetVector;

  public:  // functions which every ring must implement
    ConstRefRingElem myZero() const override  {return *myZeroPtr;}
    ConstRefRingElem myOne() const override  {return *myOnePtr;}
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
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x-y
    std::string myImplDetails() const override  {return "RingDistrMPolyInlFpPP";}
    bool myIsZero(ConstRawPtr rawx) const override;                        // x == 0
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;     // x == y

    // functions which every PolyRing must implement
    long myNumIndets() const override  {return myNumIndetsValue;}
    const ring& myCoeffRing() const override  {return myCoeffRingValue;}
    const std::vector<RingElem>& myIndets() const override  {return myIndetVector;}
    void myIndetPower(RawPtr rawf, long var, long exp) const override;

    ///@name Simple functions on polynomials
    //@{
    long myNumTerms(ConstRawPtr rawf) const override;
    bool myIsMonomial(ConstRawPtr rawf) const override;
    RingElemAlias myLC(ConstRawPtr rawf) const override;
    void myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< WEAK EXCEPTION GUARANTEE
    bool myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< WEAK EXCEPTION GUARANTEE
    //@}

    //----------------------------------------------------------------------
    // Functions which every SparsePolyRing must implement:
    //----------------------------------------------------------------------

    const PPMonoid& myPPM() const override  {return myPPMValue;}

    ///@name   Functions for creating/building polynomials
    //@{
    RingElem myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const override;
    SparsePolyIter myBeginIter(ConstRawPtr rawf) const override;
    SparsePolyIter myEndIter(ConstRawPtr rawf) const override;

    void myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const override;
    void myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const override;
    void myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const override;
    void myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const override;
    //@}

    ///@name Special functions on polynomials needed for implementing Buchberger's Algorithm
    //@{
    ConstRefPPMonoidElem myLPP(ConstRawPtr rawf) const override;
    void myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const override;
    bool myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const override; ///< f+=LM(g); g-=LM(g); assumes LPP(f)==LPP(g); returns LC(f)+LC(g)==0
    void myMoveLMToFront(RawPtr rawf, RawPtr rawg) const override;
    void myMoveLMToBack(RawPtr rawf, RawPtr rawg) const override;
    void myDeleteLM(RawPtr rawf) const override; // ????? right interface
    void myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const override; ///< lhs=div(LM(f),LM(g)); assumes f!=0,g!=0
    int  myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const override; ///< cmp(LPP(f),LPP(g)); assumes f!=0,g!=0
    void myAddClear(RawPtr rawf, RawPtr rawg) const override; ///< f+=g; g=0;
    void myAppendClear(RawPtr rawf, RawPtr rawg) const override; ///< f+=g; g=0; appends g to f with no checks

    void myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const override; ///<  f += LM(h)*g
    void myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag) const override; ///<  f += LM(h)*g
    void myReductionStep(RawPtr rawf, ConstRawPtr rawg) const override;
    // ??? aggiungere coefficiente
    void myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& FScale) const override;
    // should it all be in ReductionStep ??? ANNA
    //@}


// ***ATTENTION*** [20140726] is it worth having this special case code???
  // private: // Special homomorphism class for this type of ring.
  //   class CoeffEmbeddingHomImpl: public SparsePolyRingBase::CoeffEmbeddingHomImpl
  //   {
  //   public:
  //     CoeffEmbeddingHomImpl(const SparsePolyRing& P);
  //     virtual void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const;
  //   };

  };

  //----------------------------------------------------------------------

  typedef DistrMPolyInlFpPP::InlineFpElem_t InlineFpElem_t; // to save some typing


  inline RingDistrMPolyInlFpPPImpl::value_t& RingDistrMPolyInlFpPPImpl::import(RingElemRawPtr rawf)
  {
    return *static_cast<value_t*>(rawf.myRawPtr());
  }

  inline const RingDistrMPolyInlFpPPImpl::value_t& RingDistrMPolyInlFpPPImpl::import(RingElemConstRawPtr rawf)
  {
    return *static_cast<const value_t*>(rawf.myRawPtr());
  }


  RingDistrMPolyInlFpPPImpl::RingDistrMPolyInlFpPPImpl(const ring& R, const std::vector<symbol>& IndetNames, const PPOrdering& ord):
    myCoeffRingValue(R),
    myInlineCoeffImplValue(ConvertTo<long>(characteristic(R))),
    myPPMValue(NewPPMonoidOv(IndetNames, ord)),
    myOrdvArith(NewOrdvArith(ord)),
    myDMPPool(sizeof(DistrMPolyInlFpPP), "RingDistrMPolyInlFpPPImpl::myDMPPool"),
    myNumIndetsValue(NumIndets(ord)),
    mySummandPool(DistrMPolyInlFpPP::SummandSize(R, myOrdvArith), "RingDistrMPolyInlFpPPImpl::mySummandPool")
  {
    //    std::cout << "------RingDistrMPolyInlFpPPImpl:NewOrdvArith-called" << std::endl;
    CoCoA_ASSERT(IsCommutative(myCoeffRingValue));
    CoCoA_ASSERT(len(IndetNames) == myNumIndetsValue);
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myIndetVector.resize(myNumIndetsValue, *myZeroPtr);
    vector<long> expv(myNumIndetsValue);
    for (long i=0; i < myNumIndetsValue; ++i)
    {
      expv[i] = 1;
      myPushFront(raw(myIndetVector[i]), raw(one(R)), expv);
      expv[i] = 0;
    }
    myRefCountZero(); // otherwise it is 2 + NumIndets and won't be destroyed
  }


  RingDistrMPolyInlFpPPImpl::~RingDistrMPolyInlFpPPImpl()
  {}


  //----------------------------------------------------------------------
  // Functions which every ring must implement:
  //----------------------------------------------------------------------

  RingElemRawPtr RingDistrMPolyInlFpPPImpl::myNew() const
  {
    void* ptr = myDMPPool.alloc();
    new(ptr) DistrMPolyInlFpPP(myInlineCoeffImplValue, myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool); // placement new
    return RingElemRawPtr(ptr);
  }


  RingElemRawPtr RingDistrMPolyInlFpPPImpl::myNew(const MachineInt& n) const
  {
    if (IsZero(n)) return myNew();
    unique_ptr<DistrMPolyInlFpPP> ans(new(myDMPPool.alloc()) DistrMPolyInlFpPP(myInlineCoeffImplValue, myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool)); // placement new
    *ans = n;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDistrMPolyInlFpPPImpl::myNew(const BigInt& N) const
  {
    if (N == 0) return myNew();
    unique_ptr<DistrMPolyInlFpPP> ans(new(myDMPPool.alloc()) DistrMPolyInlFpPP(myInlineCoeffImplValue, myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool)); // placement new
    *ans = N;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDistrMPolyInlFpPPImpl::myNew(ConstRawPtr rawcopy) const
  {
    unique_ptr<DistrMPolyInlFpPP> ans(new(myDMPPool.alloc()) DistrMPolyInlFpPP(myInlineCoeffImplValue, myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool)); // placement new
    *ans = import(rawcopy);
    return RingElemRawPtr(ans.release());
  }


  void RingDistrMPolyInlFpPPImpl::myDelete(RawPtr rawx) const
  {
    import(rawx).~DistrMPolyInlFpPP();
    myDMPPool.free(rawx.myRawPtr());
  }


  void RingDistrMPolyInlFpPPImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    swap(import(rawx), import(rawy));
  }


  void RingDistrMPolyInlFpPPImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingDistrMPolyInlFpPPImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = n;
  }


  void RingDistrMPolyInlFpPPImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = N;
  }


  void RingDistrMPolyInlFpPPImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs).myAssignZero();
  }


  void RingDistrMPolyInlFpPPImpl::myRecvTwinFloat(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "RingDistrMPolyInlFpPPImpl::myRecvTwinFloat");
  }


  void RingDistrMPolyInlFpPPImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
    import(rawlhs).myNegate();
  }


  void RingDistrMPolyInlFpPPImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    add(import(rawlhs), import(rawx), import(rawy));
  }


  void RingDistrMPolyInlFpPPImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    sub(import(rawlhs), import(rawx), import(rawy));
  }


  bool RingDistrMPolyInlFpPPImpl::myIsZero(ConstRawPtr rawx) const
  {
    return IsZero(import(rawx));
  }


  bool RingDistrMPolyInlFpPPImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return IsEqual(import(rawx), import(rawy));
  }


  //----------------------------------------------------------------------
  // Functions which every PolyRing must implement
  //----------------------------------------------------------------------

  void RingDistrMPolyInlFpPPImpl::myIndetPower(RawPtr rawf, long var, long exp) const
  {
    CoCoA_ASSERT(var < myNumIndets());
    CoCoA_ASSERT(exp >= 0);
    RingElem ans(ring(this));
    vector<long> expv(myNumIndets()); // wasteful new/delete
    expv[var] = exp;
    import(raw(ans)).myPushFront(one(SmallFp), expv);
    mySwap(raw(ans), rawf); // do it this way to be exception clean
  }


  long RingDistrMPolyInlFpPPImpl::myNumTerms(ConstRawPtr rawx) const
  { return NumTerms(import(rawx)); }

  bool RingDistrMPolyInlFpPPImpl::myIsMonomial(ConstRawPtr rawf) const
  { return IsMonomial(import(rawf)); }


  RingElemAlias RingDistrMPolyInlFpPPImpl::myLC(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return RingElemAlias(myCoeffRing(), RingElemRawPtr(const_cast<InlineFpElem_t*>(&LC(import(rawf)))));
  }


  void RingDistrMPolyInlFpPPImpl::myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const // EXCEPTION SAFE
  {
    BigInt N;
    myCoeffRingValue->myIsInteger(N, rawc);  // must return true
    InlineFpElem_t c = myInlineCoeffImplValue.myReduce(N);
    import(rawf).myMulByCoeff(c);
  }


  bool RingDistrMPolyInlFpPPImpl::myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const // EXCEPTION SAFE
  {
    BigInt N;
    myCoeffRingValue->myIsInteger(N, rawc);  // must return true
    InlineFpElem_t c = myInlineCoeffImplValue.myReduce(N);
    return import(rawf).myDivByCoeff(c);
  }


  //----------------------------------------------------------------------
  // Functions which every SparsePolyRing must implement:
  //----------------------------------------------------------------------


  RingElem RingDistrMPolyInlFpPPImpl::myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(!myCoeffRing()->myIsZero(rawc));
    RingElem ans(ring(this));
    myPushFront(raw(ans), rawc, rawpp);
    return ans;
  }


  SparsePolyIter RingDistrMPolyInlFpPPImpl::myBeginIter(ConstRawPtr rawf) const
  { return SparsePolyIter(new DistrMPolyInlFpPP::iter(import(rawf))); }

  SparsePolyIter RingDistrMPolyInlFpPPImpl::myEndIter(ConstRawPtr rawf) const
  { return SparsePolyIter(new DistrMPolyInlFpPP::iter(import(rawf), 0)); }


  void RingDistrMPolyInlFpPPImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    BigInt N;
    myCoeffRingValue->myIsInteger(N, rawc);  // Necessarily returns true.
    InlineFpElem_t c = myInlineCoeffImplValue.myReduce(N);
    import(rawf).myPushFront(c, expv);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyInlFpPPImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    BigInt N;
    myCoeffRingValue->myIsInteger(N, rawc);  // Necessarily returns true.
    InlineFpElem_t c = myInlineCoeffImplValue.myReduce(N);
    import(rawf).myPushBack(c, expv);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyInlFpPPImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    BigInt N;
    myCoeffRingValue->myIsInteger(N, rawc);  // Necessarily returns true.
    InlineFpElem_t c = myInlineCoeffImplValue.myReduce(N);
    import(rawf).myPushFront(c, rawpp);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyInlFpPPImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    BigInt N;
    myCoeffRingValue->myIsInteger(N, rawc);  // Necessarily returns true.
    InlineFpElem_t c = myInlineCoeffImplValue.myReduce(N);
    import(rawf).myPushBack(c, rawpp);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  ConstRefPPMonoidElem RingDistrMPolyInlFpPPImpl::myLPP(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return LPP(import(rawf));
  }


  void RingDistrMPolyInlFpPPImpl::myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const
  { import(rawf).myMulByPP(rawpp); }

  bool RingDistrMPolyInlFpPPImpl::myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const
  { return IsZeroAddLCs(import(rawf), import(rawg)); }


  void RingDistrMPolyInlFpPPImpl::myMoveLMToFront(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    MoveLMToFront(import(rawf), import(rawg));
  }


  void RingDistrMPolyInlFpPPImpl::myMoveLMToBack(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    MoveLMToBack(import(rawf), import(rawg));
  }


  void RingDistrMPolyInlFpPPImpl::myDeleteLM(RawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    import(rawf).myDeleteLM();
  }


  void RingDistrMPolyInlFpPPImpl::myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    DivLM(import(rawlhs), import(rawf), import(rawg));
  }


  int  RingDistrMPolyInlFpPPImpl::myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    return CmpLPP(import(rawf), import(rawg));
  }


  void RingDistrMPolyInlFpPPImpl::myAddClear(RawPtr rawf, RawPtr rawg) const
  { import(rawf).myAddClear(import(rawg)); }


  void RingDistrMPolyInlFpPPImpl::myAppendClear(RawPtr rawf, RawPtr rawg) const
  {
#ifdef CoCoA_DEBUG
    if (!myIsZero(rawf) && !myIsZero(rawg))
    {
      for (SparsePolyIter it=myBeginIter(rawf); !IsEnded(it); ++it)  // INEFFICIENT   SLUG????
        CoCoA_ASSERT(PP(it) > myLPP(rawg));
    }
#endif
    import(rawf).myAppendClear(import(rawg));
  }


  void RingDistrMPolyInlFpPPImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const //???? delete me???
  { import(rawf).myAddMulLM(import(rawh), import(rawg), /* SkipLMg = */ false); }

  void RingDistrMPolyInlFpPPImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag skip) const //???? delete me???
  { import(rawf).myAddMulLM(import(rawh), import(rawg), skip==SkipLMg); }

  void RingDistrMPolyInlFpPPImpl::myReductionStep(RawPtr rawf, ConstRawPtr rawg) const
  { import(rawf).myReductionStep(import(rawg)); }

  void RingDistrMPolyInlFpPPImpl::myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& fscale) const
  { import(rawf).myReductionStepGCD(import(rawg), fscale); }


// ***ATTENTION*** CoeffEmbeddingHomImpl
// Special case code; should be a bit more efficient that general case -- is it really worth it???
  // //---------------------------------------------------------------------------
  // // Functions for the class RingDistrMPolyInlFpPPImpl::CoeffEmbeddingHom

  // RingDistrMPolyInlFpPPImpl::CoeffEmbeddingHomImpl::CoeffEmbeddingHomImpl(const SparsePolyRing& P):
  //   SparsePolyRingBase::CoeffEmbeddingHomImpl(P)
  // {}


  // void RingDistrMPolyInlFpPPImpl::CoeffEmbeddingHomImpl::myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const
  // {
  //   BigInt x;   //??? wasteful ctor/dtor???
  //   myDomain->myIsInteger(x, rawarg);  // must return true
  //   myCodomain->myAssign(rawimage, x);  // I think this is "no throw" in this case ???
  // }


  //----------------------------------------------------------------------
  // Pseudo-ctors for (sparse) polynomial rings.

  namespace  {  // anonymous
    void CheckCoeffRing(const ring& K)
    {
      if (!IsRingFp(K)) // check that CoeffRing is really a SmallFpImpl
        CoCoA_THROW_ERROR(ERR::NotQuotientRing, "NewPolyRing_DMPII pseudo ctor");
    }
    

    void CheckIndets(const ring& K, const std::vector<symbol>& IndetNames)
    {
      if (IndetNames.empty())
        CoCoA_THROW_ERROR(ERR::ReqNonEmpty, "NewPolyRing_DMPII pseudo ctor");
      if (!AreGoodIndetNames(K, IndetNames))
        CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPolyRing_DMPII pseudo ctor");
    }


    void CheckNumIndets(long n)
    {
      if (n <=0 )
        CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPolyRing_DMPII pseudo ctor");
    }


    // Apparently not used (2022-02-07)
    // void CheckNumIndets(long n, const std::vector<symbol>& IndetNames)
    // {
    //   CheckNumIndets(n);
    //   if (n != len(IndetNames))
    //     CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPolyRing_DMPII pseudo ctor");
    // }
    
  }  // anonymous namespace


  SparsePolyRing NewPolyRing_DMPII(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
    CheckCoeffRing(CoeffRing);
    CheckIndets(CoeffRing, IndetNames);
    return SparsePolyRing(new RingDistrMPolyInlFpPPImpl(CoeffRing, IndetNames, ord));
  }


  SparsePolyRing NewPolyRing_DMPII(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  {
    CheckCoeffRing(CoeffRing);
    CheckIndets(CoeffRing, IndetNames);
    return SparsePolyRing(new RingDistrMPolyInlFpPPImpl(CoeffRing, IndetNames, OrdCtor(len(IndetNames))));
  }


  SparsePolyRing NewPolyRing_DMPII(const ring& CoeffRing, const std::vector<symbol>& IndetNames)
  {
    return NewPolyRing_DMPII(CoeffRing, IndetNames, StdDegRevLex(len(IndetNames)));
  }


  SparsePolyRing NewPolyRing_DMPII(const ring& CoeffRing, long NumIndets)
  {
    CheckNumIndets(NumIndets);
    return NewPolyRing_DMPII(CoeffRing, SymbolRange("x",0,NumIndets-1), StdDegRevLex(NumIndets));
  }


  SparsePolyRing NewPolyRing_DMPII(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& OrdCtor)
  {
    CheckNumIndets(NumIndets);
    return NewPolyRing_DMPII(CoeffRing, SymbolRange("x",0,NumIndets-1), OrdCtor(NumIndets));
  }


} // end of namespace CoCoA
