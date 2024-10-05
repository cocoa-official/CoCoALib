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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingDistrMPolyInlFpPP.C,v 1.54 2024/07/15 17:10:49 bigatti Exp $
// $Log: RingDistrMPolyInlFpPP.C,v $
// Revision 1.54  2024/07/15 17:10:49  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
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
// Revision 1.50  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.49  2020/01/09 16:41:10  bigatti
// -- improved error message for empty list of IndetNames in pseudo ctor
//
// Revision 1.48  2019/03/04 10:38:11  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.47  2018/06/13 15:07:07  abbott
// Summary: Added check for empty list of symbols in pseudo-cto for DMPI and DMPII
//
// Revision 1.46  2018/05/18 12:15:56  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.45  2017/12/01 17:28:21  bigatti
// -- a commented out line (debugging info)
//
// Revision 1.44  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.43  2017/09/06 11:56:29  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.42  2015/12/01 13:11:01  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.41  2015/11/04 10:46:04  abbott
// Summary: Consequential change (after revising SmallFpImpl)
//
// Revision 1.40  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.39  2014/12/11 13:34:42  bigatti
// -- removed all input for doxygen
//
// Revision 1.38  2014/07/28 15:47:55  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble); removed lots of cruft
// Author: JAA
//
// Revision 1.37  2014/07/11 15:46:00  bigatti
// -- removed myOutputSelf (default impl) and added myImplDetails
//
// Revision 1.36  2014/07/04 13:35:42  bigatti
// -- in printing of ring, indeterminates are now more compact
//
// Revision 1.35  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.34  2014/07/02 16:33:14  bigatti
// -- new way of printing rings with ID
//
// Revision 1.33  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.32  2014/04/15 14:25:02  abbott
// Summary: Removed cruft (commented out code which is impl'd in PolyRing)
// Author: JAA
//
// Revision 1.31  2013/05/27 14:08:47  abbott
// Replaced SmallFpElem_t by SmallFpImpl::value_t.
//
// Revision 1.30  2012/10/24 12:19:38  abbott
// Changed return type of myLC.
//
// Revision 1.29  2012/10/17 09:40:16  abbott
// Replaced  RefRingElem  by  RingElem&
// (plus a few consequential changes)
//
// Revision 1.28  2012/09/26 12:29:33  abbott
// Updated to new SmallFpImpl interface.
//
// Revision 1.27  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.26  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.25  2011/11/09 14:10:49  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.24  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.23  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.22  2011/05/19 14:43:49  abbott
// Added include of ZZ.H (previously was already included in some header file).
//
// Revision 1.21  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.20  2010/10/08 14:18:58  bigatti
// -- changed printing style: RingDistrMPolyXXImpl(..) -->  RingDistrMPolyXX(..)
//
// Revision 1.19  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.18  2010/03/05 18:43:48  abbott
// Added pseudo-ctors allowing polynomial rings to be created specifying
// the ordering using a PPOrderingCtor object.
//
// Revision 1.17  2009/12/23 18:53:52  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.16  2009/12/11 11:46:32  abbott
// Changed fn  convert  into  IsConvertible.
// Added template procedure  convert.
// New version because change is not backward compatible.
//
// Revision 1.15  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.14  2009/09/28 16:19:43  bigatti
// -- unique implementation for myDeriv
//
// Revision 1.13  2009/09/25 13:02:09  bigatti
// -- myDiv with one implementation in SparsePolyRing
//
// Revision 1.12  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.11  2008/04/10 15:12:59  bigatti
// -- added myPushBack/Front(RawPtr, ConstRawPtr, PPMonoidElemConstRawPtr)
// -- minor tidying
//
// Revision 1.10  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.9  2007/10/11 16:31:42  bigatti
// -- changed  #ifdef CoCoA_ASSERT  into  CoCoA_DEBUG
//
// Revision 1.8  2007/05/31 16:34:37  bigatti
// -- Changed IsValid (now returns true of false and does not throw an error)
// -- using IsValid for sanity check in PushBack
//
// Revision 1.7  2007/05/31 15:56:33  bigatti
// -- default implementation for IamCommutative, IamIntegralDomain,
//    IamGCDDomain, IamField, myCharacteristic  in PolyRing
// -- default implementation for mySymbols  in SparsePolyRing
// -- cleaned up pseudo-ctors
//
// Revision 1.5  2007/05/21 14:50:56  bigatti
// -- myPushFront and myPushBack now accept zero coefficient
//
// Revision 1.4  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.3  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.2  2007/03/09 18:56:56  bigatti
// -- added Tmp prefix to Groebner related files
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.28  2007/03/08 17:43:10  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.27  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.26  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.25  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.24  2007/01/20 14:07:25  bigatti
// -- moved code for homomorphism into common implementation in SparsePolyRing
//
// Revision 1.23  2007/01/15 15:47:57  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.22  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.21  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.20  2006/12/06 16:19:11  cocoa
// Minor corrections to comments.
//
// Revision 1.19  2006/11/23 18:01:53  cocoa
// -- moved printing functions in unified implementation in SparsePolyRing
// -- simplified "output(f)" for debugging only
//
// Revision 1.18  2006/11/22 17:51:31  cocoa
// -- moved printing functions into unified implementation in SparsePolyRing
//
// Revision 1.17  2006/11/22 15:01:20  cocoa
// -- sorted #include
//
// Revision 1.16  2006/11/22 14:55:50  cocoa
// -- added #include "CoCoA/SparsePolyRing.H"
//
// Revision 1.15  2006/11/21 18:09:23  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.14  2006/11/14 17:18:35  cocoa
// -- added comment about myRefCountZero()
//
// Revision 1.13  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPoly
//
// Revision 1.12  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.11  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.10  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.9  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.8  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.7  2006/10/06 10:15:52  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.6  2006/09/27 14:31:59  cocoa
// -- improved checks for CoeffRing in contructor (must be RingFp)
//
// Revision 1.5  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.4  2006/07/20 17:06:08  cocoa
// -- moved myStdDeg into SparsePolyRing
//
// Revision 1.3  2006/07/17 19:49:34  cocoa
// Added some stronger arg checks (to myPushFront and myPushBack).
//
// Revision 1.2  2006/06/08 16:45:28  cocoa
// -- RingDistrMPoly*.H  have been "moved" into RingDistrMPoly*.C
// -- some coding conventions fixed in DistrMPoly*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.25  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.24  2006/05/29 13:20:33  cocoa
// -- added implementation of "intersect" for SparsePolyRing
//
// Revision 1.23  2006/05/12 17:06:09  cocoa
// -- moved myIsUnit, myGcd into SparsePolyRing (common implementation)
//
// Revision 1.22  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.21  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.20  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.19  2006/04/21 16:47:06  cocoa
// -- new syntax for ComputeGBasis by Max
//
// Revision 1.18  2006/04/05 12:43:13  cocoa
// -- fixed: RingDistrMPolyInlFpPPImpl::HomImpl::myApply
//
// Revision 1.17  2006/03/30 17:19:48  cocoa
// -- removed "static" before typedef InlineFpElem_t
//
// Revision 1.16  2006/03/30 16:59:27  cocoa
// -- changed misleading name: InlineCoeffRing --> InlineCoeffImpl
// -- new: implementation for homomorphisms
// -- rearrangement of code to mimic RingDistrMPolyInlPP
//
// Revision 1.15  2006/03/21 09:43:13  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.14  2006/03/20 17:27:42  cocoa
// -- changed in DistrMPolyInlFpPP: myMul, myDiv --> myMulByCoeff, myMulByPP, myDivByCoeff
//
// Revision 1.13  2006/03/17 18:10:27  cocoa
// -- changed: myMul --> myMulByPP
//
// Revision 1.12  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.11  2006/03/14 15:01:49  cocoa
// Improved the implementation of ring member fns for computing powers.
// Should keep Intel C++ compiler quieter too.
//
// Revision 1.10  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.9  2006/03/07 10:06:12  cocoa
// -- fixed: PPMonoidElem LPP(rawf) now returns ConstRefPPMonoidElem
//
// Revision 1.8  2006/03/07 09:57:10  cocoa
// -- changed only some commented out code in myApply (to be implemented)
//
// Revision 1.7  2006/02/20 22:41:19  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.6  2006/02/14 16:19:37  cocoa
// -- defined "IdealImpl::contains"
//
// Revision 1.5  2006/02/13 13:56:50  cocoa
// -- changed: "GReductor" --> "ComputeGBasis"
//
// Revision 1.4  2006/01/19 18:03:53  cocoa
// -- fixed coding conventions for myGBasis, myGBasisValue
// -- fixed myReduceMod (NF - tested)
//
// Revision 1.3  2006/01/19 16:34:42  cocoa
// -- added NF, myReduceMod functions (not yet tested)
//
// Revision 1.2  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.8  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.7  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.6  2005/07/15 16:34:33  cocoa
// Added iterators for sparse polynomials.
// The code compiles (and the old tests still run).
// It'd Friday evening -- I'm going home before
// getting any ideas about making the iterator code run.
//
// Revision 1.5  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.4  2005/07/07 16:40:24  cocoa
// -- added: ASSERT for gcd(0,0)
//
// Revision 1.3  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
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
// Revision 1.13  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.12  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.11  2004/11/11 14:32:01  cocoa
// -- minor changes for doxygen
// -- change: cout --> GlobalLogput()
//
// Revision 1.10  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.9  2004/11/08 13:54:48  cocoa
// -- changed calls to ZZ (after changes to ZZ.H)
//
// Revision 1.8  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.7  2004/10/29 16:14:42  cocoa
// -- changed PPOrdering calls to OrdvArith
// -- new default monoid: NewPPMonoidSafe with NewStdDegRevLexOrdering
//
// Revision 1.5  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.4  2004/07/20 15:37:08  cocoa
// Minor fix for some errors which slipped through the net...
//
// Revision 1.3  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.2  2004/07/16 10:11:34  cocoa
// -- now using the new class SmallFpImpl (or SmallFpLogImpl)
// -- updated with "my" coding convenctions
// -- NYI: LC and LCRaw
//
// Revision 1.1  2004/06/25 16:03:06  cocoa
// -- first import
//
