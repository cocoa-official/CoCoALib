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


// Source code for RingDistrMPolyCleanImpl

#include "CoCoA/RingDistrMPolyClean.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DistrMPolyClean.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <iostream>
//for operator <<
#include <memory>
using std::unique_ptr;
#include <new>
//for placement new
#include <vector>
using std::vector;


namespace CoCoA
{

  class RingDistrMPolyCleanImpl: public SparsePolyRingBase
  {
  private:
    typedef DistrMPolyClean value_t; // DistrMPolyClean is the actual type of the values in a RingDistrMPolyCleanImpl
    static value_t& import(RingElemRawPtr rawf);
    static const value_t& import(RingElemConstRawPtr rawf);

  public:
    RingDistrMPolyCleanImpl(const ring& R, const PPMonoid& PPM);
    ~RingDistrMPolyCleanImpl();

  private: // Data members of RingDistrMPolyCleanImpl
    const ring myCoeffRingValue;  ///< the coefficient ring
    const PPMonoid myPPMValue; ///< the monoid of the power-products
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
    std::string myImplDetails() const override  {return "RingDistrMPolyClean";}
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

  };

  //----------------------------------------------------------------------

  inline RingDistrMPolyCleanImpl::value_t& RingDistrMPolyCleanImpl::import(RingElemRawPtr rawf)
  {
    return *static_cast<value_t*>(rawf.myRawPtr());
  }

  inline const RingDistrMPolyCleanImpl::value_t& RingDistrMPolyCleanImpl::import(RingElemConstRawPtr rawf)
  {
    return *static_cast<const value_t*>(rawf.myRawPtr());
  }


  RingDistrMPolyCleanImpl::RingDistrMPolyCleanImpl(const ring& R, const PPMonoid& PPM):
    myCoeffRingValue(R),
    myPPMValue(PPM),
    myDMPPool(sizeof(DistrMPolyClean), "RingDistrMPolyCleanImpl::myDMPPool"),
    myNumIndetsValue(NumIndets(PPM)),
    mySummandPool(DistrMPolyClean::ourSummandSize(R, PPM), "RingDistrMPolyCleanImpl::mySummandPool")
  {
    CoCoA_ASSERT(IsCommutative(myCoeffRingValue));
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myIndetVector.resize(myNumIndetsValue, *myZeroPtr);
    // Now fill vector of "indets".
    const vector<PPMonoidElem>& x = indets(PPM);
    for (long i=0; i < myNumIndetsValue; ++i)
      myPushFront(raw(myIndetVector[i]), raw(one(R)), raw(x[i]));
    myRefCountZero(); // otherwise it is 2 + NumIndets and won't be destroyed
  }


  RingDistrMPolyCleanImpl::~RingDistrMPolyCleanImpl()
  {}


  //----------------------------------------------------------------------
  // Functions which every ring must implement:
  //----------------------------------------------------------------------

  RingElemRawPtr RingDistrMPolyCleanImpl::myNew() const
  {
    void* ptr = myDMPPool.alloc();
    new(ptr) DistrMPolyClean(myCoeffRingValue, myPPMValue, mySummandPool); // placement new
    return RingElemRawPtr(ptr);
  }


  RingElemRawPtr RingDistrMPolyCleanImpl::myNew(const MachineInt& n) const
  {
    if (IsZero(n)) return myNew();
    unique_ptr<DistrMPolyClean> ans(new(myDMPPool.alloc()) DistrMPolyClean(myCoeffRingValue, myPPMValue, mySummandPool)); // placement new
    *ans = n;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDistrMPolyCleanImpl::myNew(const BigInt& N) const
  {
    if (N == 0) return myNew();  // not really necessary
    unique_ptr<DistrMPolyClean> ans(new(myDMPPool.alloc()) DistrMPolyClean(myCoeffRingValue, myPPMValue, mySummandPool)); // placement new
    *ans = N;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDistrMPolyCleanImpl::myNew(ConstRawPtr rawcopy) const
  {
    unique_ptr<DistrMPolyClean> ans(new(myDMPPool.alloc()) DistrMPolyClean(myCoeffRingValue, myPPMValue, mySummandPool)); // placement new
    *ans = import(rawcopy);
    return RingElemRawPtr(ans.release());
  }


  void RingDistrMPolyCleanImpl::myDelete(RawPtr rawx) const
  {
    import(rawx).~DistrMPolyClean();
    myDMPPool.free(rawx.myRawPtr());
  }


  void RingDistrMPolyCleanImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    swap(import(rawx), import(rawy));
  }


  void RingDistrMPolyCleanImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingDistrMPolyCleanImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = n;
  }


  void RingDistrMPolyCleanImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = N;
  }


  void RingDistrMPolyCleanImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs).myAssignZero();
  }


  void RingDistrMPolyCleanImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myCoeffRingValue));
    RingElem tmp(myCoeffRingValue);
    myCoeffRingValue->myRecvTwinFloat(raw(tmp), rawx);
    myCoeffEmbeddingHomCtor()->myApply(rawlhs, raw(tmp));
  }


  void RingDistrMPolyCleanImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
    import(rawlhs).myNegate();
  }


  void RingDistrMPolyCleanImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    add(import(rawlhs), import(rawx), import(rawy));
  }


  void RingDistrMPolyCleanImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    sub(import(rawlhs), import(rawx), import(rawy));
  }


  bool RingDistrMPolyCleanImpl::myIsZero(ConstRawPtr rawx) const
  {
    return IsZero(import(rawx));
  }


  bool RingDistrMPolyCleanImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return IsEqual(import(rawx), import(rawy));
  }


  //----------------------------------------------------------------------
  // Functions which every PolyRing must implement
  //----------------------------------------------------------------------

  void RingDistrMPolyCleanImpl::myIndetPower(RawPtr rawf, long var, long exp) const
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(var < myNumIndets());
    RingElem ans(ring(this));
    vector<long> expv(myNumIndets()); // wasteful new/delete
    expv[var] = exp;
    import(raw(ans)).myPushFront(raw(one(myCoeffRingValue)), expv);
    mySwap(raw(ans), rawf); // do it this way to be exception clean
  }


  long RingDistrMPolyCleanImpl::myNumTerms(ConstRawPtr rawx) const
  {
    return NumTerms(import(rawx));
  }


  bool RingDistrMPolyCleanImpl::myIsMonomial(ConstRawPtr rawf) const
  {
    return IsMonomial(import(rawf));
  }


  RingElemAlias RingDistrMPolyCleanImpl::myLC(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return LC(import(rawf));
  }


  void RingDistrMPolyCleanImpl::myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const // EXCEPTION SAFE
  {
    import(rawf).myMulByCoeff(rawc);
  }


  bool RingDistrMPolyCleanImpl::myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const // EXCEPTION SAFE
  {
    return import(rawf).myDivByCoeff(rawc);
  }


  //----------------------------------------------------------------------
  // Functions which every SparsePolyRing must implement:
  //----------------------------------------------------------------------


  RingElem RingDistrMPolyCleanImpl::myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(!myCoeffRing()->myIsZero(rawc));
    RingElem ans(ring(this));
    myPushFront(raw(ans), rawc, rawpp);
    return ans;
  }


  SparsePolyIter RingDistrMPolyCleanImpl::myBeginIter(ConstRawPtr rawf) const
  {
    return SparsePolyIter(new DistrMPolyClean::iter(import(rawf)));
  }


  SparsePolyIter RingDistrMPolyCleanImpl::myEndIter(ConstRawPtr rawf) const
  {
    return SparsePolyIter(new DistrMPolyClean::iter(import(rawf), 0));
  }


  void RingDistrMPolyCleanImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    import(rawf).myPushFront(rawc, expv);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyCleanImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    import(rawf).myPushBack(rawc, expv);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyCleanImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    import(rawf).myPushFront(rawc, rawpp);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyCleanImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    import(rawf).myPushBack(rawc, rawpp);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  ConstRefPPMonoidElem RingDistrMPolyCleanImpl::myLPP(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return LPP(import(rawf));
  }


  void RingDistrMPolyCleanImpl::myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const
  {
    import(rawf).myMulByPP(rawpp);
  }


  bool RingDistrMPolyCleanImpl::myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const
  {
    return IsZeroAddLCs(import(rawf), import(rawg));
  }


  void RingDistrMPolyCleanImpl::myMoveLMToFront(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    MoveLMToFront(import(rawf), import(rawg));
  }


  void RingDistrMPolyCleanImpl::myMoveLMToBack(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    MoveLMToBack(import(rawf), import(rawg));
  }


  void RingDistrMPolyCleanImpl::myDeleteLM(RawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    import(rawf).myDeleteLM();
  }


  void RingDistrMPolyCleanImpl::myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    DivLM(import(rawlhs), import(rawf), import(rawg));
  }


  int RingDistrMPolyCleanImpl::myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    return CmpLPP(import(rawf), import(rawg));
  }


  void RingDistrMPolyCleanImpl::myAddClear(RawPtr rawf, RawPtr rawg) const
  {
    import(rawf).myAddClear(import(rawg));
  }


  void RingDistrMPolyCleanImpl::myAppendClear(RawPtr rawf, RawPtr rawg) const
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


  void RingDistrMPolyCleanImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const //???? delete me???
  {
    import(rawf).myAddMulLM(import(rawh), import(rawg), /* SkipLMg = */ false);
  }


  void RingDistrMPolyCleanImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag skip) const //???? delete me???
  {
    import(rawf).myAddMulLM(import(rawh), import(rawg), skip==SkipLMg);
  }


  void RingDistrMPolyCleanImpl::myReductionStep(RawPtr rawf, ConstRawPtr rawg) const
  {
    import(rawf).myReductionStep(import(rawg));
  }


  void RingDistrMPolyCleanImpl::myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& fscale) const
  {
    import(rawf).myReductionStepGCD(import(rawg), fscale);
  }


  //----------------------------------------------------------------------
  // Pseudo-ctors for (sparse) polynomial rings.

  namespace  {  // anonymous
    void CheckCoeffRing(const ring& K)
    {
      if (!IsCommutative(K))
        CoCoA_THROW_ERROR2(ERR::NotCommutative, "NewPolyRing_DMP pseudo ctor");
    }
    

    void CheckIndets(const ring& K, const std::vector<symbol>& IndetNames)
    {
      if (IndetNames.empty())
        CoCoA_THROW_ERROR2(ERR::ReqNonEmpty, "NewPolyRing_DMPI pseudo ctor");
      if (!AreGoodIndetNames(K, IndetNames))
        CoCoA_THROW_ERROR2(ERR::BadIndetNames, "NewPolyRing_DMPI pseudo ctor");
    }


    void CheckNumIndets(long n)
    {
      if (n <= 0)
        CoCoA_THROW_ERROR2(ERR::BadNumIndets, "NewPolyRing_DMP pseudo ctor");
    }


    void CheckNumIndets(long n, const std::vector<symbol>& IndetNames)
    {
      CheckNumIndets(n);
      if (n != len(IndetNames))
        CoCoA_THROW_ERROR2(ERR::BadNumIndets, "NewPolyRing_DMP pseudo ctor");
    }
    
  }  // anonymous namespace


  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const PPMonoid& PPM)
  {
    CheckCoeffRing(CoeffRing);
    CheckIndets(CoeffRing, symbols(PPM));
    return SparsePolyRing(new RingDistrMPolyCleanImpl(CoeffRing, PPM));
  }


  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
    CheckNumIndets(NumIndets(ord), IndetNames);
    return NewPolyRing_DMP(CoeffRing, NewPPMonoidEv(IndetNames, ord));
  }


  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  {
    return NewPolyRing_DMP(CoeffRing, NewPPMonoidEv(IndetNames, OrdCtor));
  }


  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, const std::vector<symbol>& IndetNames)
  {
    return NewPolyRing_DMP(CoeffRing, NewPPMonoidEv(IndetNames,StdDegRevLex(len(IndetNames))));
  }


  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, long NumIndets)
  {
    CheckNumIndets(NumIndets);
    return NewPolyRing_DMP(CoeffRing, NewPPMonoidEv(SymbolRange("x",0,NumIndets-1), StdDegRevLex(NumIndets)));
  }


  SparsePolyRing NewPolyRing_DMP(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& OrdCtor)
  {
    CheckNumIndets(NumIndets);
    return NewPolyRing_DMP(CoeffRing,NewPPMonoidEv(SymbolRange("x",0,NumIndets-1),OrdCtor));
  }


} // end of namespace CoCoA
