//   Copyright (c)  2001-2017,2021  John Abbott and Anna M. Bigatti
//   Author:  2005-2010  John Abbott

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


// Source code for RingDistrMPolyInlPPImpl

#include "CoCoA/RingDistrMPolyInlPP.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DistrMPolyInlPP.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoidOv.H" // includes "CoCoA/symbol.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/assert.H"
#include "CoCoA/config.H" // for myIndetPower
#include "CoCoA/convert.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::max;
using std::sort;
using std::copy;
#include <iostream>
using std::ostream;
#include <memory>
using std::unique_ptr;
#include <limits>
using std::numeric_limits; // needed only if CoCoA_DEBUG is set.
#include <new>
//for placement new
//#include <vector> // included in RingDistrMPolyInlPP.H
using std::vector;


namespace CoCoA
{

  class RingDistrMPolyInlPPImpl: public SparsePolyRingBase
  {
  private:
    typedef DistrMPolyInlPP value_t; // DistrMPolyInlPP is the actual type of the values in a RingDistrMPolyInlPPImpl
    static value_t& import(RingElemRawPtr rawf);
    static const value_t& import(RingElemConstRawPtr rawf);

  public:
//    RingDistrMPolyInlPPImpl(const ring& R, const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    RingDistrMPolyInlPPImpl(const ring& R, const PPMonoid& PPM);
    ~RingDistrMPolyInlPPImpl();

  private: // Data members of RingDistrMPolyInlPPImpl
    const ring myCoeffRingValue;  ///< the coefficient ring
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
    std::string myImplDetails() const override  {return "RingDistrMPolyInlPP";}
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
    void myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< EXCEPTION SAFE
    bool myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< EXCEPTION SAFE
    //@}
    // MOVED TO SparsePolyRing

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

  inline RingDistrMPolyInlPPImpl::value_t& RingDistrMPolyInlPPImpl::import(RingElemRawPtr rawf)
  {
    return *static_cast<value_t*>(rawf.myRawPtr());
  }

  inline const RingDistrMPolyInlPPImpl::value_t& RingDistrMPolyInlPPImpl::import(RingElemConstRawPtr rawf)
  {
    return *static_cast<const value_t*>(rawf.myRawPtr());
  }


  RingDistrMPolyInlPPImpl::RingDistrMPolyInlPPImpl(const ring& R, const PPMonoid& PPM):
    myCoeffRingValue(R),
    myPPMValue(PPM),
    myOrdvArith(NewOrdvArith(ordering(PPM))),
    myDMPPool(sizeof(DistrMPolyInlPP), "RingDistrMPolyInlPPImpl::myDMPPool"),
    myNumIndetsValue(NumIndets(PPM)),
    mySummandPool(DistrMPolyInlPP::SummandSize(R, myOrdvArith), "RingDistrMPolyInlPPImpl::mySummandPool")
  {
    //    std::cout << "------RingDistrMPolyInlPPImpl:NewOrdvArith-called" << std::endl;
    CoCoA_ASSERT(IsPPMonoidOv(PPM));
    CoCoA_ASSERT(IsCommutative(myCoeffRingValue));
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


  RingDistrMPolyInlPPImpl::~RingDistrMPolyInlPPImpl()
  {}


  //----------------------------------------------------------------------
  // Functions which every ring must implement:
  //----------------------------------------------------------------------

  RingElemRawPtr RingDistrMPolyInlPPImpl::myNew() const
  {
    void* ptr = myDMPPool.alloc();
    new(ptr) DistrMPolyInlPP(myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool); // placement new
    return RingElemRawPtr(ptr);
  }


  RingElemRawPtr RingDistrMPolyInlPPImpl::myNew(const MachineInt& n) const
  {
    if (IsZero(n)) return myNew();
    unique_ptr<DistrMPolyInlPP> ans(new(myDMPPool.alloc()) DistrMPolyInlPP(myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool)); // placement new
    *ans = n;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDistrMPolyInlPPImpl::myNew(const BigInt& N) const
  {
    if (N == 0) return myNew();
    unique_ptr<DistrMPolyInlPP> ans(new(myDMPPool.alloc()) DistrMPolyInlPP(myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool)); // placement new
    *ans = N;
    return RingElemRawPtr(ans.release());
  }


  RingElemRawPtr RingDistrMPolyInlPPImpl::myNew(ConstRawPtr rawcopy) const
  {
    unique_ptr<DistrMPolyInlPP> ans(new(myDMPPool.alloc()) DistrMPolyInlPP(myCoeffRingValue, myPPMValue, myOrdvArith, mySummandPool)); // placement new
    *ans = import(rawcopy);
    return RingElemRawPtr(ans.release());
  }


  void RingDistrMPolyInlPPImpl::myDelete(RawPtr rawx) const
  {
    import(rawx).~DistrMPolyInlPP();
    myDMPPool.free(rawx.myRawPtr());
  }


  void RingDistrMPolyInlPPImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    swap(import(rawx), import(rawy));
  }


  void RingDistrMPolyInlPPImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
  }


  void RingDistrMPolyInlPPImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    import(rawlhs) = n;
  }


  void RingDistrMPolyInlPPImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    import(rawlhs) = N;
  }


  void RingDistrMPolyInlPPImpl::myAssignZero(RawPtr rawlhs) const
  {
    import(rawlhs).myAssignZero();
  }


  void RingDistrMPolyInlPPImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myCoeffRingValue));
    RingElem tmp(myCoeffRingValue);
    myCoeffRingValue->myRecvTwinFloat(raw(tmp), rawx);
    myCoeffEmbeddingHomCtor()->myApply(rawlhs, raw(tmp));
  }


  void RingDistrMPolyInlPPImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    import(rawlhs) = import(rawx);
    import(rawlhs).myNegate();
  }


  void RingDistrMPolyInlPPImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    add(import(rawlhs), import(rawx), import(rawy));
  }


  void RingDistrMPolyInlPPImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    sub(import(rawlhs), import(rawx), import(rawy));
  }


  bool RingDistrMPolyInlPPImpl::myIsZero(ConstRawPtr rawx) const
  {
    return IsZero(import(rawx));
  }


  bool RingDistrMPolyInlPPImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return IsEqual(import(rawx), import(rawy));
  }


  //----------------------------------------------------------------------
  // Functions which every PolyRing must implement
  //----------------------------------------------------------------------

  void RingDistrMPolyInlPPImpl::myIndetPower(RawPtr rawf, long var, long exp) const
  {
    CoCoA_ASSERT(var < myNumIndets());
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(static_cast<unsigned long>(exp) <= numeric_limits<SmallExponent_t>::max());
    RingElem ans(ring(this));
    vector<long> expv(myNumIndets()); // wasteful new/delete
    expv[var] = exp;
    import(raw(ans)).myPushFront(raw(one(myCoeffRingValue)), expv);
    mySwap(raw(ans), rawf); // do it this way to be exception clean
  }


  long RingDistrMPolyInlPPImpl::myNumTerms(ConstRawPtr rawx) const
  { return NumTerms(import(rawx)); }

  bool RingDistrMPolyInlPPImpl::myIsMonomial(ConstRawPtr rawf) const
  { return IsMonomial(import(rawf)); }


  RingElemAlias RingDistrMPolyInlPPImpl::myLC(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return LC(import(rawf));
  }


  void RingDistrMPolyInlPPImpl::myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const // EXCEPTION SAFE
  { import(rawf).myMulByCoeff(rawc); }


  bool RingDistrMPolyInlPPImpl::myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const // EXCEPTION SAFE
  { return import(rawf).myDivByCoeff(rawc); }


  //----------------------------------------------------------------------
  // Functions which every SparsePolyRing must implement:
  //----------------------------------------------------------------------


  RingElem RingDistrMPolyInlPPImpl::myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(!myCoeffRing()->myIsZero(rawc));
    RingElem ans(ring(this));
    myPushFront(raw(ans), rawc, rawpp);
    return ans;
  }


  SparsePolyIter RingDistrMPolyInlPPImpl::myBeginIter(ConstRawPtr rawf) const
  { return SparsePolyIter(new DistrMPolyInlPP::iter(import(rawf))); }

  SparsePolyIter RingDistrMPolyInlPPImpl::myEndIter(ConstRawPtr rawf) const
  { return SparsePolyIter(new DistrMPolyInlPP::iter(import(rawf), 0)); }


  void RingDistrMPolyInlPPImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    import(rawf).myPushFront(rawc, expv);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyInlPPImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    import(rawf).myPushBack(rawc, expv);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyInlPPImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    import(rawf).myPushFront(rawc, rawpp);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  void RingDistrMPolyInlPPImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    import(rawf).myPushBack(rawc, rawpp);
    CoCoA_ASSERT(myIsValid(rawf));
  }


  ConstRefPPMonoidElem RingDistrMPolyInlPPImpl::myLPP(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return LPP(import(rawf));
  }


  void RingDistrMPolyInlPPImpl::myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const
  { import(rawf).myMulByPP(rawpp); }

  bool RingDistrMPolyInlPPImpl::myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const
  { return IsZeroAddLCs(import(rawf), import(rawg)); }


  void RingDistrMPolyInlPPImpl::myMoveLMToFront(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    MoveLMToFront(import(rawf), import(rawg));
  }


  void RingDistrMPolyInlPPImpl::myMoveLMToBack(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    MoveLMToBack(import(rawf), import(rawg));
  }


  void RingDistrMPolyInlPPImpl::myDeleteLM(RawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    import(rawf).myDeleteLM();
  }


void RingDistrMPolyInlPPImpl::myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    DivLM(import(rawlhs), import(rawf), import(rawg));
  }


  int  RingDistrMPolyInlPPImpl::myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    return CmpLPP(import(rawf), import(rawg));
  }


  void RingDistrMPolyInlPPImpl::myAddClear(RawPtr rawf, RawPtr rawg) const
  { import(rawf).myAddClear(import(rawg)); }


  void RingDistrMPolyInlPPImpl::myAppendClear(RawPtr rawf, RawPtr rawg) const
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


  void RingDistrMPolyInlPPImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const //???? delete me???
  { import(rawf).myAddMulLM(import(rawh), import(rawg), /* SkipLMg = */ false); }

  void RingDistrMPolyInlPPImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag skip) const //???? delete me???
  { import(rawf).myAddMulLM(import(rawh), import(rawg), skip==SkipLMg); }

void RingDistrMPolyInlPPImpl::myReductionStep(RawPtr rawf, ConstRawPtr rawg) const
  { import(rawf).myReductionStep(import(rawg)); }

  void RingDistrMPolyInlPPImpl::myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& fscale) const
  { import(rawf).myReductionStepGCD(import(rawg), fscale); }


  //----------------------------------------------------------------------
  // Pseudo-ctors for (sparse) polynomial rings.

  namespace  {  // anonymous
    void CheckCoeffRing(const ring& K)
    {
      if (!IsCommutative(K))
        CoCoA_THROW_ERROR(ERR::NotCommutative, "NewPolyRing_DMPI pseudo ctor");
    }
    

    void CheckIndets(const ring& K, const std::vector<symbol>& IndetNames)
    {
      if (IndetNames.empty())
        CoCoA_THROW_ERROR(ERR::ReqNonEmpty, "NewPolyRing_DMPI pseudo ctor");
      if (!AreGoodIndetNames(K, IndetNames))
        CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPolyRing_DMPI pseudo ctor");
    }


    void CheckNumIndets(long n)
    {
      if (n <=0 )
        CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPolyRing_DMPI pseudo ctor");
    }


    void CheckNumIndets(long n, const std::vector<symbol>& IndetNames)
    {
      CheckNumIndets(n);
      if (n != len(IndetNames))
        CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPolyRing_DMPI pseudo ctor");
    }
    
  }  // anonymous namespace


  SparsePolyRing NewPolyRing_DMPI(const ring& CoeffRing, const PPMonoid& PPM)
  {
    //    std::cout << "------NewPolyRing_DMPI PPM" << std::endl;
    CheckCoeffRing(CoeffRing);
    CheckIndets(CoeffRing, symbols(PPM));
    if (!IsPPMonoidOv(PPM))
      CoCoA_THROW_ERROR(ERR::BadArg, "NewPolyRing_DMPI(k, PPM) -- PPM is not of Ov type");
    return SparsePolyRing(new RingDistrMPolyInlPPImpl(CoeffRing, PPM));
  }


  SparsePolyRing NewPolyRing_DMPI(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
    //    std::cout << "------NewPolyRing_DMPI ord" << std::endl;
    CheckNumIndets(NumIndets(ord), IndetNames);
    return NewPolyRing_DMPI(CoeffRing, NewPPMonoidOv(IndetNames, ord));
  }


  SparsePolyRing NewPolyRing_DMPI(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  {
    //    std::cout << "------NewPolyRing_DMPI OrdCtor" << std::endl;
    CheckCoeffRing(CoeffRing);
    CheckIndets(CoeffRing, IndetNames);
    return SparsePolyRing(new RingDistrMPolyInlPPImpl(CoeffRing, NewPPMonoidOv(IndetNames, OrdCtor(len(IndetNames)))));
  }


  SparsePolyRing NewPolyRing_DMPI(const ring& CoeffRing, const std::vector<symbol>& IndetNames)
  {
    return NewPolyRing_DMPI(CoeffRing, IndetNames, StdDegRevLex(len(IndetNames)));
  }


  SparsePolyRing NewPolyRing_DMPI(const ring& CoeffRing, long NumIndets)
  {
    CheckNumIndets(NumIndets);
    return NewPolyRing_DMPI(CoeffRing, SymbolRange("x",0,NumIndets-1), StdDegRevLex(NumIndets));
  }


  SparsePolyRing NewPolyRing_DMPI(const ring& CoeffRing, long NumIndets, const PPOrderingCtor& OrdCtor)
  {
    CheckNumIndets(NumIndets);
    return NewPolyRing_DMPI(CoeffRing, SymbolRange("x",0,NumIndets-1), OrdCtor(NumIndets));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/RingDistrMPolyInlPP.C,v 1.52 2024/07/15 17:10:49 bigatti Exp $
// $Log: RingDistrMPolyInlPP.C,v $
// Revision 1.52  2024/07/15 17:10:49  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.51  2022/02/18 14:11:57  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.50  2021/10/29 19:47:58  abbott
// Summary: Added keyword override (redmine 1625)
//
// Revision 1.49  2020/06/17 15:49:26  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.48  2020/01/09 16:41:09  bigatti
// -- improved error message for empty list of IndetNames in pseudo ctor
//
// Revision 1.47  2019/03/04 10:39:07  abbott
// Summary: Changed auto_ptr into unique_ptr
//
// Revision 1.46  2018/06/13 15:07:08  abbott
// Summary: Added check for empty list of symbols in pseudo-cto for DMPI and DMPII
//
// Revision 1.45  2018/05/18 12:15:56  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.44  2017/12/01 17:27:24  bigatti
// -- updated Copyright line
// -- removed doxygen initial comment
// -- some commented out debugging info
//
// Revision 1.43  2017/11/10 16:02:27  abbott
// Summary: Removed NewLexOrdering, NewStdDegLexOrdering, NewStdDegRevLexOrdering; consequential changes
//
// Revision 1.42  2015/12/01 13:11:01  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.41  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.40  2014/07/28 15:48:07  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble); removed lots of cruft
// Author: JAA
//
// Revision 1.39  2014/07/11 15:46:14  bigatti
// -- removed myOutputSelf (default impl) and added myImplDetails
//
// Revision 1.38  2014/07/09 11:43:13  abbott
// Summary: Removed AsPolyRing from commented out code
// Author: JAA
//
// Revision 1.37  2014/07/04 13:35:42  bigatti
// -- in printing of ring, indeterminates are now more compact
//
// Revision 1.36  2014/07/04 13:08:08  bigatti
// -- RingID into RingWithID
//
// Revision 1.35  2014/07/02 16:33:13  bigatti
// -- new way of printing rings with ID
//
// Revision 1.34  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.33  2014/04/30 16:11:22  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.32  2014/04/15 14:25:02  abbott
// Summary: Removed cruft (commented out code which is impl'd in PolyRing)
// Author: JAA
//
// Revision 1.31  2014/01/30 18:08:33  abbott
// Summary: Removed commented out cruft
// Author: JAA
//
// Revision 1.30  2014/01/30 17:29:22  abbott
// Summary: Added new pseudo ctor which accepts a PPM (only of PPMonoidOv type)
// Author: JAA
//
// Revision 1.29  2013/03/08 09:24:50  abbott
// Removed an unnecessary cast, and added a missing cast (inside a CoCoA_ASSERT).
//
// Revision 1.28  2012/10/24 12:19:56  abbott
// Changed return type of myLC.
//
// Revision 1.27  2012/10/17 09:40:16  abbott
// Replaced  RefRingElem  by  RingElem&
// (plus a few consequential changes)
//
// Revision 1.26  2012/10/11 14:27:59  abbott
// Removed "semantically risky" function RefLC/RawLC from DistrMPoly*.
// Reimplemented myRecvTwinFloat in RingDistrMPoly* more cleanly (but
// new impl does make a wasteful copy of the coeff).
//
// Revision 1.25  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.24  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.23  2011/11/09 14:10:49  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.22  2011/08/14 15:52:16  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.21  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.20  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.19  2010/10/08 14:18:58  bigatti
// -- changed printing style: RingDistrMPolyXXImpl(..) -->  RingDistrMPolyXX(..)
//
// Revision 1.18  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.17  2010/03/05 18:43:48  abbott
// Added pseudo-ctors allowing polynomial rings to be created specifying
// the ordering using a PPOrderingCtor object.
//
// Revision 1.16  2009/12/23 18:53:52  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.15  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.14  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.13  2009/09/28 16:19:43  bigatti
// -- unique implementation for myDeriv
//
// Revision 1.12  2009/09/25 13:02:09  bigatti
// -- myDiv with one implementation in SparsePolyRing
//
// Revision 1.11  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.10  2008/04/10 15:12:59  bigatti
// -- added myPushBack/Front(RawPtr, ConstRawPtr, PPMonoidElemConstRawPtr)
// -- minor tidying
//
// Revision 1.9  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.8  2007/10/11 16:31:42  bigatti
// -- changed  #ifdef CoCoA_ASSERT  into  CoCoA_DEBUG
//
// Revision 1.7  2007/05/31 16:34:37  bigatti
// -- Changed IsValid (now returns true of false and does not throw an error)
// -- using IsValid for sanity check in PushBack
//
// Revision 1.6  2007/05/31 15:56:32  bigatti
// -- default implementation for IamCommutative, IamIntegralDomain,
//    IamGCDDomain, IamField, myCharacteristic  in PolyRing
// -- default implementation for mySymbols  in SparsePolyRing
// -- cleaned up pseudo-ctors
//
// Revision 1.4  2007/05/21 14:50:56  bigatti
// -- myPushFront and myPushBack now accept zero coefficient
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
// Revision 1.23  2007/03/08 17:43:10  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.22  2007/03/08 16:55:06  cocoa
// Changed name of "range" function to "SymbolRange".
//
// Revision 1.21  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.20  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.19  2007/01/20 14:07:25  bigatti
// -- moved code for homomorphism into common implementation in SparsePolyRing
//
// Revision 1.18  2007/01/15 15:47:57  cocoa
// -- added prefix "raw" to RawPtr arguments names
// -- changed rhs into rawx, n, or N
//
// Revision 1.17  2007/01/13 14:14:34  cocoa
// Overhaul of RingHom code: it nows uses SmartPtrIRC, and printing is more logical.
// Have not yet updated the documentation.
//
// Revision 1.16  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.15  2006/12/06 17:33:03  cocoa
// -- added #include "config.H"
//
// Revision 1.14  2006/11/22 17:51:31  cocoa
// -- moved printing functions into unified implementation in SparsePolyRing
//
// Revision 1.13  2006/11/21 18:09:23  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.12  2006/11/14 17:18:47  cocoa
// -- added comment about myRefCountZero()
//
// Revision 1.11  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPoly
//
// Revision 1.10  2006/11/08 16:21:59  cocoa
// Structural cleaning of RingHom; many consequential changes.
//
// Revision 1.9  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.8  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.7  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.6  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
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
// Revision 1.24  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.23  2006/05/29 13:20:33  cocoa
// -- added implementation of "intersect" for SparsePolyRing
//
// Revision 1.22  2006/05/12 17:06:09  cocoa
// -- moved myIsUnit, myGcd into SparsePolyRing (common implementation)
//
// Revision 1.21  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.20  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.19  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.18  2006/04/21 16:47:06  cocoa
// -- new syntax for ComputeGBasis by Max
//
// Revision 1.17  2006/03/30 15:56:56  cocoa
// -- changed: unsigned int --> SmallExponent_t (1 occurrence)
//
// Revision 1.16  2006/03/21 17:55:08  cocoa
// -- changed: "my.." coding conventions
//
// Revision 1.15  2006/03/21 09:43:13  cocoa
// Changed names of some member fns of ideals (dealing with setting and testing
// the flags for primeness and maximality).  Hope icc will complain less now.
//
// Revision 1.14  2006/03/17 18:10:27  cocoa
// -- changed: myMul --> myMulByPP
//
// Revision 1.13  2006/03/16 17:56:14  cocoa
// -- changed: mul, div --> myMulByCoeff, myDivByCoeff (+check IsOne(c))
// -- changed: myOutputSelf for homomorphism
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
// -- fixed: PPMonoidElem LPP(f) now returns ConstRefPPMonoidElem
//
// Revision 1.8  2006/02/20 22:41:19  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.7  2006/02/14 16:19:37  cocoa
// -- defined "IdealImpl::contains"
//
// Revision 1.6  2006/02/13 13:56:50  cocoa
// -- changed: "GReductor" --> "ComputeGBasis"
//
// Revision 1.5  2006/02/13 13:52:40  cocoa
// Fixed same bug as was in RingDistrMPoly.
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
// Revision 1.7  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.6  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.5  2005/07/15 16:34:33  cocoa
// Added iterators for sparse polynomials.
// The code compiles (and the old tests still run).
// It'd Friday evening -- I'm going home before
// getting any ideas about making the iterator code run.
//
// Revision 1.4  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
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
// Revision 1.23  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.22  2004/11/18 18:33:41  cocoa
// Now every ring know its own "one" element (as well as "zero").
// Several consequential changes.
//
// Revision 1.21  2004/11/11 14:33:38  cocoa
// -- minor changes for doxygen
//
// Revision 1.20  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.19  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.18  2004/10/29 16:13:29  cocoa
// -- new default monoid: NewPPMonoidSafe with NewStdDegRevLexOrdering
//
// Revision 1.17  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.16  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.15  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.14  2004/06/30 16:46:06  cocoa
// End of day check-in: mostly tidying up, and correcting some
// overly lax access permissions.
//
// Revision 1.13  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.12  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.11  2004/04/08 15:33:34  cocoa
// Added function IsInteger, and the related RingBase::myIsInteger
// virtual function, plus all necessary implementations.
//
// Revision 1.10  2004/03/20 17:46:10  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.9  2004/03/04 11:23:38  cocoa
// - added argument in GReductor constructor:  Reductors::DontUseBorel
//
// Revision 1.8  2004/01/30 14:07:10  cocoa
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
// Revision 1.7  2004/01/28 15:58:04  cocoa
// Added myCmpDeg member function.
//
// Revision 1.6  2003/11/14 13:06:05  cocoa
// -- New function "myIsPrintAtom" for printing polynomials and fractions
//
// Revision 1.5  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.4  2003/10/09 14:55:20  cocoa
// - minor debugging after merge
//
// Revision 1.3  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:31  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.7  2003/06/23 17:01:55  abbott
// Minor cleaning prior to public release.
// Name changes, and consequential changes.
//
// Revision 1.6  2003/05/30 12:03:32  abbott
// Fairly extensive changes: completed alignment with abstract class
// RingBase, including homomorphisms and ideals.
//
// Revision 1.5  2002/11/18 18:03:31  abbott
// COMPLETELY NEW VERSION.  All the real code has been moved to a new
// class (called DMPI), so RingDistrMPolyInlPP is just "a nice front end".
// This version appears to be a bit faster than older code, and simpler
// too.
//
// Revision 1.4  2002/11/18 17:54:37  abbott
// Checking in only to have a complete archive.  I think this code ran once.
// It will be replaced immediately by a totally different new version.
//
// Revision 1.3  2002/07/05 15:35:24  abbott
// Fixed InitSummand so that it (correctly) zeros the module position.
// Changed name of LexCmp to CmpOrdvs (reflecting the change in PPOrdering.H).
// Added definitions of member functions power and IsDivisible.
//
// Revision 1.2  2002/06/26 17:16:06  bigatti
// - corrected MoveLT
// - changed: MoveLT --> MoveLM
// - changed: some calls to MakeWritable --> ASSERT(myRefCount==1)
// - added: CmpLPP(f,g);  mul(f,pp);  reduce(f,g);  AddClear(f,g);
//
// Revision 1.1  2002/06/22 17:05:47  abbott
// Initial revision
//
//
