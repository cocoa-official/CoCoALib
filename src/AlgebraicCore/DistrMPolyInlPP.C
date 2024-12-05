//   Copyright (c)  2005-2012,2017  John Abbott, Anna M. Bigatti

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


// Source code for class DistMPolyInlPP

#include "CoCoA/DistrMPolyInlPP.H"

#include "CoCoA/BigInt.H"  // uses binomial
#include "CoCoA/error.H"


//using std::swap;
#include <iostream>
//using "<<"
//#include <vector>
using std::vector;


namespace CoCoA
{

  bool IsCompatible(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g)
  {
    return (f.myCoeffRing == g.myCoeffRing) &&
           (f.myPPM == g.myPPM) &&
           (&f.mySummandMemory == &g.mySummandMemory);
  }


  void DistrMPolyInlPP::ourSwap(DistrMPolyInlPP& f, DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    std::swap(f.mySummands, g.mySummands);
    std::swap(f.myLastSummand, g.myLastSummand);
  }


  inline DistrMPolyInlPP::AutoPtrSummand::AutoPtrSummand(const DistrMPolyInlPP& f):
      myPtr(nullptr),
      myMemMgr(f.mySummandMemory),
      myR(CoeffRing(f)),
      myOrdvArith(f.myOrdvArith)
  {}

  inline DistrMPolyInlPP::AutoPtrSummand::~AutoPtrSummand() noexcept
  {
    if (myPtr == nullptr) return;
    myR->myDelete(myPtr->myCoeff);
    myMemMgr.free(myPtr);
  }


  inline DistrMPolyInlPP::summand& DistrMPolyInlPP::AutoPtrSummand::operator*() const noexcept
  {
    CoCoA_ASSERT(myPtr != nullptr);
    return *myPtr;
  }


  inline DistrMPolyInlPP::summand* DistrMPolyInlPP::AutoPtrSummand::operator->() noexcept
  {
    CoCoA_ASSERT(myPtr != nullptr);
    return myPtr;
  }


  inline DistrMPolyInlPP::summand* DistrMPolyInlPP::AutoPtrSummand::get() const noexcept
  {
    // Deliberately do not require myPtr to be non NULL
    return myPtr;
  }


  inline DistrMPolyInlPP::summand* DistrMPolyInlPP::AutoPtrSummand::release() noexcept
  {
    CoCoA_ASSERT(myPtr != nullptr);
    summand* ans = myPtr;
    myPtr = nullptr;
    return ans;
  }


  // BUG: THIS IS NOT EXCEPTION CLEAN!!!!!
  inline void DistrMPolyInlPP::AutoPtrSummand::myRenew()
  {
    CoCoA_ASSERT(myPtr == nullptr);
    myPtr = static_cast<summand*>(myMemMgr.alloc());
    myPtr->myNext = nullptr;
    myPtr->myCoeff = myR->myNew(); // could throw
    myOrdvArith->myAssignZero(myPtr->myOrdv);
  }


  void DistrMPolyInlPP::ourDeleteSummands(DistrMPolyInlPP::summand* ptr, const ring& R, MemPool& MemMgr)
  {
    DistrMPolyInlPP::summand* next;
    while (ptr != nullptr)
    {
      next = ptr->myNext;
      R->myDelete(ptr->myCoeff);
      MemMgr.free(ptr);
      ptr = next;
    }
  }


  DistrMPolyInlPP::DistrMPolyInlPP(const ring& R, const PPMonoid& PPM, const OrdvArith::reference& OA, MemPool& MemMgr):
      myCoeffRing(R),
      myPPM(PPM),
      myOrdvArith(OA),
      mySummandMemory(MemMgr)
  {
    mySummands = nullptr;
    myLastSummand = nullptr;
  }


  DistrMPolyInlPP::~DistrMPolyInlPP()
  {
    ourDeleteSummands(mySummands, myCoeffRing, /*myPPM,*/ mySummandMemory);
  }


  DistrMPolyInlPP::DistrMPolyInlPP(const DistrMPolyInlPP& copy):
      myCoeffRing(copy.myCoeffRing),
      myPPM(copy.myPPM),
      myOrdvArith(copy.myOrdvArith),
      mySummandMemory(copy.mySummandMemory)
  {
    mySummands = nullptr;
    myLastSummand = nullptr;
    // !!! THIS LOOP IS NOT EXCEPTION CLEAN !!!
    for (summand* it = copy.mySummands; it != nullptr; it = it->myNext)
      myPushBack(myCopySummand(it));
  }


  DistrMPolyInlPP& DistrMPolyInlPP::operator=(const DistrMPolyInlPP& rhs)
  {
    if (this == &rhs) return *this;
    DistrMPolyInlPP copy(rhs);
    ourSwap(*this, copy);
    return *this;
  }


  DistrMPolyInlPP& DistrMPolyInlPP::operator=(const MachineInt& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this;
    AutoPtrSummand t(*this);
    t.myRenew();
    myCoeffRing->myAssign(t->myCoeff, rhs);
    if (!myCoeffRing->myIsZero(t->myCoeff))
    {
      mySummands = t.release();
      myLastSummand = mySummands;
    }
    return *this;
  }


  DistrMPolyInlPP& DistrMPolyInlPP::operator=(const BigInt& rhs)
  {
    myAssignZero();
    AutoPtrSummand t(*this);
    t.myRenew();
    myCoeffRing->myAssign(t->myCoeff, rhs);
    if (!myCoeffRing->myIsZero(t->myCoeff))
    {
      mySummands = t.release();
      myLastSummand = mySummands;
    }
    return *this;
  }

  DistrMPolyInlPP& DistrMPolyInlPP::operator=(const BigRat& rhs)
  {
    myAssignZero();
    AutoPtrSummand t(*this);
    t.myRenew();
    myCoeffRing->myAssign(t->myCoeff, rhs);
    if (!myCoeffRing->myIsZero(t->myCoeff))
    {
      mySummands = t.release();
      myLastSummand = mySummands;
    }
    return *this;
  }


//----------------------------------------------------------------------
// operations on summands


  DistrMPolyInlPP::summand* DistrMPolyInlPP::myCopySummand(const summand* original) const
  {
    AutoPtrSummand copy(*this);
    copy.myRenew();
    myCoeffRing->myAssign(copy->myCoeff, original->myCoeff);
    myOrdvArith->myAssign(copy->myOrdv, original->myOrdv);
    return copy.release();
  }


  void DistrMPolyInlPP::myAssignZero()
  {
    ourDeleteSummands(mySummands, myCoeffRing, /*myPPM,*/ mySummandMemory);
    mySummands = nullptr;
    myLastSummand = nullptr;
  }


  bool DistrMPolyInlPP::myIsEqual(const summand* const lhs, const summand* const rhs) const
  {
    return (myOrdvArith->myCmp(lhs->myOrdv, rhs->myOrdv) == 0 &&
            myCoeffRing->myIsEqual(lhs->myCoeff, rhs->myCoeff));
  }


  long NumTerms(const DistrMPolyInlPP& f)
  {
    long nsummands = 0;
    for (const DistrMPolyInlPP::summand* it = f.mySummands; it != nullptr; it = it->myNext)
      ++nsummands;
    return nsummands;
  }


  RingElemAlias LC(const DistrMPolyInlPP& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    return RingElemAlias(f.myCoeffRing, f.mySummands->myCoeff);
  }


  void MoveLMToFront(DistrMPolyInlPP& f, DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlPP::summand* ltg = g.mySummands;
    g.mySummands = g.mySummands->myNext;
    if (g.mySummands == nullptr) { g.myLastSummand = nullptr; }
    f.myPushFront(ltg);
  }


  void MoveLMToBack(DistrMPolyInlPP& f, DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlPP::summand* ltg = g.mySummands;
    g.mySummands = g.mySummands->myNext;
    if (g.mySummands == nullptr) { g.myLastSummand = nullptr; }
    ltg->myNext = nullptr;
    f.myPushBack(ltg);
  }


  void DistrMPolyInlPP::myDeleteLM()
  {
    CoCoA_ASSERT(!IsZero(*this));

    DistrMPolyInlPP::summand* old_lm = mySummands;
    mySummands = old_lm->myNext;
    if (mySummands == nullptr) { myLastSummand = nullptr; }
    old_lm->myNext = nullptr;
    ourDeleteSummands(old_lm, myCoeffRing, /*myPPM,*/ mySummandMemory);
  }


//   void wdeg(degree& d, const DistrMPolyInlPP& f)
//   {
//     CoCoA_ASSERT(!IsZero(f));
//     f.myOrdvArith->myWDeg(d, f.mySummands->myOrdv);
//   }


//   int CmpWDeg(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g)
//   {
//     CoCoA_ASSERT(!IsZero(f));
//     CoCoA_ASSERT(!IsZero(g));
//     return f.myOrdvArith->myCmpWDeg(f.mySummands->myOrdv, g.mySummands->myOrdv);
//   }


// This fn offers only the weak exception guarantee!!!
  void DistrMPolyInlPP::myAddMulSummand(const summand* s, const DistrMPolyInlPP& g, bool SkipLMg)  // this += s*g
  {
    CoCoA_ASSERT(IsCompatible(*this, g));

//???    const PPOrdering& ord = ordering(myPPM);
    const ring& R = myCoeffRing;

    const summand* g_smnd = g.mySummands;
    if (SkipLMg)    g_smnd = g_smnd->myNext;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;

    AutoPtrSummand tmp_smnd(*this);
    tmp_smnd.myRenew();
    RingElem tmp(myCoeffRing);
    int CMP = 0;

    //bool qIsOne = myOrdvArith->myIsZero(s->myOrdv);
    bool qIsOne = false;

    for (; f_smnd != nullptr && g_smnd != nullptr; g_smnd = g_smnd->myNext)
    {
      R->myMul(tmp_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff);
      if (R->myIsZero(tmp_smnd->myCoeff))
        continue; // skip this g_smnd as coeff mults to 0
      myOrdvArith->myMul(tmp_smnd->myOrdv, s->myOrdv, g_smnd->myOrdv);
      // tmp_smnd is now contains term to be added; it is non-zero!
      if (qIsOne)
      {
        while (f_smnd != nullptr && (CMP=myOrdvArith->myCmp(f_smnd->myOrdv, g_smnd->myOrdv)) >0)
          f_smnd = *(f_prev = &f_smnd->myNext);
        myOrdvArith->myAssign(tmp_smnd->myOrdv, g_smnd->myOrdv);
      }
      else
      {
        while (f_smnd != nullptr && (CMP=myOrdvArith->myCmp(f_smnd->myOrdv, tmp_smnd->myOrdv)) >0)
          f_smnd = *(f_prev = &f_smnd->myNext);
      }
      if (f_smnd == nullptr)
      {
//        if (!R->myIsZero(tmp_smnd->myCoeff))
//        {
          myPushBack(tmp_smnd.release());
          tmp_smnd.myRenew();
          g_smnd = g_smnd->myNext;
//        }
        break;
      }
      if (CMP == 0)
      {
//        if (R->myIsZeroAddMul(f_smnd->myCoeff, raw(tmp), s->myCoeff, g_smnd->myCoeff))
        R->myAdd(f_smnd->myCoeff, f_smnd->myCoeff, tmp_smnd->myCoeff);
        if (R->myIsZero(f_smnd->myCoeff))
          myRemoveSummand(f_prev);  // f_prev = f_prev;
        else
          f_prev = &f_smnd->myNext;
        f_smnd = *f_prev;
      }
      else // (CMP < 0)
      {
//        R->myMul(tmp_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff);
        myInsertSummand(tmp_smnd.release(), f_prev);
        tmp_smnd.myRenew();
        f_prev = &(*f_prev)->myNext;
        // f_smnd = f_smnd;
      }
    }
    for (; g_smnd != nullptr; g_smnd = g_smnd->myNext)
    {
      R->myMul(tmp_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff);
      if (R->myIsZero(tmp_smnd->myCoeff)) continue;
      myOrdvArith->myMul(tmp_smnd->myOrdv, s->myOrdv, g_smnd->myOrdv);
      myPushBack(tmp_smnd.release());
      tmp_smnd.myRenew();
    }

//      clog << "AddMul: produced f=";output(clog,f);clog<<std::endl;
//      clog << "------------------------------------------------------"<<std::endl;

  }


  void DistrMPolyInlPP::myAddMulLM(const DistrMPolyInlPP& h, const DistrMPolyInlPP& g, bool SkipLMg)
  {
    if ( IsZero(h) )  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    //    if ( !IsMonomial(h) )
    //    CoCoA_THROW_ERROR(ERR::NYI, "DistrMPolyInlPP::myAddMul first argument len>1");
    myAddMulSummand(h.mySummands, g, SkipLMg);
  }


//   void DistrMPolyInlPP::myWeylAddMulSummand(const summand* s, const DistrMPolyInlPP& g, bool SkipLMg)
//   {
//     CoCoA_ASSERT(IsCompatible(*this, g));

// //???    const PPOrdering& ord = ordering(myPPM);
//     const ring& R = myCoeffRing;
//     const size_t nvars = NumIndets(myPPM);
// //???    clog << "AddMul: Doing funny product of the following two polys" << std::endl;
// //???    output(clog, g);
// //???    CoCoA_ASSERT(myRefCount == 1);
//     //????    MakeWritable(f);

//     const summand* g_term = g.mySummands;
//     if (SkipLMg)    g_term = g_term->myNext;
// //???    summand** f_prev = &mySummands;
// //???    summand*  f_term = *f_prev;

//     DistrMPolyInlPP ppg = g;
//     vector<long> expv(nvars);
//     myOrdvArith->myComputeExpv(expv, s->myOrdv);
// //???    clog << "expv: "; for (int i=0; i<myNumIndets;++i) clog << expv[i] << "  "; clog << std::endl;
//     for (size_t indet = nvars/2; indet < nvars; ++indet)
//     {
//       long n = expv[indet];
//       if (n == 0) continue;
// //???      clog << "AddMul: doing D variable with index " << indet - myNumIndets/2 << std::endl;
//       DistrMPolyInlPP der = ppg;

// //???      ppg *= IndetPower(myPPM, indet, n);
//       ppg.myMulByPP(raw(IndetPower(myPPM, indet, n)));
// //      mul(raw(ppg), raw(ppg), raw(IndetPower(indet, n)));

//       for (long i=1; i <= n; ++i)
//       {
//         deriv(der, der, indet-nvars/2);
// //???        deriv(raw(der), raw(der), indet-nvars/2);
// //???        clog << "der(" << i << ")="; output(clog, raw(der)); clog << std::endl;

// //        ppg += binomial(n, i)*der*IndetPower(myPPM, indet, n-i); // *IndetPower(myPPM, h, 2*i); // for homog case
//         AutoPtrSummand jaa(*this);
//         jaa.myRenew();
//         R->myAssign(jaa->myCoeff, binomial(n, i));
//         myOrdvArith->myMulIndetPower(jaa->myOrdv, indet, n-i);
// // if (HOMOG_CASE)        ord->mul(jaa->myOrdv, jaa->myOrdv, IndetPower(myPPM, h, 2*i));
//         ppg.myAddMulSummand(jaa.get(), der, false);
//       }
//     }
//     { // f *= indet^deg(pp, indet); for the early vars
//       for (size_t indet = nvars/2; indet < nvars; ++indet)
//         expv[indet] = 0;
//       AutoPtrSummand jaa(*this);
//       jaa.myRenew();
//       R->myAssign(jaa->myCoeff, s->myCoeff);
// //???      ord->assign(jaa->myOrdv, expv????);
//       myOrdvArith->myAssignFromExpv(jaa->myOrdv, expv);
//       myAddMulSummand(jaa.get(), ppg, false);
//     }
//   }




  void DistrMPolyInlPP::myReductionStep(const DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(&g!=this);
    CoCoA_ASSERT(mySummands != nullptr);
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlPP tmp_poly(myCoeffRing, myPPM, myOrdvArith, mySummandMemory);

    DivLM(tmp_poly, *this, g);
    tmp_poly.myNegate();
    myDeleteLM();
    myAddMulLM(tmp_poly, g, /*SkipLMg = */ true );
  }


  static void ComputeFScaleAndGScale(ConstRefRingElem LCf,
                                     ConstRefRingElem LCg,
                                     RingElem& fscale,
                                     RingElem& gscale)
  {
//???    CoCoA_ASSERT(!R->myIsEqual(LCg,-1));
    RingElem gcdLC = gcd(LCf, LCg);
    gscale = -LCf/gcdLC;
    if (LCg == gcdLC)
    {
      fscale = 1;
      return;
    }
    fscale = LCg/gcdLC;
    if (!IsInvertible(fscale)) return;
    gscale /= fscale;
    fscale = 1;
  }


  void DistrMPolyInlPP::myReductionStepGCD(const DistrMPolyInlPP& g, RingElem& fscale)
  {
    CoCoA_ASSERT(&g!=this);
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlPP tmp_poly(myCoeffRing, myPPM, myOrdvArith, mySummandMemory);
    RingElem gscale(myCoeffRing);

    ComputeFScaleAndGScale(LC(*this), LC(g), fscale, gscale);

    if ( !IsOne(fscale) )    myMulByCoeff(raw(fscale));
    DivLM(tmp_poly, *this, g);
    tmp_poly.myNegate();
    myDeleteLM();
    myAddMulLM(tmp_poly, g, /*SkipLMg = */ true);
  }


  void DistrMPolyInlPP::myAddClear(DistrMPolyInlPP& g) // sets g to 0 as a side-effect
  {
    if (g.mySummands == nullptr) return;

    const ring& R = myCoeffRing;
//???    const PPOrdering& ord = myOrdvArith;
    typedef DistrMPolyInlPP::summand summand;

    summand*  g_smnd = g.mySummands;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;
    int CMP=0;

    if (mySummands != nullptr &&
        myOrdvArith->myCmp(myLastSummand->myOrdv, g.mySummands->myOrdv) <= 0)
    {
      while ( f_smnd != nullptr && g_smnd != nullptr )
      {
        while (f_smnd != nullptr &&
               (CMP = myOrdvArith->myCmp(f_smnd->myOrdv, g_smnd->myOrdv)) >0)
          f_smnd = *(f_prev = &f_smnd->myNext);
        if (f_smnd == nullptr)  break;
        g.mySummands = g.mySummands->myNext;
        g_smnd->myNext = nullptr;
        if (CMP == 0)
        {
          R->myAdd(f_smnd->myCoeff, f_smnd->myCoeff, g_smnd->myCoeff);
          if (R->myIsZero(f_smnd->myCoeff))
            myRemoveSummand(f_prev);
          ourDeleteSummands(g_smnd, myCoeffRing, /*myPPM,*/ mySummandMemory);
        }
        else // (CMP < 0)
        {
          myInsertSummand(g_smnd, f_prev);
          f_prev = &(*f_prev)->myNext;
        }
        f_smnd = *f_prev;
        g_smnd = g.mySummands;
      }
    }
    if (g.mySummands != nullptr)
    {
      if (myLastSummand == nullptr)
        mySummands = g.mySummands;
      else
        myLastSummand->myNext = g.mySummands;
      myLastSummand = g.myLastSummand;
      g.mySummands = nullptr;
    }
    g.myLastSummand = nullptr;
  }


  void DistrMPolyInlPP::myAppendClear(DistrMPolyInlPP& g)  noexcept // sets g to 0
  {
    if (g.mySummands == nullptr) return;
    if (myLastSummand == nullptr)
      mySummands = g.mySummands;
    else
      myLastSummand->myNext = g.mySummands;
    myLastSummand = g.myLastSummand;
    g.mySummands = nullptr;
    g.myLastSummand = nullptr;
  }


  ConstRefPPMonoidElem LPP(const DistrMPolyInlPP& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    return ConstRefPPMonoidElem(f.myPPM, PPMonoidElemConstRawPtr(f.mySummands->myOrdv));
  }


  void DivLM(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, const DistrMPolyInlPP& g) // lhs = LM(f)/LM(g)
  {
    CoCoA_ASSERT(IsCompatible(f, g) && IsCompatible(lhs, f));
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    //    clog << "DivLM" << std::endl;
    const ring& R = f.myCoeffRing;    // shorthand
    typedef DistrMPolyInlPP::summand summand;            // shorthand
    const summand* const LMf = f.mySummands;  // shorthand
    const summand* const LMg = g.mySummands;  // shorthand

    DistrMPolyInlPP::AutoPtrSummand SpareSummand(lhs);
    SpareSummand.myRenew();
    CoCoA_ASSERT(R->myIsDivisible(SpareSummand->myCoeff, LMf->myCoeff, LMg->myCoeff));
    R->myDiv(SpareSummand->myCoeff, LMf->myCoeff, LMg->myCoeff);
    f.myOrdvArith->myDiv(SpareSummand->myOrdv, LMf->myOrdv, LMg->myOrdv);
    lhs.myAssignZero();  // CANNOT myAssignZero() EARLIER in case lhs aliases f or g.
    lhs.myPushBack(SpareSummand.release());
  }



  // This is lengthier than you might expect because it is exception safe.
  // It also works fine if c aliases a coefficient of the polynomial.
  void DistrMPolyInlPP::myMulByCoeff(RingElemConstRawPtr rawc)
  {
    CoCoA_ASSERT(!myCoeffRing->myIsZero(rawc));
    if (myCoeffRing->myIsOne(rawc)) return;
    const long n = NumTerms(*this);
    vector<RingElem> NewCoeffs(n, zero(myCoeffRing));
    long k=0;
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
    {
      myCoeffRing->myMul(raw(NewCoeffs[k]), it->myCoeff, rawc);
      ++k;
    }
    // The new coeffs are in NewCoeffs now swap them into the poly,
    // but if a coeff is 0 (!zero-divisors!), delete the corresponding term.
    summand** SummandPtr = &mySummands;
    summand* it = mySummands;
    for (long k=0; k < n; ++k)
    {
      CoCoA_ASSERT(it != nullptr);
      CoCoA_ASSERT(*SummandPtr == it);
      if (!IsZero(NewCoeffs[k]))
      {
        myCoeffRing->mySwap(raw(NewCoeffs[k]), it->myCoeff);
        SummandPtr = &it->myNext;
        it = it->myNext;
        continue;
      }
      summand* next = it->myNext;
      myRemoveSummand(SummandPtr);
      it = next;
    }
  }


  // This is lengthier than you might expect because it is exception safe.
  // It also works fine if c aliases a coefficient of the polynomial.
  bool DistrMPolyInlPP::myDivByCoeff(RingElemConstRawPtr rawc)
  {
    CoCoA_ASSERT(!myCoeffRing->myIsZero(rawc));
    if (myCoeffRing->myIsOne(rawc)) return true;
    vector<RingElem> NewCoeffs(NumTerms(*this), zero(myCoeffRing));
    long n=0;
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
    {
      if (!myCoeffRing->myIsDivisible(raw(NewCoeffs[n]), it->myCoeff, rawc))
        return false;
      ++n;
    }
    // The new coeffs are in NewCoeffs now swap them into the poly.
    n = 0;
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
    {
      myCoeffRing->mySwap(raw(NewCoeffs[n]), it->myCoeff);
      ++n;
    }
    return true;
  }


  void DistrMPolyInlPP::myMulByPP(PPMonoidElemConstRawPtr rawpp)
  {
    vector<OrdvArith::OrdvElem> ordv(OrdvWords(myOrdvArith));
    {
      vector<long> expv(NumIndets(myPPM));
      myPPM->myExponents(expv, rawpp);
      myOrdvArith->myAssignFromExpv(&ordv[0], expv);
    }

    for (summand* f_smnd = mySummands; f_smnd != nullptr; f_smnd = f_smnd->myNext)
      myOrdvArith->myMul(f_smnd->myOrdv, f_smnd->myOrdv, &ordv[0]);
  }

// ??? WANT DIVISION BY A PP TOO???


//   void DistrMPolyInlPP::myWeylMul(PPMonoidElemConstRawPtr rawpp)
//   {
//     vector<OrdvArith::OrdvElem> ordv(OrdvWords(myOrdvArith));
//     {
//       vector<long> expv(NumIndets(myPPM));
//       myPPM->myExponents(expv, rawpp);
//       myOrdvArith->myAssignFromExpv(&ordv[0], expv);
//     }

//     for (summand* f_smnd = mySummands; f_smnd != nullptr; f_smnd = f_smnd->myNext)
//       myOrdvArith->myMul(f_smnd->myOrdv, f_smnd->myOrdv, &ordv[0]);
//   }


  void DistrMPolyInlPP::myPushFront(RingElemConstRawPtr rawc, const std::vector<long>& expv)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    AutoPtrSummand tmp(*this);
    tmp.myRenew();
    myCoeffRing->myAssign(tmp->myCoeff, rawc);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushFront(tmp.release()); 
  }


  void DistrMPolyInlPP::myPushBack(RingElemConstRawPtr rawc, const std::vector<long>& expv)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    AutoPtrSummand tmp(*this);
    tmp.myRenew();
    myCoeffRing->myAssign(tmp->myCoeff, rawc);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushBack(tmp.release());
  }


  void DistrMPolyInlPP::myPushFront(RingElemConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    AutoPtrSummand tmp(*this);
    tmp.myRenew();
    myCoeffRing->myAssign(tmp->myCoeff, rawc);
    //    myOrdvArith->myAssign(tmp->myOrdv, rawpp);
    vector<long> expv(NumIndets(myPPM));
    myPPM->myExponents(expv, rawpp);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushFront(tmp.release()); 
  }


  void DistrMPolyInlPP::myPushBack(RingElemConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    AutoPtrSummand tmp(*this);
    tmp.myRenew();
    myCoeffRing->myAssign(tmp->myCoeff, rawc);
    //    myOrdvArith->myAssign(tmp->myOrdv, rawpp);
    vector<long> expv(NumIndets(myPPM));
    myPPM->myExponents(expv, rawpp);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushBack(tmp.release());
  }


  void DistrMPolyInlPP::myPushFront(summand* t)
  {
    // CoCoA_ASSERT(t->myNext == nullptr);  // deal with t->myNext elsewhere
    CoCoA_ASSERT(IsZero(*this) || myOrdvArith->myCmp(t->myOrdv, mySummands->myOrdv) > 0);
    t->myNext = mySummands;
    mySummands = t;
    if (myLastSummand == nullptr) myLastSummand = t;
  }


  void DistrMPolyInlPP::myPushBack(summand* t)
  {
    CoCoA_ASSERT(t->myNext == nullptr);
    // Cannot easily check that t really is smaller than smallest PP in the polynomial.
    if (mySummands == nullptr) { myPushFront(t); return; }
    myLastSummand->myNext = t;
    myLastSummand = t;
  }


  // Delete summand *prev_link
  void DistrMPolyInlPP::myRemoveSummand(summand** prev_link)
  {

    summand* DeleteMe = *prev_link;
    CoCoA_ASSERT(DeleteMe != nullptr);
    if (DeleteMe->myNext == nullptr)
    {
      if (*prev_link == mySummands) myLastSummand = nullptr;
      else
      {
        for (myLastSummand = mySummands; myLastSummand->myNext != *prev_link; myLastSummand = myLastSummand->myNext) {}
      }
    }

    *prev_link = DeleteMe->myNext; // updates mySummands if deleting 1st summand
    DeleteMe->myNext = nullptr;
    ourDeleteSummands(DeleteMe, myCoeffRing, /*myPPM,*/ mySummandMemory);
  }


  // insert summand s right before summand *prev_link;
  // to insert after last term **prev_link is the myNext field of last summand.
  void DistrMPolyInlPP::myInsertSummand(summand* s, summand** prev_link)
  {
    // Need to update myLastSummand in 2 cases: poly was 0, or prev_link is myNext of last summand
    if (myLastSummand == nullptr ||
        prev_link == &(myLastSummand->myNext))
      myLastSummand = s;
    s->myNext = (*prev_link);
    (*prev_link) = s;
  }


  bool IsZeroAddLCs(DistrMPolyInlPP& f, DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    CoCoA_ASSERT( LPP(f) == LPP(g) );
    f.myCoeffRing->myAdd(f.mySummands->myCoeff, f.mySummands->myCoeff, g.mySummands->myCoeff);
    g.myDeleteLM();
    if (!f.myCoeffRing->myIsZero(f.mySummands->myCoeff)) return false;
    f.myDeleteLM();
    return true;
  }


  void DistrMPolyInlPP::myNegate()
  {
    for (summand* iter = mySummands; iter != nullptr; iter = iter->myNext)
      myCoeffRing->myNegate(iter->myCoeff, iter->myCoeff);
  }


  void add(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h)
  {
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(g, h));
    const ring& R = lhs.myCoeffRing;
    DistrMPolyInlPP ans(lhs.myCoeffRing, lhs.myPPM, lhs.myOrdvArith, lhs.mySummandMemory);

    typedef DistrMPolyInlPP::summand summand;
    const summand* gterm = g.mySummands;
    const summand* hterm = h.mySummands;

    if (&lhs==&g && IsMonomial(h))
    {
      lhs.myAddMonomial(h);
      return;
    }
    if (&lhs==&h && IsMonomial(g))
    {
      lhs.myAddMonomial(g);
      return;
    }

    DistrMPolyInlPP::AutoPtrSummand SpareSummand(ans);
    SpareSummand.myRenew();
    while (gterm != nullptr && hterm != nullptr)
    {
      int cmp = lhs.myOrdvArith->myCmp(gterm->myOrdv, hterm->myOrdv);

      if (cmp < 0)
      {
	summand* hcopy = ans.myCopySummand(hterm);
	ans.myPushBack(hcopy);
	hterm = hterm->myNext;
	continue;
      }

      if (cmp > 0)
      {
	summand* gcopy = ans.myCopySummand(gterm);
	ans.myPushBack(gcopy);
	gterm = gterm->myNext;
	continue;
      }

      // Must have cmp == 0 here.
      // The leading PPs are the same, so we sum the coeffs.
      R->myAdd(SpareSummand->myCoeff, gterm->myCoeff, hterm->myCoeff);
      if (!R->myIsZero(SpareSummand->myCoeff))
      {
	lhs.myOrdvArith->myAssign(SpareSummand->myOrdv, gterm->myOrdv); // set PP ordv
	ans.myPushBack(SpareSummand.release());
	SpareSummand.myRenew();
      }
      gterm = gterm->myNext;
      hterm = hterm->myNext;
    }
    while (gterm != nullptr)
    {
      summand* gcopy = ans.myCopySummand(gterm);
      ans.myPushBack(gcopy);
      gterm = gterm->myNext;
    }
    while (hterm != nullptr)
    {
      summand* hcopy = ans.myCopySummand(hterm);
      ans.myPushBack(hcopy);
      hterm = hterm->myNext;
    }

    swap(lhs, ans); // really an assignment
  }


  // EXCEPTION SAFE
  void DistrMPolyInlPP::myAddMonomial(const DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(IsCompatible(*this, g));
    CoCoA_ASSERT(NumTerms(g)==1);

    typedef DistrMPolyInlPP::summand summand;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;
    AutoPtrSummand APs(*this);
    APs.myRenew();
    myCoeffRing->myAssign(APs->myCoeff, (g.mySummands)->myCoeff);
    myOrdvArith->myAssign(APs->myOrdv, (g.mySummands)->myOrdv);
    //    auto_ptr<summand> s(myCopySummand(g.mySummands));
    int CMP;

    while (f_smnd != nullptr &&
           (CMP = myOrdvArith->myCmp(f_smnd->myOrdv, APs->myOrdv)) >0)
      f_smnd = *(f_prev = &f_smnd->myNext);
    if (f_smnd == nullptr)  myPushBack(APs.release());
    else 
      if (CMP == 0)
      {
        myCoeffRing->myAdd(APs->myCoeff, APs->myCoeff, f_smnd->myCoeff);
        if (myCoeffRing->myIsZero(APs->myCoeff))
          myRemoveSummand(f_prev);
        else 
          myCoeffRing->mySwap(APs->myCoeff, f_smnd->myCoeff);
      }
      else // (CMP < 0)
      {
        myInsertSummand(APs.release(), f_prev);
        f_prev = &(*f_prev)->myNext;
      }
  }


  void sub(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h)
  {
    // This code is almost a copy of add(...).
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(g, h));
    const ring& R = lhs.myCoeffRing;
    DistrMPolyInlPP ans(lhs.myCoeffRing, lhs.myPPM, lhs.myOrdvArith, lhs.mySummandMemory);

    typedef DistrMPolyInlPP::summand summand;
    const summand* gterm = g.mySummands;
    const summand* hterm = h.mySummands;
    DistrMPolyInlPP::AutoPtrSummand SpareSummand(ans);
    SpareSummand.myRenew();
    while (gterm != nullptr && hterm != nullptr)
    {
      int ord = lhs.myOrdvArith->myCmp(gterm->myOrdv, hterm->myOrdv);

      if (ord < 0)
      {
	summand* hcopy = ans.myCopySummand(hterm);
	R->myNegate(hcopy->myCoeff, hcopy->myCoeff);  //??? can this throw?
	ans.myPushBack(hcopy);
	hterm = hterm->myNext;
	continue;
      }

      if (ord > 0)
      {
	summand* gcopy = ans.myCopySummand(gterm);
	ans.myPushBack(gcopy);
	gterm = gterm->myNext;
	continue;
      }

      // The leading PPs are the same, so we subtract the coeffs.
      R->mySub(SpareSummand->myCoeff, gterm->myCoeff, hterm->myCoeff);
      if (!R->myIsZero(SpareSummand->myCoeff))
      {
	lhs.myOrdvArith->myAssign(SpareSummand->myOrdv, gterm->myOrdv); // set PP ordv
	ans.myPushBack(SpareSummand.release());
	SpareSummand.myRenew();
      }
      gterm = gterm->myNext;
      hterm = hterm->myNext;
    }
    while (gterm != nullptr)
    {
      summand* gcopy = ans.myCopySummand(gterm);
      ans.myPushBack(gcopy);
      gterm = gterm->myNext;
    }
    while (hterm != nullptr)
    {
      summand* hcopy = ans.myCopySummand(hterm);
      R->myNegate(hcopy->myCoeff, hcopy->myCoeff);  //??? can this throw?
      ans.myPushBack(hcopy);
      hterm = hterm->myNext;
    }

    swap(lhs, ans); // really an assignment
  }


  bool div(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h)  // result is true iff quotient is exact.
  {
    PPMonoid PPM = lhs.myPPM;
//???    const PPOrdering ord = ordering(PPM);
    const ring& R = lhs.myCoeffRing;
    const DistrMPolyInlPP::summand* LMh = h.mySummands;
    const PPMonoidElem LPPden(LPP(h));
    DistrMPolyInlPP ans(lhs.myCoeffRing, lhs.myPPM, lhs.myOrdvArith, lhs.mySummandMemory);
    DistrMPolyInlPP dividend(g);
    while (!IsZero(dividend))
    {
      const DistrMPolyInlPP::summand* LMdividend = dividend.mySummands;
      DistrMPolyInlPP::AutoPtrSummand qterm(lhs);
      qterm.myRenew();
      if (!R->myIsDivisible(qterm->myCoeff, LMdividend->myCoeff, LMh->myCoeff)) return false;
      {
        //???        PPMonoidElem LPPnum(LPP(dividend));
        if (!IsDivisible(LPP(dividend), LPPden)) return false;
      }
      g.myOrdvArith->myDiv(qterm->myOrdv, LMdividend->myOrdv, LMh->myOrdv);  //??? check whether this overflows?
      R->myNegate(qterm->myCoeff, qterm->myCoeff);
      dividend.myAddMulSummand(qterm.get(), h, false);
      R->myNegate(qterm->myCoeff, qterm->myCoeff);
      ans.myPushBack(qterm.release());
    }
    swap(lhs, ans); // really an assignment
    return true;
  }


  void output(std::ostream& out, const DistrMPolyInlPP& f)  // for debugging only
  {
    if (!out) return;  // short-cut for bad ostreams

    if (IsZero(f)) { out << "0"; return; }
    const ring R = f.myCoeffRing;
    const PPMonoid PPM = f.myPPM;
    for (DistrMPolyInlPP::summand* it = f.mySummands; it != nullptr; it = it->myNext)
    {
      out << " +(";
      R->myOutput(out, it->myCoeff);
      out << ")*" << ConstRefPPMonoidElem(PPM, PPMonoidElemConstRawPtr(it->myOrdv));
    }
  }


//   bool IsConstant(const DistrMPolyInlPP& f)
//   {
//     if (IsZero(f)) return true;
//     if (!IsMonomial(f)) return false;
//     return IsOne(LPP(f));
//   }


  bool IsZero(const DistrMPolyInlPP& f)
  {
    return (f.mySummands == nullptr);
  }


//   bool IsOne(const DistrMPolyInlPP& f)
//   {
//     if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
//     if (!f.myOrdvArith->myIsZero(f.mySummands->myOrdv)) return false;
//     return f.myCoeffRing->myIsOne(f.mySummands->myCoeff);
//   }


//   bool IsMinusOne(const DistrMPolyInlPP& f)
//   {
//     if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
//     if (!f.myOrdvArith->myIsZero(f.mySummands->myOrdv)) return false;
//     return f.myCoeffRing->myIsMinusOne(f.mySummands->myCoeff);
//   }


//   bool IsConstant(const DistrMPolyInlPP& f)
//   {
//     if (IsZero(f)) return true;
//     if (f.mySummands->myNext != nullptr) return false; // NULL ptr
//     return f.myOrdvArith->myIsZero(f.mySummands->myOrdv);
//   }


//   bool IsIndet(std::size_t& index, const DistrMPolyInlPP& f)
//   {
//     if (IsZero(f)) return false;
//     if (f.mySummands->myNext != nullptr) return false; // NULL ptr
//     if (!f.myCoeffRing->myIsOne(f.mySummands->myCoeff)) return false;
//     return f.myOrdvArith->myIsIndet(index, f.mySummands->myOrdv);
//   }


  bool IsMonomial(const DistrMPolyInlPP& f)
  {
    if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
    return true;
  }


  bool IsEqual(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    if (&f == &g) return true;
    const DistrMPolyInlPP::summand* fterm = f.mySummands;
    const DistrMPolyInlPP::summand* gterm = g.mySummands;
    while (fterm != nullptr && gterm != nullptr)
    {
      if (!f.myIsEqual(fterm, gterm)) return false;
      fterm = fterm->myNext;
      gterm = gterm->myNext;
    }
    return fterm == gterm; // either both are 0 (when the polys are equal), or only one is 0
  }


//   bool IsEqual(const DistrMPolyInlPP& f, long n)
//   {
//     if (n == 0) return IsZero(f);
//     if (IsZero(f)) return IsZero(RingElem(f.myCoeffRing, n));
//     // From here on the polynomial is known to be non-zero
//     if (f.mySummands->myNext != nullptr) return false;
//     if (!f.myOrdvArith->myIsZero(f.mySummands->myOrdv)) return false;
//     return f.myCoeffRing->myIsEqual(f.mySummands->myCoeff, n);
//   }


//   void WeylMul(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h)
//   {

//   }


  void WeylDiv(DistrMPolyInlPP& /*lhs*/, const DistrMPolyInlPP& /*g*/, const DistrMPolyInlPP& /*h*/)
  {
  }


//   void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, std::size_t IndetIndex)
//   {
//     return deriv(lhs, f, indet(f.myPPM, IndetIndex));
//   }


//   void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, ConstRefPPMonoidElem x)
//   {
//     if (IsOne(x)) { lhs = f; return; }
//     const size_t nvars = NumIndets(owner(x));
// //???    const PPOrdering ord = ordering(owner(x));
//     const ring& R = f.myCoeffRing;
//     vector<long> expv(nvars);
//     exponents(expv, x);
//     vector<OrdvArith::OrdvElem> ordvx(OrdvWords(f.myOrdvArith));
//     f.myOrdvArith->myAssignFromExpv(&ordvx[0], expv);
// //clog<<"differentiating wrt expv: [";for(size_t i=0;i<nvars;++i)clog<<expv[i]<<" ";clog<<"]"<<std::endl;
//     DistrMPolyInlPP ans(f.myCoeffRing, f.myPPM, f.myOrdvArith, f.mySummandMemory);

//     for (const DistrMPolyInlPP::summand* f_term = f.mySummands; f_term != nullptr; f_term = f_term->myNext)
//     {
// //clog<<"LOOPHEAD\n";
//       BigInt scale(1);
//       for (size_t indet=0; indet < nvars; ++indet)
//       {
//         if (expv[indet] == 0) continue;
//         long d = f.myOrdvArith->myExponent(f_term->myOrdv, indet);
// //clog<<"log is "<<d<<" wrt var "<<indet<<std::endl;
//         if (d < expv[indet]) { scale = 0; break; }
//         scale *= RangeFactorial(d-expv[indet]+1, d);
//       }
// //if(IsZero(scale))clog<<"skipping term\n";
//       if (IsZero(scale)) continue;
// //clog<<"rescaling term by "<<scale<<std::endl;
//       DistrMPolyInlPP::AutoPtrSummand tmp(f);
//       tmp.myRenew();
//       R->myAssign(tmp->myCoeff, scale);
//       R->myMul(tmp->myCoeff, tmp->myCoeff, f_term->myCoeff);
//       if (R->myIsZero(tmp->myCoeff)) continue;
// //clog<<"dividing ordv [";for(size_t i=0;i<2;++i)clog<<f_term->myOrdv[i]<<" ";clog<<"]\n";
// //clog<<"by       ordv [";for(size_t i=0;i<2;++i)clog<<ordvx[i]<<" ";clog<<"]\n";
//       f.myOrdvArith->myDiv(tmp->myOrdv, f_term->myOrdv, &ordvx[0]);
// //clog<<"Quotient is   [";for(size_t i=0;i<2;++i)clog<<tmp->myOrdv[i]<<" ";clog<<"]\n";
//       ans.myPushBack(tmp.release());
//     }
//     swap(lhs, ans); // really an assignment
//   }


} // end of namespace CoCoA
