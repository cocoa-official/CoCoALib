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


  // EXEPTION SAFE
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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DistrMPolyInlPP.C,v 1.39 2024/08/02 16:09:14 bigatti Exp $
// $Log: DistrMPolyInlPP.C,v $
// Revision 1.39  2024/08/02 16:09:14  bigatti
// Summary: using CoCoA_THROW_ERROR1/2
//
// Revision 1.38  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.37  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.36  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.35  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.34  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.33  2020/02/11 16:12:17  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.32  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.31  2017/12/18 12:07:26  abbott
// Summary: Correct handling for polys with coeffs in rings with zero-divisors
//
// Revision 1.30  2016/05/18 12:23:56  abbott
// Summary: Major change: replaced myEnd with myLastSummand
//
// Revision 1.29  2015/04/27 10:30:51  bigatti
// -- commented out excessive assert in myPushFront
//
// Revision 1.28  2015/04/27 10:06:04  bigatti
// -- simplified myAppendClear
// -- fixed myMoveLMToBack
// -- added CoCoA_ASSERT to myMoveLMToBack and myMoveLMToFront
//
// Revision 1.27  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.26  2015/04/17 16:20:36  abbott
// Summary: Corrected myAddMul (but see redmine #686)
// Author: JAA
//
// Revision 1.25  2014/04/30 16:05:55  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.24  2014/01/28 10:57:51  bigatti
// -- removed useless ==1 on boolean
//
// Revision 1.23  2012/10/24 12:13:51  abbott
// Changed return type of LC.
//
// Revision 1.22  2012/10/16 10:28:40  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.21  2012/10/11 14:27:59  abbott
// Removed "semantically risky" function RefLC/RawLC from DistrMPoly*.
// Reimplemented myRecvTwinFloat in RingDistrMPoly* more cleanly (but
// new impl does make a wasteful copy of the coeff).
//
// Revision 1.20  2012/10/05 15:35:13  bigatti
// -- changed renew --> myRenew
// -- added myAddMonomial
//
// Revision 1.19  2012/01/25 13:44:40  bigatti
// -- moved up: IsCompatible, ourSwap
// -- minor cleaning in includes
//
// Revision 1.18  2011/11/09 14:03:59  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.17  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.16  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.15  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.14  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.13  2010/12/20 15:19:29  bigatti
// -- modified IsZeroAddMul with temporary variable (slug found with cyclotomic)
//
// Revision 1.12  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.11  2009/09/28 17:14:41  bigatti
// -- commented out unused functions (div, deriv, *Weyl*)
//
// Revision 1.10  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.9  2008/04/10 15:15:32  bigatti
// -- added  void myPushFront(rawc, rawpp)
//
// Revision 1.8  2007/12/05 12:11:07  bigatti
// -- cleaning (mostly removing unused code)
//
// Revision 1.7  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.6  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.5  2007/09/28 14:19:49  bigatti
// -- added "my" to some member functions
// -- added some comments
// -- commented out some unused functions (implemented in SparsePolyRing)
//
// Revision 1.4  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.3  2007/05/21 14:50:56  bigatti
// -- myPushFront and myPushBack now accept zero coefficient
//
// Revision 1.2  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.16  2007/03/08 18:22:30  cocoa
// Just whitespace cleaning.
//
// Revision 1.15  2007/03/07 13:42:45  bigatti
// -- removed useless argument and other minor changes
//
// Revision 1.14  2007/01/15 13:34:30  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.13  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.12  2006/11/24 17:01:43  cocoa
// -- reorganized includes of header files
//
// Revision 1.11  2006/11/23 18:01:53  cocoa
// -- moved printing functions in unified implementation in SparsePolyRing
// -- simplified "output(f)" for debugging only
//
// Revision 1.10  2006/11/22 15:11:36  cocoa
// -- added #include "CoCoA/symbol.H"
//
// Revision 1.9  2006/11/21 18:09:23  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.8  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.7  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.6  2006/10/06 10:01:21  cocoa
// Consequential changes from the modifications to the header files.
//
// Revision 1.5  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.4  2006/07/20 17:06:08  cocoa
// -- moved myStdDeg into SparsePolyRing
//
// Revision 1.3  2006/06/22 14:07:18  cocoa
// Minor cleaning and elimination of useless #includes.
//
// Revision 1.2  2006/06/08 16:45:28  cocoa
// -- RingDistrMPoly*.H  have been "moved" into RingDistrMPoly*.C
// -- some coding conventions fixed in DistrMPoly*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.14  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.13  2006/04/27 12:19:26  cocoa
// -- coding conventions (myAddMulSummand)
//
// Revision 1.12  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.11  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.10  2006/03/21 17:55:08  cocoa
// -- changed: "my.." coding conventions
//
// Revision 1.9  2006/03/20 17:28:37  cocoa
// -- commented out RawLC
//
// Revision 1.8  2006/03/17 18:08:51  cocoa
// -- changed: mul --> myMulByPP
//
// Revision 1.7  2006/03/16 17:56:14  cocoa
// -- changed: mul, div --> myMulByCoeff, myDivByCoeff (+check IsOne(c))
// -- changed: myOutputSelf for homomorphism
//
// Revision 1.6  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.5  2006/03/07 16:27:18  cocoa
// -- fixed: LPP now works correctly
//
// Revision 1.4  2006/03/07 10:06:12  cocoa
// -- fixed: PPMonoidElem LPP(f) now returns ConstRefPPMonoidElem
//
// Revision 1.3  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.2  2006/02/13 13:17:40  cocoa
// -- fixed: "const PPMonoidElem&" --> "ConstRefPPMonoidElem"
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.8  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.7  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
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
// Revision 1.4  2005/07/08 15:09:29  cocoa
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
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:04  cocoa
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
// Revision 1.18  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.17  2004/11/11 13:32:03  cocoa
// -- minor changes for doxygen
// -- change: cout --> GlobalLogput()
//
// Revision 1.16  2004/11/04 18:47:43  cocoa
// (1) Ring member functions which previously expected mpz_t args
//     now expect ZZ args.  Numerous minor consequential changes.
// (2) Renamed function which gives access to the mpz_t value inside
//     a ZZ object: previously was raw(...), now is mpzref(...).
//     Plenty of calls had to be altered.
//
// Revision 1.15  2004/11/02 18:21:21  cocoa
// -- changed: myGetExpvBuffer --> myExpvBufferRef
//
// Revision 1.14  2004/11/02 15:52:06  cocoa
// -- changed LPP body: now it uses myGetExpvBuffer
//
// Revision 1.13  2004/11/02 15:18:13  cocoa
// -- call to myPrintVarName --> myOutputIndetName
//
// Revision 1.12  2004/10/29 15:26:18  cocoa
// -- code fixed for compatibility with OrdvArith
//
// Revision 1.10  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.9  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.8  2004/03/20 17:46:11  cocoa
// Check in prior to departure to RWCA
//
// Revision 1.7  2004/01/28 16:27:00  cocoa
// Added the necessary for CmpDeg to work.
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
// Revision 1.3  2003/10/09 12:48:18  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:32  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.9  2003/09/17 15:17:47  abbott
// Added missing "reference" to third arg of deriv.
//
// Revision 1.8  2003/06/23 17:03:51  abbott
// Minor cleaning prior to public release.
// Name change, and consequences.
// Corrected division which blindy assumed inputs were divisible.
//
// Revision 1.7  2003/05/28 08:50:36  bigatti
// - fixed: output(...) for PolyRing coefficients
//
// Revision 1.6  2003/05/27 15:13:08  bigatti
// - fixed: output(...): -1*x[0]
//
// Revision 1.5  2003/05/27 12:35:58  bigatti
// - improved output(ostream& out, const DMPI& f)
//
// Revision 1.4  2003/05/27 12:18:42  abbott
// Added a proper implementation of div.
//
// Revision 1.3  2003/05/08 10:41:48  abbott
// Checking in prior to a probable major overhaul (this code is a mess).
// Numerous changes, mostly due to name changes for rings and PPMonoids.
// Added function deriv to compute derivatives.
// At least it works as it stands; well simple tests work OK.
//
// Revision 1.2  2002/11/18 17:54:52  bigatti
// - renamed: AddMulSummand
//
// Revision 1.1  2002/11/15 15:50:39  abbott
// Initial revision
//
//
