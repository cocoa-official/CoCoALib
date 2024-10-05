//   Copyright (c)  2005-2009,2015,2017  John Abbott, Anna M. Bigatti

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


// Source code for class DistMPoly

#include "CoCoA/DistrMPolyClean.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"

#include <algorithm>
//using std::swap; in ourSwap
#include <iostream>
//using << in output
//#include <vector> ---  included in DistrMPolyClean.H
using std::vector;


namespace CoCoA
{

  bool IsCompatible(const DistrMPolyClean& f, const DistrMPolyClean& g)
  {
    return (f.myCoeffRing == g.myCoeffRing) &&
           (f.myPPM == g.myPPM) &&
           (&f.myMemMgr == &g.myMemMgr);
  }


  void DistrMPolyClean::ourSwap(DistrMPolyClean& f, DistrMPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    std::swap(f.mySummands, g.mySummands);
    std::swap(f.myEnd, g.myEnd);
    if (f.mySummands == nullptr) f.myEnd = &f.mySummands;  // CAREFUL if f or g is zero!
    if (g.mySummands == nullptr) g.myEnd = &g.mySummands;  //
  }


  // MemPool DistrMPolyClean::summand::ourMemMgr(sizeof(DistrMPolyClean::summand),
  //                                        "DistrMPolyClean::summand::ourMemMgr");


  void DistrMPolyClean::ourDeleteSummands(DistrMPolyClean::summand* ptr/*,const ring& R, const PPMonoid& M*/, MemPool& MemMgr)
  {
    while (ptr != nullptr)
    {
      DistrMPolyClean::summand* next = ptr->myNext;
      ptr->~summand();
      MemMgr.free(ptr);
      ptr = next;
    }
  }


  DistrMPolyClean::DistrMPolyClean(const ring& R, const PPMonoid& PPM, MemPool& MemMgr):
      myCoeffRing(R),
      myPPM(PPM),
      myMemMgr(MemMgr)
  {
    mySummands = nullptr;
    myEnd = &mySummands;
  }


  DistrMPolyClean::~DistrMPolyClean()
  {
    ourDeleteSummands(mySummands/*, myCoeffRing, myPPM*/, myMemMgr);
  }


  DistrMPolyClean::DistrMPolyClean(const DistrMPolyClean& copy):
      myCoeffRing(copy.myCoeffRing),
      myPPM(copy.myPPM),
      myMemMgr(copy.myMemMgr)
  {
    mySummands = nullptr;
    myEnd = &mySummands;
    // !!! THIS LOOP IS NOT EXCEPTION CLEAN !!!
    for (summand* it = copy.mySummands; it != nullptr; it = it->myNext)
      myPushBack(myCopySummand(it));
  }


  DistrMPolyClean& DistrMPolyClean::operator=(const DistrMPolyClean& rhs)
  {
    if (this == &rhs) return *this;
    DistrMPolyClean copy(rhs);
    ourSwap(*this, copy);
    return *this;
  }


  DistrMPolyClean& DistrMPolyClean::operator=(const MachineInt& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this; // to avoid needless alloc/free
    NewSummandPtr t(*this); t.myRenew();
    myCoeffRing->myAssign(raw(t->myCoeff), rhs);
    if (!IsZero(t->myCoeff))
    {
      mySummands = t.relinquish();
      myEnd = &mySummands->myNext;
    }
    return *this;
  }


  DistrMPolyClean& DistrMPolyClean::operator=(const BigInt& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this; // to avoid needless alloc/free
    NewSummandPtr t(*this); t.myRenew();
    myCoeffRing->myAssign(raw(t->myCoeff), rhs);
    if (!IsZero(t->myCoeff))
    {
      mySummands = t.relinquish();
      myEnd = &mySummands->myNext;
    }
    return *this;
  }

  DistrMPolyClean& DistrMPolyClean::operator=(const BigRat& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this; // to avoid needless alloc/free
    NewSummandPtr t(*this); t.myRenew();
    myCoeffRing->myAssign(raw(t->myCoeff), rhs);
    if (!IsZero(t->myCoeff))
    {
      mySummands = t.relinquish();
      myEnd = &mySummands->myNext;
    }
    return *this;
  }


//----------------------------------------------------------------------
// operations on summands


  DistrMPolyClean::summand* DistrMPolyClean::myCopySummand(const summand* original) const
  {
    NewSummandPtr copy(*this);
    copy.myRenew();

//    copy->myCoeff = original->myCoeff;
//    copy->myPP = original->myPP;
    myCoeffRing->myAssign(raw(copy->myCoeff), raw(original->myCoeff));
    myPPM->myAssign(raw(copy->myPP), raw(original->myPP));
    return copy.relinquish();
  }


  void DistrMPolyClean::myAssignZero()
  {
    ourDeleteSummands(mySummands/*, myCoeffRing, myPPM*/, myMemMgr);
    mySummands = nullptr;
    myEnd = &mySummands;
  }


  bool DistrMPolyClean::myIsEqual(const summand* const lhs, const summand* const rhs) const
  {
    return (myPPM->myIsEqual(raw(lhs->myPP), raw(rhs->myPP)) &&
            myCoeffRing->myIsEqual(raw(lhs->myCoeff), raw(rhs->myCoeff)));
  }


  long NumTerms(const DistrMPolyClean& f)
  {
    long nsummands = 0;
    for (const DistrMPolyClean::summand* it = f.mySummands; it != nullptr; it = it->myNext)
      ++nsummands;
    return nsummands;
  }


  ConstRefRingElem LC(const DistrMPolyClean& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    return f.mySummands->myCoeff;
  }


  void MoveLMToFront(DistrMPolyClean& f, DistrMPolyClean& g)
  {
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyClean::summand* ltg = g.mySummands;
    g.mySummands = g.mySummands->myNext;
    if (g.mySummands == nullptr) g.myEnd = &(g.mySummands);
    ltg->myNext = nullptr;
    f.myPushFront(ltg);
  }


  void MoveLMToBack(DistrMPolyClean& f, DistrMPolyClean& g)
  {
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyClean::summand* ltg = g.mySummands;
    g.mySummands = g.mySummands->myNext;
    if (g.mySummands == nullptr) g.myEnd = &(g.mySummands);
    ltg->myNext = nullptr;
    f.myPushBack(ltg);
  }


  void DistrMPolyClean::myDeleteLM()
  {
    CoCoA_ASSERT(!IsZero(*this));

    DistrMPolyClean::summand* old_lm = mySummands;
    mySummands = old_lm->myNext;
    if (mySummands == nullptr) myEnd = &mySummands;
    old_lm->myNext = nullptr;
    ourDeleteSummands(old_lm/*, myCoeffRing, myPPM*/, myMemMgr);
  }


//   void wdeg(degree& d, const DistrMPolyClean& f)
//   {
//     CoCoA_ASSERT(!IsZero(f));
//     f.myPPM->myWDeg(d, raw(f.mySummands->myPP));
//   }


//   int CmpWDeg(const DistrMPolyClean& f, const DistrMPolyClean& g)
//   {
//     CoCoA_ASSERT(!IsZero(f));
//     CoCoA_ASSERT(!IsZero(g));
//     if (!DistrMPolyClean::compatible(f, g))
//       CoCoA_THROW_ERROR("Incompatible polynomials", "CmpDeg(DistrMPolyClean,DistrMPolyClean)");
//     return f.myPPM->myCmpWDeg(raw(f.mySummands->myPP), raw(g.mySummands->myPP));
//   }


// This fn offers only the weak exception guarantee!!!
  void DistrMPolyClean::myAddMulSummand(const summand* s, const DistrMPolyClean& g, bool SkipLMg)  // this += s*g
  {
    CoCoA_ASSERT(IsCompatible(*this, g));
    CoCoA_ASSERT(!IsZero(s->myCoeff));

    const ring& R = myCoeffRing;

    const summand* g_smnd = g.mySummands;
    if (SkipLMg)    g_smnd = g_smnd->myNext;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;

    NewSummandPtr tmp_smnd(*this); tmp_smnd.myRenew();
    RingElem tmp(myCoeffRing);

    int CMP = 0;

    //    bool qIsOne = IsOne(s->myPP);
    bool qIsOne = false;

    for (; f_smnd != nullptr && g_smnd != nullptr; g_smnd = g_smnd->myNext)
    {
      R->myMul(raw(tmp_smnd->myCoeff), raw(s->myCoeff), raw(g_smnd->myCoeff));
      if (R->myIsZero(raw(tmp_smnd->myCoeff)))
        continue; // skip this g_smnd as coeff mults to 0
      myPPM->myMul(raw(tmp_smnd->myPP), raw(s->myPP), raw(g_smnd->myPP));
      if (qIsOne)
      {
        while (f_smnd != nullptr && (CMP=myPPM->myCmp(raw(f_smnd->myPP), raw(g_smnd->myPP))) >0)
          f_smnd = *(f_prev = &f_smnd->myNext);
        myPPM->myAssign(raw(tmp_smnd->myPP), raw(g_smnd->myPP));
      }
      else
      {
//        myPPM->myMul(raw(tmp_smnd->myPP), raw(s->myPP), raw(g_smnd->myPP));
        while (f_smnd != nullptr && (CMP=myPPM->myCmp(raw(f_smnd->myPP), raw(tmp_smnd->myPP))) >0)
          f_smnd = *(f_prev = &f_smnd->myNext);
      }
      if (f_smnd == nullptr)
      {
//        R->myMul(raw(tmp_smnd->myCoeff), raw(s->myCoeff), raw(g_smnd->myCoeff));
        myPushBack(tmp_smnd.relinquish());
        tmp_smnd.myRenew();
        g_smnd = g_smnd->myNext;
        break;
      }
      if (CMP == 0)
      {
//        if (R->myIsZeroAddMul(raw(f_smnd->myCoeff), raw(tmp), raw(s->myCoeff), raw(g_smnd->myCoeff)))
        R->myAdd(raw(f_smnd->myCoeff), raw(f_smnd->myCoeff), raw(tmp_smnd->myCoeff));
        if (R->myIsZero(raw(f_smnd->myCoeff)))
          myRemoveSummand(f_prev);  // f_prev = f_prev;
        else
          f_prev = &f_smnd->myNext;
        f_smnd = *f_prev;
      }
      else // (CMP < 0)
      {
//        R->myMul(raw(tmp_smnd->myCoeff), raw(s->myCoeff), raw(g_smnd->myCoeff));
        myInsertSummand(tmp_smnd.relinquish(), f_prev);
        tmp_smnd.myRenew();
        f_prev = &(*f_prev)->myNext;
      }
    }
    for (;g_smnd != nullptr; g_smnd = g_smnd->myNext)
    {
      R->myMul(raw(tmp_smnd->myCoeff), raw(s->myCoeff), raw(g_smnd->myCoeff));
      if (R->myIsZero(raw(tmp_smnd->myCoeff))) continue;
      myPPM->myMul(raw(tmp_smnd->myPP), raw(s->myPP), raw(g_smnd->myPP));
      myPushBack(tmp_smnd.relinquish());
      tmp_smnd.myRenew();
    }

//      clog << "AddMul: produced f=";output(clog,f);clog<<std::endl;
//      clog << "------------------------------------------------------"<<std::endl;

  }


  void DistrMPolyClean::myAddMulLM(const DistrMPolyClean& h, const DistrMPolyClean& g, bool SkipLMg)
  {                                                 //???
    if (IsZero(h)) CoCoA_THROW_ERROR(ERR::ReqNonZero, "DistrMPolyClean::myAddMul");
    myAddMulSummand(h.mySummands, g, SkipLMg);     //???
  }                                                 //???


//   void DistrMPolyClean::myWeylAddMulSummand(const summand* s, const DistrMPolyClean& g, bool SkipLMg)
//   {
//     CoCoA_ASSERT(IsCompatible(*this, g));

//     //    const PPOrdering& ord = ordering(myPPM);
//     const ring& R = myCoeffRing;
//     const size_t nvars = NumIndets(myPPM);
// //???    clog << "AddMul: Doing funny product of the following two polys" << std::endl;
// //???    output(clog, g);
// //???    CoCoA_ASSERT(myRefCount == 1);
//     //???    MakeWritable(f);

//     const summand* g_term = g.mySummands;
//     if (SkipLMg)    g_term = g_term->myNext;
// //???    summand** f_prev = &mySummands;
// //???    summand*  f_term = *f_prev;

//     DistrMPolyClean ppg = g;
//     vector<long> expv(nvars);
//     myPPM->myExponents(expv, raw(s->myPP));
// //???    clog << "expv: "; for (int i=0; i<myNumIndets;++i) clog << expv[i] << "  "; clog << std::endl;
//     for (size_t indet = nvars/2; indet < nvars; ++indet)
//     {
//       long n = expv[indet];
//       if (n == 0) continue;
// //???      clog << "AddMul: doing D variable with index " << indet - myNumIndets/2 << std::endl;
//       DistrMPolyClean der = ppg;

// //???      ppg *= IndetPower(myPPM, indet, n);
//       ppg.myMulByPP(raw(IndetPower(myPPM, indet, n)));
// //      mul(raw(ppg), raw(ppg), raw(IndetPower(indet, n)));

//       for (long i=1; i <= n; ++i)
//       {
//         deriv(der, der, indet-nvars/2);
// //???        deriv(raw(der), raw(der), indet-nvars/2);
// //???        clog << "der(" << i << ")="; output(clog, raw(der)); clog << std::endl;

// //        ppg += binomial(n, i)*der*IndetPower(myPPM, indet, n-i); // *IndetPower(myPPM, h, 2*i); // for homog case
//         auto_ptr<summand> jaa(new summand(myCoeffRing, myPPM));
//         R->myAssign(raw(jaa->myCoeff), binomial(n, i));
//         myPPM->myMulIndetPower(raw(jaa->myPP), indet, n-i);
// // if (HOMOG_CASE)        myPPM->mul(jaa->myPP, jaa->myPP, IndetPower(myPPM, h, 2*i));
//         ppg.myAddMulSummand(jaa.get(), der, false);
//       }
//     }
//     { // f *= indet^deg(pp, indet); for the early vars
//       for (size_t indet = nvars/2; indet < nvars; ++indet)
//         expv[indet] = 0;
//       auto_ptr<summand> jaa(new summand(myCoeffRing, myPPM));
//       R->myAssign(raw(jaa->myCoeff), raw(s->myCoeff));
//       myPPM->myAssign(raw(jaa->myPP), expv);
//       myAddMulSummand(jaa.get(), ppg, false);
//     }
//   }




  void DistrMPolyClean::myReductionStep(const DistrMPolyClean& g)
  {
    CoCoA_ASSERT(&g!=this);
    CoCoA_ASSERT(mySummands != nullptr);
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyClean tmp_poly(myCoeffRing, myPPM, myMemMgr);

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
    const RingElem gcdLC = gcd(LCf, LCg);
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


  void DistrMPolyClean::myReductionStepGCD(const DistrMPolyClean& g, RingElem& fscale)
  {
    CoCoA_ASSERT(&g!=this);
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyClean tmp_poly(myCoeffRing, myPPM, myMemMgr);
    RingElem gscale(myCoeffRing);

    ComputeFScaleAndGScale(LC(*this), LC(g), fscale, gscale);

    if ( !IsOne(fscale) )    myMulByCoeff(raw(fscale));
    DivLM(tmp_poly, *this, g);
    tmp_poly.myNegate();
    myDeleteLM();
    myAddMulLM(tmp_poly, g, /*SkipLMg = */ true);
  }


  void DistrMPolyClean::myAddClear(DistrMPolyClean& g) // sets g to 0 as a side-effect
  {
    const ring& R = myCoeffRing;
    //    const PPOrdering& ord = ordering(myPPM);
    typedef DistrMPolyClean::summand summand;

    summand*  g_smnd = g.mySummands;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;
    int CMP=0;
//???    CoCoA_ASSERT(*(G.myEnd)==nullptr);//BUG HUNTING  ???

   //    clog << "input f = "; output(clog, *this) ;clog << std::endl;
    while ( f_smnd!=nullptr && g_smnd!=nullptr )
    {
      while (f_smnd!=nullptr &&
             (CMP = myPPM->myCmp(raw(f_smnd->myPP), raw(g_smnd->myPP))) >0)
        f_smnd = *(f_prev = &f_smnd->myNext);
      if (f_smnd == nullptr)  break;
      //clog <<   "(myAddClear error: should never happen for Basic Reduction)" << std::endl;
      g.mySummands = g.mySummands->myNext;
      g_smnd->myNext = nullptr;
      if (CMP == 0)
      {
        R->myAdd(raw(f_smnd->myCoeff), raw(f_smnd->myCoeff), raw(g_smnd->myCoeff));
        if (IsZero(f_smnd->myCoeff))
          myRemoveSummand(f_prev);
        ourDeleteSummands(g_smnd/*, myCoeffRing, myPPM*/, myMemMgr);
      }
      else // (CMP < 0)
      {
        myInsertSummand(g_smnd, f_prev);
        f_prev = &(*f_prev)->myNext;
      }
      f_smnd = *f_prev;
      g_smnd = g.mySummands;
    }
    if (g.mySummands != nullptr)
    {
      *myEnd = g.mySummands;
      myEnd = g.myEnd;
      g.mySummands = nullptr;
    }
    g.myEnd = &g.mySummands;
    //    if (rare) {clog << "f2 = "; output(clog, f) ;clog << std::endl;}
  }


  void DistrMPolyClean::myAppendClear(DistrMPolyClean& g)  // sets g to 0; no throw guarantee!
  {
    if (g.mySummands == nullptr) return;
    *(myEnd) = g.mySummands;
    myEnd = g.myEnd;
    g.mySummands = nullptr;
    g.myEnd = &g.mySummands;
  }


  ConstRefPPMonoidElem LPP(const DistrMPolyClean& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    return f.mySummands->myPP;
  }


  void DivLM(DistrMPolyClean& lhs, const DistrMPolyClean& f, const DistrMPolyClean& g) // lhs = LM(f)/LM(g)
  {
    CoCoA_ASSERT(IsCompatible(f, g) && IsCompatible(lhs, f));
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    //    clog << "DivLM" << std::endl;
    const ring& R = f.myCoeffRing;    // shorthand
    typedef DistrMPolyClean::summand summand;            // shorthand
    const summand* const LMf = f.mySummands;  // shorthand
    const summand* const LMg = g.mySummands;  // shorthand

    DistrMPolyClean::NewSummandPtr SpareSummand(lhs); SpareSummand.myRenew();
    CoCoA_ASSERT( R->myIsDivisible(raw(SpareSummand->myCoeff),
                                   raw(LMf->myCoeff), raw(LMg->myCoeff)) );
    R->myDiv(raw(SpareSummand->myCoeff), raw(LMf->myCoeff), raw(LMg->myCoeff));
    (f.myPPM)->myDiv(raw(SpareSummand->myPP), raw(LMf->myPP), raw(LMg->myPP));
    lhs.myAssignZero();  // CANNOT myAssignZero() EARLIER in case lhs aliases f or g.
    lhs.myPushBack(SpareSummand.relinquish());
  }



  // This is lengthier than you might expect because it is exception safe.
  // It also works fine if c aliases a coefficient of the polynomial.
  void DistrMPolyClean::myMulByCoeff(RingElemConstRawPtr rawc)
  {
    CoCoA_ASSERT(!myCoeffRing->myIsZero(rawc));
    if (myCoeffRing->myIsOne(rawc)) return;
    const long n = NumTerms(*this);
    vector<RingElem> NewCoeffs(n, zero(myCoeffRing));
    long k=0;
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
    {
      myCoeffRing->myMul(raw(NewCoeffs[k]), raw(it->myCoeff), rawc);
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
        myCoeffRing->mySwap(raw(NewCoeffs[k]), raw(it->myCoeff));
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
  // It also works fine if rawc aliases a coefficient of the polynomial.
  bool DistrMPolyClean::myDivByCoeff(RingElemConstRawPtr rawc)
  {
    CoCoA_ASSERT(!myCoeffRing->myIsZero(rawc));
    if (myCoeffRing->myIsOne(rawc)) return true;
    vector<RingElem> NewCoeffs(NumTerms(*this), zero(myCoeffRing));
    long n=0;
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
    {
      if (!myCoeffRing->myIsDivisible(raw(NewCoeffs[n]), raw(it->myCoeff), rawc))
        return false;
      ++n;
    }
    // The new coeffs are in NewCoeffs now swap them into the poly.
    n = 0;
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
    {
      myCoeffRing->mySwap(raw(NewCoeffs[n]), raw(it->myCoeff));
      ++n;
    }
    return true;
  }


  // This would not be exception safe if PP multiplication might allocate memory.
  void DistrMPolyClean::myMulByPP(PPMonoidElemConstRawPtr rawpp)
  {
    for (summand* f_smnd = mySummands; f_smnd != nullptr ; f_smnd = f_smnd->myNext)
      myPPM->myMul(raw(f_smnd->myPP), raw(f_smnd->myPP), rawpp);
  }

// ??? WANT DIVISION BY A PP TOO???


//   void DistrMPolyClean::myWeylMul(PPMonoidElemConstRawPtr rawpp)
//   {
//     for (summand* f_smnd = mySummands; f_smnd != nullptr ; f_smnd = f_smnd->myNext)
//       myPPM->myMul(raw(f_smnd->myPP), raw(f_smnd->myPP), rawpp);
//   }


  void DistrMPolyClean::myPushFront(RingElemConstRawPtr rawc, const std::vector<long>& expv)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    NewSummandPtr tmp(*this); tmp.myRenew();
    myCoeffRing->myAssign(raw(tmp->myCoeff), rawc);
    myPPM->myAssign(raw(tmp->myPP), expv);
    myPushFront(tmp.relinquish()); 
  }


  void DistrMPolyClean::myPushBack(RingElemConstRawPtr rawc, const std::vector<long>& expv)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    NewSummandPtr tmp(*this); tmp.myRenew();
    myCoeffRing->myAssign(raw(tmp->myCoeff), rawc);
    myPPM->myAssign(raw(tmp->myPP), expv);
    myPushBack(tmp.relinquish());
  }


  void DistrMPolyClean::myPushFront(RingElemConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    NewSummandPtr tmp(*this); tmp.myRenew();
    myCoeffRing->myAssign(raw(tmp->myCoeff), rawc);
    myPPM->myAssign(raw(tmp->myPP), rawpp);
    myPushFront(tmp.relinquish()); 
  }


  void DistrMPolyClean::myPushBack(RingElemConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp)
  {
    if (myCoeffRing->myIsZero(rawc)) return;

    NewSummandPtr tmp(*this); tmp.myRenew();
    myCoeffRing->myAssign(raw(tmp->myCoeff), rawc);
    myPPM->myAssign(raw(tmp->myPP), rawpp);
    myPushBack(tmp.relinquish());
  }


  void DistrMPolyClean::myPushFront(summand* t)
  {
    CoCoA_ASSERT(t->myNext == nullptr);
    CoCoA_ASSERT(IsZero(*this) || myPPM->myCmp(raw(t->myPP), raw(mySummands->myPP)) > 0);
    t->myNext = mySummands;
    mySummands = t;
    if (myEnd == &mySummands) myEnd = &t->myNext;
  }


  void DistrMPolyClean::myPushBack(summand* t)
  {
    CoCoA_ASSERT(t->myNext == nullptr);
    // Cannot easily check that t really is smaller than smallest PP in the polynomial.
#ifdef CoCoA__DEBUG
    // Costly check that t really is smaller than smallest PP in the DistrMPolyInlFpPP.
    if (!IsZero(*this))
    {
      summand* LastSummand = mySummands;
      while (LastSummand->myNext != nullptr)
        LastSummand = LastSummand->myNext;
      CoCoA_ASSERT(cmp(t->myPP, LastSummand->myPP) < 0);
      CoCoA_ASSERT(myPPM->myCmp(t->myPP, LastSummand->myPP) < 0);
    }
#endif
    *myEnd = t;
    myEnd = &t->myNext;
  }


  void DistrMPolyClean::myRemoveSummand(summand** prev_link)
  {
    summand* DeleteMe = *prev_link;
    CoCoA_ASSERT(DeleteMe != nullptr);
    if (DeleteMe->myNext == nullptr)
      myEnd = prev_link;

    *prev_link = DeleteMe->myNext;
    DeleteMe->myNext = nullptr;
    ourDeleteSummands(DeleteMe/*, myCoeffRing, myPPM*/, myMemMgr);
  }


  void DistrMPolyClean::myInsertSummand(summand* s, summand** prev_link)
  {
    s->myNext = (*prev_link);
    (*prev_link) = s;
    if (myEnd == prev_link) myEnd = &(s->myNext);
  }


  bool IsZeroAddLCs(DistrMPolyClean& f, DistrMPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    CoCoA_ASSERT( LPP(f)==LPP(g) );
    f.myCoeffRing->myAdd(raw(f.mySummands->myCoeff), raw(f.mySummands->myCoeff), raw(g.mySummands->myCoeff));
    g.myDeleteLM();
    if (!IsZero(f.mySummands->myCoeff)) return false;
    f.myDeleteLM();
    return true;
  }


  void DistrMPolyClean::myNegate()
  {
    for (summand* iter = mySummands; iter != nullptr; iter = iter->myNext)
      myCoeffRing->myNegate(raw(iter->myCoeff), raw(iter->myCoeff));   // MIGHT THROW???
  }


  void add(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h)
  {
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(g, h));
    const ring& R = lhs.myCoeffRing;
    DistrMPolyClean ans(lhs.myCoeffRing, lhs.myPPM, lhs.myMemMgr);

    typedef DistrMPolyClean::summand summand;
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

    DistrMPolyClean::NewSummandPtr SpareSummand(lhs); SpareSummand.myRenew();
    while (gterm != nullptr && hterm != nullptr)
    {
      int cmp = (lhs.myPPM)->myCmp(raw(gterm->myPP), raw(hterm->myPP));

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
      R->myAdd(raw(SpareSummand->myCoeff), raw(gterm->myCoeff), raw(hterm->myCoeff));
      if (!IsZero(SpareSummand->myCoeff))
      {
	(lhs.myPPM)->myAssign(raw(SpareSummand->myPP), raw(gterm->myPP)); // set PP ordv
	ans.myPushBack(SpareSummand.relinquish());
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
  void DistrMPolyClean::myAddMonomial(const DistrMPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(*this, g));
    CoCoA_ASSERT(NumTerms(g)==1);

    typedef DistrMPolyClean::summand summand;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;
    NewSummandPtr s(*this); s.myRenew();
    s->myCoeff = g.mySummands->myCoeff; s->myPP = g.mySummands->myPP;

    int CMP;
    while (f_smnd!=nullptr &&
           (CMP = myPPM->myCmp(raw(f_smnd->myPP), raw(s->myPP))) >0)
      f_smnd = *(f_prev = &f_smnd->myNext);
    if (f_smnd == nullptr)  { myPushBack(s.relinquish()); return; }

    if (CMP < 0)
    {
      myInsertSummand(s.relinquish(), f_prev);
      f_prev = &(*f_prev)->myNext;
      return;
    }
      
    myCoeffRing->myAdd(raw(s->myCoeff), raw(f_smnd->myCoeff), raw(s->myCoeff));
    if (IsZero(s->myCoeff))
      myRemoveSummand(f_prev);
    else 
      myCoeffRing->mySwap(raw(s->myCoeff), raw(f_smnd->myCoeff));
  }
  

  void sub(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h)
  {
    // This code is almost a copy of add(...).
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(g, h));
    const ring& R = lhs.myCoeffRing;
    DistrMPolyClean ans(lhs.myCoeffRing, lhs.myPPM, lhs.myMemMgr);

    typedef DistrMPolyClean::summand summand;
    const summand* gterm = g.mySummands;
    const summand* hterm = h.mySummands;
    DistrMPolyClean::NewSummandPtr SpareSummand(lhs); SpareSummand.myRenew();
    while (gterm != nullptr && hterm != nullptr)
    {
      int ord = (lhs.myPPM)->myCmp(raw(gterm->myPP), raw(hterm->myPP));

      if (ord < 0)
      {
	summand* hcopy = ans.myCopySummand(hterm);
	R->myNegate(raw(hcopy->myCoeff), raw(hcopy->myCoeff));  //??? can this throw?
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
      R->mySub(raw(SpareSummand->myCoeff), raw(gterm->myCoeff), raw(hterm->myCoeff));
      if (!IsZero(SpareSummand->myCoeff))
      {
	(lhs.myPPM)->myAssign(raw(SpareSummand->myPP), raw(gterm->myPP)); // set PP ordv
	ans.myPushBack(SpareSummand.relinquish());
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
      R->myNegate(raw(hcopy->myCoeff), raw(hcopy->myCoeff));  //??? can this throw?
      ans.myPushBack(hcopy);
      hterm = hterm->myNext;
    }

    swap(lhs, ans); // really an assignment
  }


//   bool div(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h)  // result is true iff quotient is exact.
//   {
//     if (IsZero(h)) return false; //???CoCoA_THROW_ERROR(ERR::DivByZero,"div(DMP,DMP,DMP)");
//     PPMonoid PPM = lhs.myPPM;
//     //    const PPOrdering ord = ordering(PPM);
//     const ring& R = lhs.myCoeffRing;
//     const DistrMPolyClean::summand* LMh = h.mySummands;
//     const PPMonoidElem LPPden(LMh->myPP);
//     DistrMPolyClean ans(lhs.myCoeffRing, lhs.myPPM, lhs.myMemMgr);
//     DistrMPolyClean dividend(g);
//     while (!IsZero(dividend))
//     {
//       const DistrMPolyClean::summand* LMdividend = dividend.mySummands;
//       auto_ptr<DistrMPolyClean::summand> qterm(new DistrMPolyClean::summand(R, lhs.myPPM));
//       if (!R->myIsDivisible(raw(qterm->myCoeff), raw(LMdividend->myCoeff), raw(LMh->myCoeff))) return false;
//       {
//         PPMonoidElem LPPnum(LMdividend->myPP);
//         if (!IsDivisible(LPPnum, LPPden)) return false;
//       }
//       PPM->myDiv(raw(qterm->myPP), raw(LMdividend->myPP), raw(LMh->myPP));  //??? check whether this fails!
//       R->myNegate(raw(qterm->myCoeff), raw(qterm->myCoeff));
//       dividend.myAddMulSummand(qterm.get(), h, false);
//       R->myNegate(raw(qterm->myCoeff), raw(qterm->myCoeff));
//       ans.myPushBack(qterm.release());
//     }
//     swap(lhs, ans); // really an assignment
//     return true;
//   }


  void output(std::ostream& out, const DistrMPolyClean& f)  // for debugging only
  {
    if (!out) return;  // short-cut for bad ostreams
    if (IsZero(f)) { out << "0"; return; }
    for (DistrMPolyClean::summand* it = f.mySummands; it != nullptr; it = it->myNext)
      out << " +(" << it->myCoeff << ")*" << it->myPP;
  }


//   bool IsConstant(const DistrMPolyClean& f)
//   {
//     if (IsZero(f)) return true;
//     if (!IsMonomial(f)) return false;
//     return IsOne(LPP(f));
//   }


  bool IsZero(const DistrMPolyClean& f)
  {
    return (f.mySummands == nullptr);
  }


//   bool IsOne(const DistrMPolyClean& f)
//   {
//     if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
//     if (!IsOne(f.mySummands->myPP)) return false;
//     return IsOne(f.mySummands->myCoeff);
//   }


//   bool IsMinusOne(const DistrMPolyClean& f)
//   {
//     if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
//     if (!IsOne(f.mySummands->myPP)) return false;
//     return IsMinusOne(f.mySummands->myCoeff);
//   }


//   bool IsConstant(const DistrMPolyClean& f)
//   {
//     if (IsZero(f)) return true;
//     if (f.mySummands->myNext != nullptr) return false; // NULL ptr
//     return IsOne(f.mySummands->myPP);
//   }


//   bool IsIndet(std::size_t& IndetIndex, const DistrMPolyClean& f)
//   {
//     if (IsZero(f)) return false;
//     if (f.mySummands->myNext != nullptr) return false; // NULL ptr
//     if (!IsOne(f.mySummands->myCoeff)) return false;
//     return IsIndet(IndetIndex, f.mySummands->myPP);
//   }


  bool IsMonomial(const DistrMPolyClean& f)
  {
    if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
    return true;
  }


  bool IsEqual(const DistrMPolyClean& f, const DistrMPolyClean& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    if (&f == &g) return true;
    const DistrMPolyClean::summand* fterm = f.mySummands;
    const DistrMPolyClean::summand* gterm = g.mySummands;
    while (fterm != nullptr && gterm != nullptr)
    {
      if (!f.myIsEqual(fterm, gterm)) return false;
      fterm = fterm->myNext;
      gterm = gterm->myNext;
    }
    return fterm == gterm; // either both are nullptr (when the polys are equal), or only one is nullptr
  }


//   void WeylMul(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h)
//   {

//   }


//   void WeylDiv(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h)
//   {
//     CoCoA_THROW_ERROR(ERR::NYI, "WeylDiv (DistrMPolyClean)");
//   }


//   void deriv(DistrMPolyClean& lhs, const DistrMPolyClean& f, std::size_t IndetIndex)
//   {
//     deriv(lhs, f, indet(f.myPPM, IndetIndex));
//   }


//   void deriv(DistrMPolyClean& lhs, const DistrMPolyClean& f, ConstRefPPMonoidElem x)
//   {
//     if (IsOne(x)) { lhs = f; return; }
//     const size_t nvars = NumIndets(owner(x));
//     //    const PPOrdering ord = ordering(owner(x));
//     const PPMonoid PPM = owner(x);
//     const ring& R = f.myCoeffRing;
//     vector<long> expv(nvars);
//     exponents(expv, x);
//     //    vector<PPOrderingBase::OrdvElem> ordvx(OrdvWords(ord));
//     //    PPM->myInit(&ordvx[0], &expv[0]);
// //clog<<"differentiating wrt expv: [";for(size_t i=0;i<nvars;++i)clog<<expv[i]<<" ";clog<<"]"<<std::endl;
//     DistrMPolyClean ans(f.myCoeffRing, f.myPPM, f.myMemMgr);

//     for (const DistrMPolyClean::summand* f_term = f.mySummands; f_term != nullptr; f_term = f_term->myNext)
//     {
// //clog<<"LOOPHEAD\n";
//       BigInt scale(1);
//       for (size_t indet=0; indet < nvars; ++indet)
//       {
//         if (expv[indet] == 0) continue;
//         long d = PPM->myExponent(raw(f_term->myPP), indet);
// //clog<<"log is "<<d<<" wrt var "<<indet<<std::endl;
//         if (d < expv[indet]) { scale = 0; break; }
//         scale *= RangeFactorial(d-expv[indet]+1, d);
//       }
// //if(IsZero(scale))clog<<"skipping term\n";
//       if (IsZero(scale)) continue;
// //clog<<"rescaling term by "<<scale<<std::endl;
//       auto_ptr<DistrMPolyClean::summand> tmp(new DistrMPolyClean::summand(R, PPM));
//       R->myAssign(raw(tmp->myCoeff), scale);
//       R->myMul(raw(tmp->myCoeff), raw(tmp->myCoeff), raw(f_term->myCoeff));
//       if (IsZero(tmp->myCoeff)) continue;
// //clog<<"dividing ordv [";for(size_t i=0;i<2;++i)clog<<f_term->myPP[i]<<" ";clog<<"]\n";
// //clog<<"by       ordv [";for(size_t i=0;i<2;++i)clog<<ordvx[i]<<" ";clog<<"]\n";
//       PPM->myDiv(raw(tmp->myPP), raw(f_term->myPP), raw(x));
// //clog<<"Quotient is   [";for(size_t i=0;i<2;++i)clog<<tmp->myPP[i]<<" ";clog<<"]\n";
//       ans.myPushBack(tmp.release());
//     }
//     swap(lhs, ans); // really an assignment
//   }



//  //----------------------------------------------------------------------

//    PolyIter::PolyIter(const PolyRing& Rx, DistrMPolyClean::summand** ptrptr):
//      myPolyRing(Rx),
//      myPtrPtr(ptrptr),
//      myExpv(NumIndets(Rx))
//    {}


//    PolyIter::PolyIter(const PolyIter& copy):
//      myPolyRing(copy.myPolyRing),
//      myPtrPtr(copy.myPtrPtr),
//      myExpv(NumIndets(myPolyRing))
//    {}


//    PolyIter& PolyIter::operator=(const PolyIter& rhs)
//    {
//      CoCoA_ASSERT(&myPolyRing == &rhs.myPolyRing);
//      myPtrPtr = rhs.myPtrPtr;
//      return *this;
//    }


//    PolyIter& PolyIter::operator++()
//    {
//      CoCoA_ASSERT(*myPtrPtr != nullptr);
//      myPtrPtr = &((*myPtrPtr)->myNext);
//      return *this;
//    }


//    PolyIter PolyIter::operator++(int)
//    {
//      CoCoA_ASSERT(*myPtrPtr != nullptr);
//      PolyIter copy(*this);
//      myPtrPtr = &((*myPtrPtr)->myNext);
//      return copy;
//    }


//    bool operator==(const PolyIter& lhs, const PolyIter& rhs)
//    {
//      return lhs.myPtrPtr == rhs.myPtrPtr;
//    }


//    bool operator!=(const PolyIter& lhs, const PolyIter& rhs)
//    {
//      return lhs.myPtrPtr != rhs.myPtrPtr;
//    }


//    RingElemRawPtr RawCoeff(PolyIter& i)
//    {
//      CoCoA_ASSERT(*i.myPtrPtr != nullptr);
//      return (*i.myPtrPtr)->myCoeff;
//    }


//    const RingElem coeff(const PolyIter& i)
//    {
//      CoCoA_ASSERT(*i.myPtrPtr != nullptr);
//      return RingElem(CoeffRing(i.myPolyRing), (*i.myPtrPtr)->myCoeff);
//    }


//    PPMonoidElem pp(PolyIter& i)
//    {
//      CoCoA_ASSERT(*i.myPtrPtr != nullptr);
//      return PPMonoidElem(PPM(i.myPolyRing), PPMonoidElem::FromOrdv, (*i.myPtrPtr)->myPP);
//    }


//    PPOrdering::OrdvElem* ordv(PolyIter& i)
//    {
//      CoCoA_ASSERT(*i.myPtrPtr != nullptr);
//      return (*i.myPtrPtr)->myPP;
//    }



} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DistrMPolyClean.C,v 1.32 2024/07/02 09:57:57 bigatti Exp $
// $Log: DistrMPolyClean.C,v $
// Revision 1.32  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.31  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.30  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.29  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.28  2020/02/11 16:12:17  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.27  2019/10/15 11:54:08  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.26  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.25  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.24  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.23  2017/12/18 12:07:26  abbott
// Summary: Correct handling for polys with coeffs in rings with zero-divisors
//
// Revision 1.22  2015/11/24 17:46:56  abbott
// Summary: Corrected MoveLMToFront
//
// Revision 1.21  2015/11/04 10:42:05  abbott
// Summary: Replaced auto_ptr<summand>  with  NewSummandPtr
//
// Revision 1.20  2015/04/27 10:30:51  bigatti
// -- commented out excessive assert in myPushFront
//
// Revision 1.19  2015/04/27 10:06:04  bigatti
// -- simplified myAppendClear
// -- fixed myMoveLMToBack
// -- added CoCoA_ASSERT to myMoveLMToBack and myMoveLMToFront
//
// Revision 1.18  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.17  2015/04/17 16:24:07  abbott
// Summary: Added check that first arg is non-zero
// Author: JAA
//
// Revision 1.16  2014/04/30 16:05:21  abbott
// Summary: Replaced size_t by long
// Author: JAA
//
// Revision 1.15  2014/01/28 10:57:24  bigatti
// -- removed useless ==1 on boolean
//
// Revision 1.14  2012/10/24 12:13:34  abbott
// Added keyword const to local variable in ComputeFScaleAndGScale.
//
// Revision 1.13  2012/10/16 10:27:34  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.12  2012/10/11 14:27:59  abbott
// Removed "semantically risky" function RefLC/RawLC from DistrMPoly*.
// Reimplemented myRecvTwinFloat in RingDistrMPoly* more cleanly (but
// new impl does make a wasteful copy of the coeff).
//
// Revision 1.11  2012/10/05 15:33:28  bigatti
// -- added myAddMonomial
//
// Revision 1.10  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.9  2012/01/25 13:43:45  bigatti
// -- moved up: IsCompatible, ourSwap
//
// Revision 1.8  2011/11/09 14:03:59  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.7  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.6  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.5  2011/07/15 16:54:04  abbott
// Minor correction to a comment.
//
// Revision 1.4  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.3  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.2  2010/12/20 15:19:29  bigatti
// -- modified IsZeroAddMul with temporary variable (slug found with cyclotomic)
//
// Revision 1.1  2010/10/08 08:05:28  bigatti
// -- renamed (Ring)DistrMPoly --> (Ring)DistrMPolyClean
//
// Revision 1.12  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.11  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.10  2009/09/28 17:14:41  bigatti
// -- commented out unused functions (div, deriv, *Weyl*)
//
// Revision 1.9  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.8  2008/04/10 15:15:32  bigatti
// -- added  void myPushFront(rawc, rawpp)
//
// Revision 1.7  2007/12/05 12:11:07  bigatti
// -- cleaning (mostly removing unused code)
//
// Revision 1.6  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.5  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
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
// Revision 1.17  2007/03/08 18:22:30  cocoa
// Just whitespace cleaning.
//
// Revision 1.16  2007/03/07 13:42:45  bigatti
// -- removed useless argument and other minor changes
//
// Revision 1.15  2007/01/15 13:34:30  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.14  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPolyClean(..)
//
// Revision 1.13  2006/11/23 18:01:53  cocoa
// -- moved printing functions in unified implementation in SparsePolyRing
// -- simplified "output(f)" for debugging only
//
// Revision 1.12  2006/11/22 15:11:36  cocoa
// -- added #include "CoCoA/symbol.H"
//
// Revision 1.11  2006/11/21 18:09:24  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPolyClean(..) and RingDistrMPolyClean(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.10  2006/11/14 17:17:13  cocoa
// -- fixed coding convention "our"
//
// Revision 1.9  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPolyClean
//
// Revision 1.8  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.7  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyCleanInlPP.
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
// -- RingDistrMPolyClean*.H  have been "moved" into RingDistrMPolyClean*.C
// -- some coding conventions fixed in DistrMPolyClean*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.10  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPolyClean* and DistrMPolyClean* have been disabled
//
// Revision 1.9  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.8  2006/03/17 18:08:51  cocoa
// -- changed: mul --> myMulByPP
//
// Revision 1.7  2006/03/16 17:52:16  cocoa
// -- changed: mul, div --> myMulByCoeff, myDivByCoeff
//
// Revision 1.6  2006/03/16 13:13:20  cocoa
// -- added: CoCoA_ASSERT in DivLM
//
// Revision 1.5  2006/03/12 21:28:34  cocoa
// Major check in after many changes
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
// Revision 1.6  2005/10/12 15:52:09  cocoa
// Completed test-RingFp1 and corrected/cleaned the SmallFp*
// and RingFp* files.
//
// Some minor tidying elsewhere.
//
// Revision 1.5  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.4  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.3  2005/07/01 16:08:16  cocoa
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
// Revision 1.13  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.12  2004/11/12 16:25:57  cocoa
// -- printing aligned with DistrMPolyCleanInlFpPP and DistrMPolyCleanInlPP
//
// Revision 1.11  2004/11/11 13:30:52  cocoa
// -- minor changes for doxygen
// -- changed: cout --> GlobalLogput()
//
// Revision 1.10  2004/11/11 11:56:09  cocoa
// (1) Tidied makefiles, and introduced common.mki
// (2) Improved several tests, and cleaned them so that they
//     handle sanely any otherwise unhandled exceptions.
//
// Revision 1.9  2004/11/08 13:56:02  cocoa
// -- changed calls to ZZ (after changes to ZZ.H)
//
// Revision 1.8  2004/11/02 15:59:11  cocoa
// -- fixed: memory leak on summand (auto_ptr)
//
// Revision 1.7  2004/10/28 16:02:08  cocoa
// -- changed one last PPOrdering::ExpvElem into SmallExponent_t
//
// Revision 1.6  2004/10/21 17:16:37  cocoa
// Fairly major change: new OrdvArith namspace with various members,
//   new global typedef  SmallExponent_t (defined in config.H).
//
// Revision 1.5  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.4  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.3  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.2  2004/01/28 16:27:00  cocoa
// Added the necessary for CmpDeg to work.
//
// Revision 1.1  2003/11/21 14:33:55  cocoa
// -- First Import
//
