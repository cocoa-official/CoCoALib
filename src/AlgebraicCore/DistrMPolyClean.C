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


  // EXCEPTION SAFE
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
