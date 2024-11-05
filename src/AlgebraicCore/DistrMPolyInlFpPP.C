//   Copyright (c)  2005-2012  Anna Bigatti

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


// Source code for class DistMPolyInlFpPP

#include "CoCoA/DistrMPolyInlFpPP.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/error.H"


//using std::swap;
#include <iostream>
//using "<<"
//#include <vector>
using std::vector;


namespace CoCoA
{

  bool IsCompatible(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g)
  {
    return //(f.myILCoeffImpl == g.myILCoeffImpl) && // ANNA: this should work!!!
      (f.myPPM == g.myPPM) &&
      (&f.mySummandMemory == &g.mySummandMemory);
  }


  // I have to make my own swap function as I cannot declare the template
  // specialization (of std::swap) to be a friend; std::swap will call this fn.
  void DistrMPolyInlFpPP::ourSwap(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    std::swap(f.mySummands, g.mySummands);
    std::swap(f.myEnd, g.myEnd);
    if (f.mySummands == nullptr) f.myEnd = &f.mySummands;  // CAREFUL if f or g is zero!
    if (g.mySummands == nullptr) g.myEnd = &g.mySummands;  //
  }


  void DistrMPolyInlFpPP::ourDeleteSummands(DistrMPolyInlFpPP::summand* ptr, MemPool& MemMgr)
  {
    while (ptr != nullptr)
    {
      DistrMPolyInlFpPP::summand* next = ptr->myNext;
      ptr->~summand();
//???      M.KillOrdv(curr->myOrdv);
      MemMgr.free(ptr);
      ptr = next;
    }
  }


  DistrMPolyInlFpPP::DistrMPolyInlFpPP(const InlineFpImpl& arith, const ring& R, const PPMonoid& PPM, const OrdvArith::reference& OA, MemPool& MemMgr):
      myILCoeffImpl(arith),
      myCoeffRing(R),
      myPPM(PPM),
      myOrdvArith(OA),
      mySummandMemory(MemMgr)
  {
    mySummands = nullptr;
    myEnd = &mySummands;
  }


  DistrMPolyInlFpPP::~DistrMPolyInlFpPP()
  {
    ourDeleteSummands(mySummands, /*myILCoeffImpl, myPPM,*/ mySummandMemory);
  }


  DistrMPolyInlFpPP::DistrMPolyInlFpPP(const DistrMPolyInlFpPP& copy):
      myILCoeffImpl(copy.myILCoeffImpl),
      myCoeffRing(copy.myCoeffRing),
      myPPM(copy.myPPM),
      myOrdvArith(copy.myOrdvArith),
      mySummandMemory(copy.mySummandMemory)
  {
    mySummands = nullptr;
    myEnd = &mySummands;
    // !!! THIS LOOP IS NOT EXCEPTION CLEAN !!!
    for (summand* it = copy.mySummands; it != nullptr; it = it->myNext)
      myPushBack(myCopySummand(it));
  }


  DistrMPolyInlFpPP& DistrMPolyInlFpPP::operator=(const DistrMPolyInlFpPP& rhs)
  {
    if (this == &rhs) return *this;
    DistrMPolyInlFpPP copy(rhs);
    ourSwap(*this, copy);
    return *this;
  }


  DistrMPolyInlFpPP& DistrMPolyInlFpPP::operator=(const MachineInt& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this; // to avoid needless alloc/free
    NewSummandPtr t(*this);
    t.myRenew();
    t->myCoeff = myILCoeffImpl.myReduce(rhs);
    if (!IsZero(t->myCoeff))
    {
      mySummands = grab(t);
      myEnd = &mySummands->myNext;
    }
    return *this;
  }


  DistrMPolyInlFpPP& DistrMPolyInlFpPP::operator=(const BigInt& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this; // to avoid needless alloc/free
    NewSummandPtr t(*this);
    t.myRenew();
    t->myCoeff = myILCoeffImpl.myReduce(rhs);
    if (!IsZero(t->myCoeff))
    {
      mySummands = grab(t);
      myEnd = &mySummands->myNext;
    }
    return *this;
  }

  DistrMPolyInlFpPP& DistrMPolyInlFpPP::operator=(const BigRat& rhs)
  {
    myAssignZero();
    if (IsZero(rhs)) return *this; // to avoid needless alloc/free
    NewSummandPtr t(*this);
    t.myRenew();
    t->myCoeff = myILCoeffImpl.myReduce(rhs);
    if (!IsZero(t->myCoeff))
    {
      mySummands = grab(t);
      myEnd = &mySummands->myNext;
    }
    return *this;
  }


//----------------------------------------------------------------------
// operations on summands
//----------------------------------------------------------------------


  DistrMPolyInlFpPP::summand* DistrMPolyInlFpPP::myCopySummand(const summand* original) const
  {
    NewSummandPtr copy(*this);
    copy.myRenew();

    copy->myCoeff = original->myCoeff;
    myOrdvArith->myAssign(copy->myOrdv, original->myOrdv);
    return grab(copy);
  }


// NEVER USED???
//    void DistrMPolyInlFpPP::SetSummandOrdv(summand* dest, const summand* src) const
//    {
//      size_t len = OrdvWords(ordering(myPPM));
//      for (size_t i=0; i < len; ++i)
//        dest->myOrdv[i] = src->myOrdv[i];
//    }


  void DistrMPolyInlFpPP::myAssignZero()
  {
    ourDeleteSummands(mySummands, /*myILCoeffImpl, myPPM,*/ mySummandMemory);
    mySummands = nullptr;
    myEnd = &mySummands;
  }


  bool DistrMPolyInlFpPP::myIsEqual(const summand* const lhs, const summand* const rhs) const
  {
    return (myOrdvArith->myCmp(lhs->myOrdv, rhs->myOrdv) == 0 &&
            (lhs->myCoeff == rhs->myCoeff));
  }


  long NumTerms(const DistrMPolyInlFpPP& f)
  {
    long nsummands = 0;
    for (const DistrMPolyInlFpPP::summand* it = f.mySummands; it != nullptr; it = it->myNext)
      ++nsummands;
    return nsummands;
  }


  const DistrMPolyInlFpPP::InlineFpElem_t& LC(const DistrMPolyInlFpPP& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    return f.mySummands->myCoeff;
  }


  void MoveLMToFront(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlFpPP::summand* ltg = g.mySummands;
    g.mySummands = g.mySummands->myNext;
    if (g.mySummands == nullptr) g.myEnd = &(g.mySummands);
    ltg->myNext = nullptr; // not really needed, but checked when debugging is active
    f.myPushFront(ltg);
  }


  void MoveLMToBack(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlFpPP::summand* ltg = g.mySummands;
    g.mySummands = g.mySummands->myNext;
    if (g.mySummands == nullptr) g.myEnd = &(g.mySummands);
    ltg->myNext = nullptr;
    f.myPushBack(ltg);
  }


  void DistrMPolyInlFpPP::myDeleteLM()
  {
    CoCoA_ASSERT(!IsZero(*this));

    DistrMPolyInlFpPP::summand* old_lm = mySummands;
    mySummands = old_lm->myNext;
    if (mySummands == nullptr) myEnd = &mySummands;
    old_lm->myNext = nullptr;
    ourDeleteSummands(old_lm, /*myILCoeffImpl, myPPM,*/ mySummandMemory);
  }


//   void wdeg(degree& d, const DistrMPolyInlFpPP& f)
//   {
//     CoCoA_ASSERT(!IsZero(f));
//     f.myOrdvArith->myWDeg(d, f.mySummands->myOrdv);
//   }


//   int CmpWDeg(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g)
//   {
//     CoCoA_ASSERT(!IsZero(f));
//     CoCoA_ASSERT(!IsZero(g));
//     return f.myOrdvArith->myCmpWDeg(f.mySummands->myOrdv, g.mySummands->myOrdv);
//   }


// This fn offers only the weak exception guarantee!!!
  void DistrMPolyInlFpPP::myAddMulSummand(const summand* s, const DistrMPolyInlFpPP& g, bool SkipLMg)  // this += s*g
  {
    CoCoA_ASSERT(IsCompatible(*this, g));

    const InlineFpImpl& Fp = myILCoeffImpl;

    const summand* g_smnd = g.mySummands;
    if (SkipLMg)    g_smnd = g_smnd->myNext;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;

    NewSummandPtr tmp_smnd(*this);
    tmp_smnd.myRenew();

    int CMP = 0;

    //    bool qIsOne = myOrdvArith->myIsZero(s->myOrdv);  // becomes slower!!

    for (; f_smnd != nullptr && g_smnd != nullptr; g_smnd = g_smnd->myNext)
    {
//       if (qIsOne)
//       {
//         while (f_smnd != nullptr && (CMP=myOrdvArith->myCmp(f_smnd->myOrdv, g_smnd->myOrdv)) >0)
//           f_smnd = *(f_prev = &f_smnd->myNext);
//         myOrdvArith->myAssign(tmp_smnd->myOrdv, g_smnd->myOrdv);
//       }
//       else
      {
        myOrdvArith->myMul(tmp_smnd->myOrdv, s->myOrdv, g_smnd->myOrdv);
        while (f_smnd != nullptr && (CMP=myOrdvArith->myCmp(f_smnd->myOrdv, tmp_smnd->myOrdv)) >0)
          f_smnd = *(f_prev = &f_smnd->myNext);
      }
      if (f_smnd == nullptr)
      {
        tmp_smnd->myCoeff = Fp.myMul(s->myCoeff, g_smnd->myCoeff);
        myPushBack(grab(tmp_smnd));
        tmp_smnd.myRenew();
        g_smnd = g_smnd->myNext;
        break;
      }
      if (CMP == 0)
      {
        if (Fp.myIsZeroAddMul(f_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff))
          myRemoveSummand(f_prev);  // f_prev = f_prev;
        else
          f_prev = &f_smnd->myNext;
        f_smnd = *f_prev;
      }
      else // (CMP < 0)
      {
        tmp_smnd->myCoeff = Fp.myMul(s->myCoeff, g_smnd->myCoeff);
        myInsertSummand(grab(tmp_smnd), f_prev);
        tmp_smnd.myRenew();
        f_prev = &(*f_prev)->myNext;
        // f_smnd = f_smnd;
      }
    }
    for (;g_smnd != nullptr; g_smnd = g_smnd->myNext)
    {
      myOrdvArith->myMul(tmp_smnd->myOrdv, s->myOrdv, g_smnd->myOrdv);
      tmp_smnd->myCoeff = Fp.myMul(s->myCoeff, g_smnd->myCoeff);
      myPushBack(grab(tmp_smnd));
      tmp_smnd.myRenew();
    }

//      clog << "AddMul: produced f=";output(clog,f);clog<<std::endl;
//      clog << "------------------------------------------------------"<<std::endl;

  }


  void DistrMPolyInlFpPP::myAddMulLM(const DistrMPolyInlFpPP& h, const DistrMPolyInlFpPP& g, bool SkipLMg)
  {                                                 //???
    if (IsZero(h))  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    myAddMulSummand(h.mySummands, g, SkipLMg);     //???
  }                                                 //???


//   void DistrMPolyInlFpPP::myWeylAddMulSummand(const summand* s, const DistrMPolyInlFpPP& g, bool SkipLMg)
//   {
//     CoCoA_ASSERT(IsCompatible(*this, g));

//     const InlineFpImpl& Fp = myILCoeffImpl;
//     const size_t nvars = NumIndets(myPPM);
// //???    clog << "AddMul: Doing funny product of the following two polys" << std::endl;
// //???    output(clog, g);
// //???    CoCoA_ASSERT(myRefCount == 1);
//     //????    MakeWritable(f);

//     const summand* g_term = g.mySummands;
//     if (SkipLMg)    g_term = g_term->myNext;
// //???    summand** f_prev = &mySummands;
// //???    summand*  f_term = *f_prev;

//     DistrMPolyInlFpPP ppg = g;
//     vector<long> expv(nvars);
//     myOrdvArith->myComputeExpv(expv, s->myOrdv);
// //???    clog << "expv: "; for (int i=0; i<myNumIndets;++i) clog << expv[i] << "  "; clog << std::endl;
//     for (size_t indet = nvars/2; indet < nvars; ++indet)
//     {
//       long n = expv[indet];
//       if (n == 0) continue;
// //???      clog << "AddMul: doing D variable with index " << indet - myNumIndets/2 << std::endl;
//       DistrMPolyInlFpPP der = ppg;

// //???      ppg *= IndetPower(myPPM, indet, n);
//       ppg.myMulByPP(raw(IndetPower(myPPM, indet, n)));

//       for (long i=1; i <= n; ++i)
//       {
//         deriv(der, der, indet-nvars/2);
// //???        deriv(raw(der), raw(der), indet-nvars/2);
// //???        clog << "der(" << i << ")="; output(clog, raw(der)); clog << std::endl;

// //        ppg += binomial(n, i)*der*IndetPower(myPPM, indet, n-i); // *IndetPower(myPPM, h, 2*i); // for homog case
//         NewSummandPtr jaa(*this);
//         jaa.myRenew();
//         Fp.myAssign(jaa->myCoeff, binomial(n, i));
//         myOrdvArith->myMulIndetPower(jaa->myOrdv, indet, n-i);
// // if (HOMOG_CASE)        myOrdvArith->mul(jaa->myOrdv, jaa->myOrdv, IndetPower(myPPM, h, 2*i));
//         ppg.myAddMulSummand(jaa.get(), der, false);
//       }
//     }
//     { // f *= indet^deg(pp, indet); for the early vars
//       for (size_t indet = nvars/2; indet < nvars; ++indet)
//         expv[indet] = 0;
//       NewSummandPtr jaa(*this);
//       jaa.myRenew();
//       Fp.myAssign(jaa->myCoeff, s->myCoeff);
// //???      myOrdvArith->assign(jaa->myOrdv, expv????);
//       myOrdvArith->myAssignFromExpv(jaa->myOrdv, expv);
//       myAddMulSummand(jaa.get(), ppg, false);
//     }
//   }




  void DistrMPolyInlFpPP::myReductionStep(const DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(&g!=this);
    CoCoA_ASSERT(mySummands != nullptr);
    CoCoA_ASSERT(!IsZero(g));

    DistrMPolyInlFpPP tmp_poly(myILCoeffImpl, myCoeffRing, myPPM, myOrdvArith, mySummandMemory);

    DivLM(tmp_poly, *this, g);
    tmp_poly.myNegate();
    myDeleteLM();
    myAddMulLM(tmp_poly, g, /*SkipLMg = */ true );
  }


//   static void ComputeFScaleAndGScale(const ring& R,
//                                      const InlineFpImpl::RawValue& LCf,
//                                      const InlineFpImpl::RawValue& LCg,
//                                      InlineFpImpl::RawValue& fscale,
//                                      InlineFpImpl::RawValue& gscale)
//   {
//     CoCoA_THROW_ERROR("It does not make sense", "DistrMPolyInlFpPP::ComputeFScaleAndGScale");
//   }


  void DistrMPolyInlFpPP::myReductionStepGCD(const DistrMPolyInlFpPP& /*g*/, RingElem& /*FScale*/)
  { CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere); }


  void DistrMPolyInlFpPP::myAddClear(DistrMPolyInlFpPP& g) // sets g to 0 as a side-effect
  {
    const InlineFpImpl& Fp = myILCoeffImpl;
    typedef DistrMPolyInlFpPP::summand summand;

    summand*  g_smnd = g.mySummands;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;
    int CMP=0;
//???????    CoCoA_ASSERT(*(G.myEnd)==nullptr);//BUG HUNTING  ???

   //    clog << "input f = "; output(clog, *this) ;clog << std::endl;
    while ( f_smnd!=nullptr && g_smnd!=nullptr )
    {
      while (f_smnd!=nullptr &&
             (CMP = myOrdvArith->myCmp(f_smnd->myOrdv, g_smnd->myOrdv)) >0)
        f_smnd = *(f_prev = &f_smnd->myNext);
      if (f_smnd == nullptr)  break;
      //clog <<   "(AddClear error: should never happen for Basic Reduction)" << std::endl;
      g.mySummands = g.mySummands->myNext;
      g_smnd->myNext = nullptr;
      if (CMP == 0)
      {
        f_smnd->myCoeff = Fp.myAdd(f_smnd->myCoeff, g_smnd->myCoeff);
        if (IsZero(f_smnd->myCoeff))
          myRemoveSummand(f_prev);
        ourDeleteSummands(g_smnd, /*myILCoeffImpl, myPPM,*/ mySummandMemory);
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


  void DistrMPolyInlFpPP::myAppendClear(DistrMPolyInlFpPP& g)  // sets g to 0; no throw guarantee!
  {
    if (g.mySummands == nullptr) return;
    *(myEnd) = g.mySummands;
    myEnd = g.myEnd;
    g.mySummands = nullptr;
    g.myEnd = &g.mySummands;
  }


  ConstRefPPMonoidElem LPP(const DistrMPolyInlFpPP& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    return ConstRefPPMonoidElem(f.myPPM, PPMonoidElemConstRawPtr(f.mySummands->myOrdv));
  }


  void DivLM(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g) // lhs = LM(f)/LM(g)
  {
    CoCoA_ASSERT(IsCompatible(f, g) && IsCompatible(lhs, f));
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    //    clog << "DivLM" << std::endl;
    const DistrMPolyInlFpPP::InlineFpImpl& Fp = f.myILCoeffImpl;    // shorthand
    typedef DistrMPolyInlFpPP::summand summand;            // shorthand
    const summand* const LMf = f.mySummands;  // shorthand
    const summand* const LMg = g.mySummands;  // shorthand

    CoCoA_ASSERT(IsDivisible(LPP(f), LPP(g)));
    DistrMPolyInlFpPP::NewSummandPtr SpareSummand(lhs);
    SpareSummand.myRenew();
    SpareSummand->myCoeff = Fp.myDiv(LMf->myCoeff, LMg->myCoeff);
    f.myOrdvArith->myDiv(SpareSummand->myOrdv, LMf->myOrdv, LMg->myOrdv);
    lhs.myAssignZero();  // CANNOT myAssignZero() EARLIER in case lhs aliases f or g.
    lhs.myPushBack(grab(SpareSummand));
  }


  //??? THIS IS WRONG IF THERE ARE zero-divisors
  // BUT do we allow zero-divisors in this case???
  void DistrMPolyInlFpPP::myMulByCoeff(const InlineFpElem_t c)
  {
    CoCoA_ASSERT(!IsZero(c));
    for (summand* it = mySummands; it != nullptr; it = it->myNext)
      it->myCoeff = myILCoeffImpl.myMul(it->myCoeff, c);
  }


  bool DistrMPolyInlFpPP::myDivByCoeff(const InlineFpElem_t c)
  {
    if (IsZero(c)) return false;
    // modified by JAA 2013-03-24
    myMulByCoeff(myILCoeffImpl.myRecip(c));
    return true;
//     for (summand* it = mySummands; it != nullptr; it = it->myNext)
//       it->myCoeff = myILCoeffImpl.myDiv(it->myCoeff, c);
//     return true;
  }


  // ANNA: this can be improved... is it worth?
  void DistrMPolyInlFpPP::myMulByPP(PPMonoidElemConstRawPtr rawpp)
  {
    vector<OrdvArith::OrdvElem> ordv(OrdvWords(myOrdvArith));
    {
      vector<long> expv(NumIndets(myPPM));
      myPPM->myExponents(expv, rawpp);
      myOrdvArith->myAssignFromExpv(&ordv[0], expv);
    }

    for (summand* f_smnd = mySummands; f_smnd != nullptr ; f_smnd = f_smnd->myNext)
      myOrdvArith->myMul(f_smnd->myOrdv, f_smnd->myOrdv, &ordv[0]);
  }

// ??? WANT DIVISION BY A PP TOO???


//   void DistrMPolyInlFpPP::myWeylMul(PPMonoidElemConstRawPtr rawpp)
//   {
//     vector<OrdvArith::OrdvElem> ordv(OrdvWords(myOrdvArith));
//     {
//       vector<long> expv(NumIndets(myPPM));
//       myPPM->myExponents(expv, rawpp);
//       myOrdvArith->myAssignFromExpv(&ordv[0], expv);
//     }

//     for (summand* f_smnd = mySummands; f_smnd != nullptr ; f_smnd = f_smnd->myNext)
//       myOrdvArith->myMul(f_smnd->myOrdv, f_smnd->myOrdv, &ordv[0]);
//   }


  void DistrMPolyInlFpPP::myPushFront(const InlineFpElem_t c, const std::vector<long>& expv)
  {
    if (IsZero(c)) return;

    NewSummandPtr tmp(*this);
    tmp.myRenew();
    tmp->myCoeff = c;///???     tmp->myCoeff =  myILCoeffImpl.myReduce(c);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushFront(grab(tmp)); 
  }


  void DistrMPolyInlFpPP::myPushBack(const InlineFpElem_t c, const std::vector<long>& expv)
  {
    if (IsZero(c)) return;

    NewSummandPtr tmp(*this);
    tmp.myRenew();
    tmp->myCoeff = c;///???     tmp->myCoeff =  myILCoeffImpl.myReduce(c);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushBack(grab(tmp));
  }


  void DistrMPolyInlFpPP::myPushFront(const InlineFpElem_t c, PPMonoidElemConstRawPtr rawpp)
  {
    if (IsZero(c)) return;

    NewSummandPtr tmp(*this);
    tmp.myRenew();
    tmp->myCoeff = c;
    //    myOrdvArith->myAssign(tmp->myOrdv, rawpp);
    vector<long> expv(NumIndets(myPPM));
    myPPM->myExponents(expv, rawpp);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushFront(grab(tmp));
  }


  void DistrMPolyInlFpPP::myPushBack(const InlineFpElem_t c, PPMonoidElemConstRawPtr rawpp)
  {
    if (IsZero(c)) return;

    NewSummandPtr tmp(*this);
    tmp.myRenew();
    tmp->myCoeff = c;
    //    myOrdvArith->myAssign(tmp->myOrdv, rawpp);
    vector<long> expv(NumIndets(myPPM));
    myPPM->myExponents(expv, rawpp);
    myOrdvArith->myAssignFromExpv(tmp->myOrdv, expv);
    myPushBack(grab(tmp));
  }


  void DistrMPolyInlFpPP::myPushFront(summand* t)
  {
    CoCoA_ASSERT(t->myNext == nullptr);
    CoCoA_ASSERT(IsZero(*this) || myOrdvArith->myCmp(t->myOrdv, mySummands->myOrdv) > 0);
    t->myNext = mySummands;
    mySummands = t;
    if (myEnd == &mySummands) myEnd = &t->myNext;
  }


  void DistrMPolyInlFpPP::myPushBack(summand* t)
  {
    CoCoA_ASSERT(t->myNext == nullptr);
#ifdef CoCoA__DEBUG
    // Costly check that t really is smaller than smallest PP in the DistrMPolyInlFpPP.
    if (!IsZero(*this))
    {
      summand* LastSummand = mySummands;
      while (LastSummand->myNext != nullptr)
        LastSummand = LastSummand->myNext;
      CoCoA_ASSERT(myOrdvArith->myCmp(t->myOrdv, LastSummand->myOrdv) < 0);
    }
#endif
    *myEnd = t;
    myEnd = &t->myNext;
  }


  void DistrMPolyInlFpPP::myRemoveSummand(summand** prev_link)
  {
    summand* DeleteMe = *prev_link;
    CoCoA_ASSERT(DeleteMe != nullptr);
    if (DeleteMe->myNext == nullptr)
      myEnd = prev_link;

    *prev_link = DeleteMe->myNext;
    DeleteMe->myNext = nullptr;
    ourDeleteSummands(DeleteMe, /*myILCoeffImpl, myPPM,*/ mySummandMemory);
  }


  void DistrMPolyInlFpPP::myInsertSummand(summand* s, summand** prev_link)
  {
    s->myNext = (*prev_link);
    (*prev_link) = s;
    if (myEnd == prev_link) myEnd = &(s->myNext);
  }


  bool IsZeroAddLCs(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    CoCoA_ASSERT( LPP(f) == LPP(g) );
    f.mySummands->myCoeff = f.myILCoeffImpl.myAdd(f.mySummands->myCoeff, g.mySummands->myCoeff);
    g.myDeleteLM();
    if (!IsZero(f.mySummands->myCoeff)) return false;
    f.myDeleteLM();
    return true;
  }


  void DistrMPolyInlFpPP::myNegate()
  {
    for (summand* iter = mySummands; iter != nullptr; iter = iter->myNext)
      iter->myCoeff = myILCoeffImpl.myNegate(iter->myCoeff);
  }


  void add(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h)
  {
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(g, h));
    const DistrMPolyInlFpPP::InlineFpImpl& Fp = lhs.myILCoeffImpl;
    DistrMPolyInlFpPP ans(lhs.myILCoeffImpl, lhs.myCoeffRing, lhs.myPPM, lhs.myOrdvArith, lhs.mySummandMemory);

    typedef DistrMPolyInlFpPP::summand summand;
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

    DistrMPolyInlFpPP::NewSummandPtr SpareSummand(ans);
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
      SpareSummand->myCoeff = Fp.myAdd(gterm->myCoeff, hterm->myCoeff);
      if (!IsZero(SpareSummand->myCoeff))
      {
	lhs.myOrdvArith->myAssign(SpareSummand->myOrdv, gterm->myOrdv); // set PP ordv
	ans.myPushBack(grab(SpareSummand));
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
  void DistrMPolyInlFpPP::myAddMonomial(const DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(IsCompatible(*this, g));
    CoCoA_ASSERT(NumTerms(g)==1);

    typedef DistrMPolyInlFpPP::summand summand;
    summand** f_prev = &mySummands;
    summand*  f_smnd = *f_prev;
    DistrMPolyInlFpPP::NewSummandPtr NewTerm(*this);
    NewTerm.myRenew();
    NewTerm->myCoeff = (g.mySummands)->myCoeff;
    myOrdvArith->myAssign(NewTerm->myOrdv, (g.mySummands)->myOrdv);

    int CMP;
    while (f_smnd != nullptr &&
           (CMP = myOrdvArith->myCmp(f_smnd->myOrdv, NewTerm->myOrdv)) > 0)
      f_smnd = *(f_prev = &f_smnd->myNext);

    if (f_smnd == nullptr) { myPushBack(grab(NewTerm)); return; }

    if (CMP < 0) // g comes before f_smnd, so insert it
    {
      myInsertSummand(grab(NewTerm), f_prev);
//?????JAA: USELESS???        f_prev = &(*f_prev)->myNext;
      return;
    }

    // PPs are the same, so add coeffs
    NewTerm->myCoeff = myILCoeffImpl.myAdd(NewTerm->myCoeff, f_smnd->myCoeff);
    if (IsZero(NewTerm->myCoeff))
      myRemoveSummand(f_prev);
    else 
      f_smnd->myCoeff = NewTerm->myCoeff;
  }


  void sub(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h)
  {
    // This code is almost a copy of add(...).
    CoCoA_ASSERT(IsCompatible(lhs, g) && IsCompatible(g, h));
    const DistrMPolyInlFpPP::InlineFpImpl& Fp = lhs.myILCoeffImpl;
    DistrMPolyInlFpPP ans(lhs.myILCoeffImpl, lhs.myCoeffRing, lhs.myPPM, lhs.myOrdvArith, lhs.mySummandMemory);

    typedef DistrMPolyInlFpPP::summand summand;
    const summand* gterm = g.mySummands;
    const summand* hterm = h.mySummands;
    DistrMPolyInlFpPP::NewSummandPtr SpareSummand(ans);
    SpareSummand.myRenew();
    while (gterm != nullptr && hterm != nullptr)
    {
      int ord = lhs.myOrdvArith->myCmp(gterm->myOrdv, hterm->myOrdv);

      if (ord < 0)
      {
	summand* hcopy = ans.myCopySummand(hterm);
	hcopy->myCoeff = Fp.myNegate(hcopy->myCoeff);
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
      SpareSummand->myCoeff = Fp.mySub(gterm->myCoeff, hterm->myCoeff);
      if (!IsZero(SpareSummand->myCoeff))
      {
	lhs.myOrdvArith->myAssign(SpareSummand->myOrdv, gterm->myOrdv); // set PP ordv
	ans.myPushBack(grab(SpareSummand));
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
      hcopy->myCoeff = Fp.myNegate(hcopy->myCoeff);
      ans.myPushBack(hcopy);
      hterm = hterm->myNext;
    }

    swap(lhs, ans); // really an assignment
  }


  bool div(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h)  // result is true iff quotient is exact.
  {
    CoCoA_ASSERT(!IsZero(h));
    PPMonoid PPM = lhs.myPPM;
    const PPOrdering ord = ordering(PPM);
    const DistrMPolyInlFpPP::InlineFpImpl& Fp = lhs.myILCoeffImpl;
    const DistrMPolyInlFpPP::summand* LMh = h.mySummands;
    const PPMonoidElem LPPden(LPP(h));
    const SmallFpImpl::value RecipLCh = Fp.myRecip(LMh->myCoeff);
    DistrMPolyInlFpPP ans(lhs.myILCoeffImpl, lhs.myCoeffRing, lhs.myPPM, lhs.myOrdvArith, lhs.mySummandMemory);
    DistrMPolyInlFpPP dividend(g);
    while (!IsZero(dividend))
    {
      const DistrMPolyInlFpPP::summand* LMdividend = dividend.mySummands;
      DistrMPolyInlFpPP::NewSummandPtr qterm(lhs);
      if (!IsDivisible(LPP(dividend), LPPden)) return false;
      qterm.myRenew();
      qterm->myCoeff = Fp.myMul(LMdividend->myCoeff, RecipLCh); // JAA 2013-03-24
//      qterm->myCoeff = Fp.myDiv(LMdividend->myCoeff, LMh->myCoeff);
      g.myOrdvArith->myDiv(qterm->myOrdv, LMdividend->myOrdv, LMh->myOrdv);  //??? check whether this overflows?
      qterm->myCoeff = Fp.myNegate(qterm->myCoeff);
      dividend.myAddMulSummand(qterm.get(), h, false);
      qterm->myCoeff = Fp.myNegate(qterm->myCoeff);
      ans.myPushBack(grab(qterm));
    }
    swap(lhs, ans); // really an assignment
    return true;
  }


  void output(std::ostream& out, const DistrMPolyInlFpPP& f)  // for debugging only
  {
    if (!out) return;  // short-cut for bad ostreams

    if (IsZero(f)) { out << "0"; return; }
    const DistrMPolyInlFpPP::InlineFpImpl& Fp = f.myILCoeffImpl;
    const PPMonoid PPM = f.myPPM;
    for (DistrMPolyInlFpPP::summand* it = f.mySummands; it != nullptr; it = it->myNext)
    {
      out << " +(" << Fp.myExport(it->myCoeff) << ")*"
          << ConstRefPPMonoidElem(PPM, PPMonoidElemConstRawPtr(it->myOrdv));
    }
  }


//   bool IsConstant(const DistrMPolyInlFpPP& f)
//   {
//     if (IsZero(f)) return true;
//     if (!IsMonomial(f)) return false;
//     return IsOne(LPP(f));
//   }


  bool IsZero(const DistrMPolyInlFpPP& f)
  {
    return (f.mySummands == nullptr);
  }


//   bool IsOne(const DistrMPolyInlFpPP& f)
//   {
//     if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
//     if (!f.myOrdvArith->myIsZero(f.mySummands->myOrdv)) return false;
//     return f.myILCoeffImpl.myIsOne(f.mySummands->myCoeff);
//   }


//   bool IsMinusOne(const DistrMPolyInlFpPP& f)
//   {
//     if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
//     if (!f.myOrdvArith->myIsZero(f.mySummands->myOrdv)) return false;
//     return f.myILCoeffImpl.myIsMinusOne(f.mySummands->myCoeff);
//   }


//   bool IsConstant(const DistrMPolyInlFpPP& f)
//   {
//     if (IsZero(f)) return true;
//     if (f.mySummands->myNext != nullptr) return false;
//     return f.myOrdvArith->myIsZero(f.mySummands->myOrdv);
//   }


//   bool IsIndet(std::size_t& index, const DistrMPolyInlFpPP& f)
//   {
//     if (IsZero(f)) return false;
//     if (f.mySummands->myNext != nullptr) return false;
//     if (!f.myILCoeffImpl.myIsOne(f.mySummands->myCoeff)) return false;
//     return f.myOrdvArith->myIsIndet(index, f.mySummands->myOrdv);
//   }


  bool IsMonomial(const DistrMPolyInlFpPP& f)
  {
    if (IsZero(f) || f.mySummands->myNext != nullptr) return false;
    return true;
  }


  bool IsEqual(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(IsCompatible(f, g));
    if (&f == &g) return true;
    const DistrMPolyInlFpPP::summand* fterm = f.mySummands;
    const DistrMPolyInlFpPP::summand* gterm = g.mySummands;
    while (fterm != nullptr && gterm != nullptr)
    {
      if (!f.myIsEqual(fterm, gterm)) return false;
      fterm = fterm->myNext;
      gterm = gterm->myNext;
    }
    return fterm == gterm; // either both are nullptr (when the polys are equal), or only one is nullptr
  }


//   bool IsEqual(const DistrMPolyInlFpPP& f, long n)
//   {
//     if (n == 0) return IsZero(f);
//     if (IsZero(f)) return IsZero(RingElem(f.myILCoeffImpl, n));
//     // From here on the polynomial is known to be non-zero
//     if (f.mySummands->myNext != nullptr) return false;
//     if (!f.myOrdvArith->myIsZero(f.mySummands->myOrdv)) return false;
//     return f.myILCoeffImpl.myIsEqual(f.mySummands->myCoeff, n);
//   }


//   void WeylMul(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h)
//   {

//   }


//   void WeylDiv(DistrMPolyInlFpPP& /*lhs*/, const DistrMPolyInlFpPP& /*g*/, const DistrMPolyInlFpPP& /*h*/)
//   {
//   }


//   void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, long IndetIndex)
//   {
//     deriv(lhs, f, indet(f.myPPM, IndetIndex));
//   }


//   void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, ConstRefPPMonoidElem x)
//   {
//     if (IsOne(x)) { lhs = f; return; }
//     const long nvars = NumIndets(owner(x));
// //???    const PPOrdering ord = ordering(owner(x));
//     const DistrMPolyInlFpPP::InlineFpImpl& Fp = f.myILCoeffImpl;
//     vector<long> expv(nvars);
//     exponents(expv, x);
//     vector<OrdvArith::OrdvElem> ordvx(OrdvWords(f.myOrdvArith));
//     f.myOrdvArith->myAssignFromExpv(&ordvx[0], expv);
// //clog<<"differentiating wrt expv: [";for(long i=0;i<nvars;++i)clog<<expv[i]<<" ";clog<<"]"<<std::endl;
//     DistrMPolyInlFpPP ans(f.myILCoeffImpl, f.myCoeffRing, f.myPPM, f.myOrdvArith, f.mySummandMemory);

//     for (const DistrMPolyInlFpPP::summand* f_term = f.mySummands; f_term != nullptr; f_term = f_term->myNext)
//     {
// //clog<<"LOOPHEAD\n";
//       BigInt scale(1);
//       for (long indet=0; indet < nvars; ++indet)
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
//       DistrMPolyInlFpPP::NewSummandPtr tmp(f);
//       tmp.myRenew();
//       Fp.myAssign(tmp->myCoeff, scale);
//       Fp.myMul(tmp->myCoeff, tmp->myCoeff, f_term->myCoeff);
//       if (Fp.myIsZero(tmp->myCoeff)) continue;
// //clog<<"dividing ordv [";for(long i=0;i<2;++i)clog<<f_term->myOrdv[i]<<" ";clog<<"]\n";
// //clog<<"by       ordv [";for(long i=0;i<2;++i)clog<<ordvx[i]<<" ";clog<<"]\n";
//       f.myOrdvArith->myDiv(tmp->myOrdv, f_term->myOrdv, &ordvx[0]);
// //clog<<"Quotient is   [";for(long i=0;i<2;++i)clog<<tmp->myOrdv[i]<<" ";clog<<"]\n";
//       ans.myPushBack(grab(tmp));
//     }
//     swap(lhs, ans); // really an assignment
//   }



} // end of namespace CoCoA


// Source code for class DistrMPolyInlFpPPtMPolyInlPP

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/DistrMPolyInlFpPP.C,v 1.40 2024/08/02 16:09:14 bigatti Exp $
// $Log: DistrMPolyInlFpPP.C,v $
// Revision 1.40  2024/08/02 16:09:14  bigatti
// Summary: using CoCoA_THROW_ERROR1/2
//
// Revision 1.39  2024/07/02 09:57:57  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.38  2022/02/18 14:11:53  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.37  2020/06/17 15:49:22  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.36  2020/02/11 16:56:40  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.35  2020/02/11 16:12:17  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.34  2019/10/15 11:54:08  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.33  2019/03/19 11:07:07  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.32  2018/05/22 14:16:39  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.31  2018/05/18 12:13:36  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.30  2017/09/06 11:56:28  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.29  2015/11/04 10:32:16  abbott
// Summary: Consequential changes (after revising SmallFpImpl)
//
// Revision 1.28  2015/04/27 10:31:41  bigatti
// -- commented out excessive assert in myPushFront
//
// Revision 1.27  2015/04/27 10:06:04  bigatti
// -- simplified myAppendClear
// -- fixed myMoveLMToBack
// -- added CoCoA_ASSERT to myMoveLMToBack and myMoveLMToFront
//
// Revision 1.26  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.25  2015/04/17 16:24:07  abbott
// Summary: Added check that first arg is non-zero
// Author: JAA
//
// Revision 1.24  2014/11/18 16:07:08  abbott
// Summary: Replaced calls ModP.myDiv(1,x) by ModP.myRecip(x)
// Author: JAA
//
// Revision 1.23  2014/04/30 16:05:40  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.22  2014/01/28 10:57:38  bigatti
// -- removed useless ==1 on boolean
//
// Revision 1.21  2013/03/25 17:26:30  abbott
// Changed return type of NumTerms (now long, was size_t).
// Cleaned impl of 2 fns.
// There is A LOT of cruft here -- major cleaning needed!!
//
// Revision 1.20  2012/10/16 10:28:18  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.19  2012/10/05 15:35:43  bigatti
// -- added myAddMonomial
//
// Revision 1.18  2012/10/02 15:28:52  abbott
// Updated two assertions (they wanted to the use the obsolete myIsZero fn)
//
// Revision 1.17  2012/09/26 12:27:51  abbott
// Updated to new SmallFpImpl interface.
//
// Revision 1.16  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.15  2011/11/09 14:03:59  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.14  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.13  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.12  2011/05/20 19:26:05  abbott
// Updated SmallFp*Impl: removed all output-related fns (must use myExport instead).
//
// Revision 1.11  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.10  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.9  2009/09/28 17:14:41  bigatti
// -- commented out unused functions (div, deriv, *Weyl*)
//
// Revision 1.8  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.7  2008/04/10 15:15:32  bigatti
// -- added  void myPushFront(rawc, rawpp)
//
// Revision 1.6  2007/12/05 12:11:07  bigatti
// -- cleaning (mostly removing unused code)
//
// Revision 1.5  2007/12/04 14:27:07  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.4  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
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
// Revision 1.15  2007/03/08 18:22:30  cocoa
// Just whitespace cleaning.
//
// Revision 1.14  2007/03/07 13:42:45  bigatti
// -- removed useless argument and other minor changes
//
// Revision 1.13  2007/01/15 13:34:30  cocoa
// -- added prefix "raw" to RawPtr arguments names
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
// Revision 1.9  2006/11/21 18:09:24  cocoa
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
// Revision 1.11  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.10  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.9  2006/03/30 16:59:27  cocoa
// -- changed misleading name: InlineCoeffRing --> InlineCoeffImpl
// -- new: implementation for homomorphisms
// -- rearrangement of code to mimic RingDistrMPolyInlPP
//
// Revision 1.8  2006/03/27 12:21:25  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.7  2006/03/20 17:27:42  cocoa
// -- changed in DistrMPolyInlFpPP: myMul, myDiv --> myMulByCoeff, myMulByPP, myDivByCoeff
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
// Revision 1.5  2005/07/08 15:09:29  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.4  2005/07/01 16:09:14  cocoa
// Degrees may now have negative components.
//
// Revision 1.3  2005/06/30 16:09:58  cocoa
// -- in "div" arguments of IsDivisible were swapped
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
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.9  2004/11/11 13:22:24  cocoa
// -- minor changes for doxygen
// -- changed: cout --> GlobalLogput()
//
// Revision 1.8  2004/11/08 13:56:02  cocoa
// -- changed calls to ZZ (after changes to ZZ.H)
//
// Revision 1.7  2004/11/02 18:21:21  cocoa
// -- changed: myGetExpvBuffer --> myExpvBufferRef
//
// Revision 1.6  2004/11/02 15:52:06  cocoa
// -- changed LPP body: now it uses myGetExpvBuffer
//
// Revision 1.5  2004/10/29 15:26:18  cocoa
// -- code fixed for compatibility with OrdvArith
//
// Revision 1.3  2004/07/20 15:37:08  cocoa
// Minor fix for some errors which slipped through the net...
//
// Revision 1.2  2004/07/16 10:11:34  cocoa
// -- now using the new class SmallFpImpl (or SmallFpLogImpl)
// -- updated with "my" coding convenctions
// -- NYI: LC and LCRaw
//
// Revision 1.1  2004/06/25 16:03:58  cocoa
// -- first import
//

