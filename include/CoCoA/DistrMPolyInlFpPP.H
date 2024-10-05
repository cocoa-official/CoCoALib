#ifndef CoCoA_DistrMPolyInlFpPP_H
#define CoCoA_DistrMPolyInlFpPP_H

//   Copyright (c)  2005  Anna Bigatti

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


// Header file for the class DistrMPolyInlFpPP


#include "CoCoA/assert.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/OrdvArith.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SmallFpLogImpl.H"


// #include <cstddef> --- already included in MemPool.H
// using std::size_t;
#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;


namespace CoCoA
{

/*-- class DistrMPolyInlFpPP --------------------------------------*/
/**
  \brief Implementation of polynomials for RingDistrMPolyInlFpPP

  \include DistrMPolyInlFpPP.txt
  */
/*-----------------------------------------------------------------*/

  class DistrMPolyInlFpPP
  {
  public:
    typedef SmallFpImpl InlineFpImpl;
    //typedef SmallFpLogImpl InlineFpImpl;
    typedef InlineFpImpl::value  InlineFpElem_t;

  private:

    //    class NewSummandPtr;
    friend class NewSummandPtr;

    struct summand
    {
    public:                    // These two lines are here
      static void useless() {} // just so gcc doesn't print annoying warnings.
      summand(): myNext(nullptr), myCoeff(zero(SmallFp)) {};
      // default dtor is OK
    public:
      summand* myNext;
      InlineFpElem_t myCoeff;
      OrdvArith::OrdvElem myOrdv[1]; // ASSUMES the compiler will allocate this at the end,
                                     // so I can address elements with index greater than
                                     // 1 with impunity if I know the space exists.
    public: // disable copy ctor & assignment
      summand(const summand&) = delete;
      summand& operator=(const summand&) = delete;
    };
    static void ourDeleteSummands(summand* ptr, MemPool& MemMgr);

    // This class should permit exception-clean implementations.
    class NewSummandPtr
    {
    public:
      explicit NewSummandPtr(const DistrMPolyInlFpPP& f): myPtr(nullptr), /*myFp(f.myILCoeffImpl),*/ /*myPPM(f.myPPM),*/ myOrdvArith(f.myOrdvArith.myRawPtr()), myMemMgr(f.mySummandMemory) { }
      ~NewSummandPtr() { if (myPtr == nullptr) return; myPtr->~summand(); myMemMgr.free(myPtr); }
      const summand* get() const { return myPtr; }
      summand* operator->() { return myPtr; }
      const summand* operator->() const { return myPtr; }
      summand* relinquish()  { CoCoA_ASSERT(myPtr != nullptr); summand* ans = myPtr; myPtr = nullptr; return ans; }
      void myRenew() { CoCoA_ASSERT(myPtr == nullptr); alloc(); myAssignZero(); }
      bool myCoeffIsZero() { CoCoA_ASSERT(myPtr != nullptr); return IsZero(myPtr->myCoeff); }
    private: // data members
      summand* myPtr; // either nullptr or points to a summand
/*      const InlineFpImpl& myFp;*/
//???      const PPMonoid myPPM;
      const OrdvArith::base* myOrdvArith;
      MemPool& myMemMgr;
      void alloc() { CoCoA_ASSERT(myPtr == nullptr); myPtr = static_cast<summand*>(myMemMgr.alloc()); }
      void myAssignZero() { new(myPtr) summand(); myOrdvArith->myAssignZero(myPtr->myOrdv); };

    public: // disable copy ctor & assignment
      NewSummandPtr(const NewSummandPtr&) = delete;
      NewSummandPtr& operator=(const NewSummandPtr&) = delete;
    };
    friend inline summand* grab(NewSummandPtr& s) { CoCoA_ASSERT(!s.myCoeffIsZero()); return s.relinquish(); }


  private: // data members
    const InlineFpImpl& myILCoeffImpl;
    const ring myCoeffRing;
    const PPMonoid myPPM;
    OrdvArith::reference myOrdvArith;
    MemPool& mySummandMemory; ///< the MemPool which manages the memory occupied by the summands
    summand** myEnd;          ///< pointer to facilitate implementing the one-past-the-last iterator
    summand* mySummands;      ///< pointer to the first summand

  public:
    DistrMPolyInlFpPP(const InlineFpImpl&, const ring&, const PPMonoid&, const OrdvArith::reference&, MemPool&);
    DistrMPolyInlFpPP(const DistrMPolyInlFpPP& copy);
    ~DistrMPolyInlFpPP();
    DistrMPolyInlFpPP& operator=(const DistrMPolyInlFpPP& rhs);
    DistrMPolyInlFpPP& operator=(const MachineInt& rhs);
    DistrMPolyInlFpPP& operator=(const BigInt& rhs);
    DistrMPolyInlFpPP& operator=(const BigRat& rhs);

    friend bool IsCompatible(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
    static std::size_t SummandSize(const ring&, const OrdvArith::reference& OA);
    static void ourSwap(DistrMPolyInlFpPP&, DistrMPolyInlFpPP&);
    summand* myCopySummand(const summand* original) const;
    // functions depending solely on the structure
    void myAssignZero();

    bool myIsEqual(const summand* const lhs, const summand*const rhs) const;
    friend long NumTerms(const DistrMPolyInlFpPP& f);
    friend const InlineFpElem_t& LC(const DistrMPolyInlFpPP& f);
    friend ConstRefPPMonoidElem LPP(const DistrMPolyInlFpPP& f);
    friend int CmpLPP(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
    friend void MoveLMToFront(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g);
    friend void MoveLMToBack(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g);
    void myDeleteLM();
//     friend void wdeg(degree& d, const DistrMPolyInlFpPP& f);
//     friend int CmpWDeg(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
    void myAddMulLM(const DistrMPolyInlFpPP& h, const DistrMPolyInlFpPP& g, bool SkipLMg);
//    void myWeylAddMul(const DistrMPolyInlFpPP& h, const DistrMPolyInlFpPP& g);
    void myReductionStep(const DistrMPolyInlFpPP& g);
    void myReductionStepGCD(const DistrMPolyInlFpPP& g, RingElem& fscale);
    void myAddClear(DistrMPolyInlFpPP& g);
    void myAppendClear(DistrMPolyInlFpPP& g);
    friend void DivLM(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
    void myMulByCoeff(InlineFpElem_t c);  ///< NO THROW GUARANTEE (???)
    bool myDivByCoeff(InlineFpElem_t c);  ///< NO THROW GUARANTEE (???)
    void myMulByPP(PPMonoidElemConstRawPtr rawpp);  ///< EXCEPTION SAFE
    //    void myWeylMul(PPMonoidElemConstRawPtr rawpp);
    void myPushFront(const InlineFpElem_t c, const std::vector<long>& expv);
    void myPushBack(const InlineFpElem_t c, const std::vector<long>& expv);
    void myPushFront(const InlineFpElem_t c, PPMonoidElemConstRawPtr rawpp);
    void myPushBack(const InlineFpElem_t c, PPMonoidElemConstRawPtr rawpp);
    void myPushFront(summand* t);
    void myPushBack(summand* t);
    void myAddMonomial(const DistrMPolyInlFpPP& g); ///< EXCEPTION SAFE
    void myRemoveSummand(summand** prev_link);
    void myInsertSummand(summand* s, summand** prev_link);
    friend bool IsZeroAddLCs(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g);
    void myNegate();

    friend void add(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
    friend void sub(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
    friend bool div(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
    friend void output(std::ostream& out, const DistrMPolyInlFpPP& f);
    friend bool IsZero(const DistrMPolyInlFpPP& f);
//     friend bool IsOne(const DistrMPolyInlFpPP& f);
//     friend bool IsMinusOne(const DistrMPolyInlFpPP& f);
//     friend bool IsConstant(const DistrMPolyInlFpPP& f);
//     friend bool IsIndet(long& index, const DistrMPolyInlFpPP& f);
    friend bool IsMonomial(const DistrMPolyInlFpPP& f);
    friend bool IsEqual(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
//     friend void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP&, long IndetIndex);
//     friend void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, ConstRefPPMonoidElem x);

//     friend void WeylMul(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
//     friend void WeylDiv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
  private:
     void myAddMulSummand(const summand* s, const DistrMPolyInlFpPP& g, bool SkipLMg);
    //     void myWeylAddMulSummand(const summand* s, const DistrMPolyInlFpPP& g, bool SkipLMg);

  public:
    // iter implements a SparsePolyIter for DistrMPolyInlPP objects
    class iter: public SparsePolyIterBase
    {
    public:
      explicit iter(const DistrMPolyInlFpPP& f);
      iter(const DistrMPolyInlFpPP& f, const summand* TermPtr);
      // Default copy ctor and assignment are fine.
      // Default dtor is fine.
      iter* myClone() const override;
      bool IamEnded() const override;
      void myForward() override;
      iter& operator++() /*override*/;
      iter operator++(int) /*override*/;
      RingElemAlias myCoeff() const override;
      ConstRefPPMonoidElem myPP() const override;
//???      bool operator==(const iter& i2) const override;
//???      bool operator!=(const iter& i2) const override;
    protected:
      const void* myPolyPtr() const override;
      const void* myTermPtr() const override;
    private: // data members
      const DistrMPolyInlFpPP& myPoly;
      const summand* mySummand;
    };

  };

  void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP&, long IndetIndex);
  void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, ConstRefPPMonoidElem x);


  // Friend functions - must also be declared outside the class (how very tedious).
  bool IsCompatible(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
  long NumTerms(const DistrMPolyInlFpPP& f);
  const DistrMPolyInlFpPP::InlineFpElem_t& LC(const DistrMPolyInlFpPP& f);
  ConstRefPPMonoidElem LPP(const DistrMPolyInlFpPP& f);
  int CmpLPP(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
  void MoveLMToFront(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g);
  void MoveLMToBack(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g);
  long log(const DistrMPolyInlFpPP& f, long IndetIndex);
  void DivLM(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
  bool IsZeroAddLCs(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g);
  void add(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
  void sub(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
  bool div(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
  void output(std::ostream& out, const DistrMPolyInlFpPP& f);
  bool IsZero(const DistrMPolyInlFpPP& f);
  //  bool IsOne(const DistrMPolyInlFpPP& f);
  //  bool IsMinusOne(const DistrMPolyInlFpPP& f);
  //  bool IsConstant(const DistrMPolyInlFpPP& f);
  //  bool IsIndet(long& index, const DistrMPolyInlFpPP& f);
  bool IsEqual(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g);
  void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP&, long IndetIndex);
  void deriv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& f, ConstRefPPMonoidElem x);

//   void WeylMul(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);
//   void WeylDiv(DistrMPolyInlFpPP& lhs, const DistrMPolyInlFpPP& g, const DistrMPolyInlFpPP& h);

  //---------------------------------------------------------------------------
  // Inline functions.

  inline std::size_t DistrMPolyInlFpPP::SummandSize(const ring&, const OrdvArith::reference& OA)
  {
    return sizeof(summand) + sizeof(OrdvArith::OrdvElem)*(OrdvWords(OA)-1);
  }


  inline int CmpLPP(const DistrMPolyInlFpPP& f, const DistrMPolyInlFpPP& g)
  {
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    return f.myOrdvArith->myCmp(f.mySummands->myOrdv, g.mySummands->myOrdv);
  }


  inline DistrMPolyInlFpPP::iter::iter(const DistrMPolyInlFpPP& f):
      myPoly(f),
      mySummand(f.mySummands)
  {}


  inline DistrMPolyInlFpPP::iter::iter(const DistrMPolyInlFpPP& f, const summand* TermPtr):
      myPoly(f),
      mySummand(TermPtr) // TermPtr *MUST* point to a term in f, or be nullptr.
  {}


  inline DistrMPolyInlFpPP::iter* DistrMPolyInlFpPP::iter::myClone() const
  {
    return new iter(*this);
  }


  inline bool DistrMPolyInlFpPP::iter::IamEnded() const
  {
    return mySummand == nullptr;
  }

  inline bool IsEnded(const DistrMPolyInlFpPP::iter& i)
  {
    return i.IamEnded();
  }


  inline void DistrMPolyInlFpPP::iter::myForward()
  {
    CoCoA_ASSERT(!IamEnded());
    mySummand = mySummand->myNext;
  }

  inline DistrMPolyInlFpPP::iter& DistrMPolyInlFpPP::iter::operator++()
  {
    if (IamEnded())  CoCoA_THROW_ERROR2(ERR::PolyIterEnded, "DistrMPolyInlFpPP::iter");
    myForward(); //???mySummand = mySummand->myNext;
    return *this;
  }


  inline DistrMPolyInlFpPP::iter DistrMPolyInlFpPP::iter::operator++(int)
  {
    if (IamEnded())  CoCoA_THROW_ERROR2(ERR::PolyIterEnded, "DistrMPolyInlFpPP::iter");
    iter ans(*this);
    myForward(); //???mySummand = mySummand->myNext;
    return ans;
  }


  inline RingElemAlias DistrMPolyInlFpPP::iter::myCoeff() const
  {
    CoCoA_ASSERT(!IamEnded());
    return RingElemAlias(myPoly.myCoeffRing, RingElemConstRawPtr(&(mySummand->myCoeff)));
  }

  inline RingElemAlias coeff(const DistrMPolyInlFpPP::iter& i)
  {
////    if (IsEnded(i))  CoCoA_THROW_ERROR(ERR::PolyIterEnded, "coeff(PolyIter)");
    if (IsEnded(i))  CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myCoeff();
  }


  inline ConstRefPPMonoidElem DistrMPolyInlFpPP::iter::myPP() const
  {
    CoCoA_ASSERT(!IamEnded());
    return ConstRefPPMonoidElem(myPoly.myPPM, PPMonoidElemConstRawPtr(mySummand->myOrdv));
  }

  inline ConstRefPPMonoidElem PP(const DistrMPolyInlFpPP::iter& i)
  {
////    if (IsEnded(i)) CoCoA_THROW_ERROR(ERR::PolyIterEnded, "PP(PolyIter)");
    if (IsEnded(i))  CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myPP();
  }


//   inline bool DistrMPolyInlFpPP::iter::operator==(const iter& i2) const
//   {
//     if (&myPoly != &(i2.myPoly)) CoCoA_THROW_ERROR(ERR::MixedPolyIters, "PolyIter == PolyIter");
//     return mySummand == i2.mySummand;
//   }

//   inline bool DistrMPolyInlFpPP::iter::operator!=(const iter& i2) const
//   {
//     if (&myPoly != &(i2.myPoly)) CoCoA_THROW_ERROR(ERR::MixedPolyIters, "PolyIter != PolyIter");
//     return mySummand != i2.mySummand;
//   }


  inline const void* DistrMPolyInlFpPP::iter::myPolyPtr() const
  {
    return &myPoly;
  }


  inline const void* DistrMPolyInlFpPP::iter::myTermPtr() const
  {
    return mySummand;
  }


  inline void swap(DistrMPolyInlFpPP& f, DistrMPolyInlFpPP& g)
  {
    DistrMPolyInlFpPP::ourSwap(f, g);
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/DistrMPolyInlFpPP.H,v 1.34 2024/07/31 07:52:46 abbott Exp $
// $Log: DistrMPolyInlFpPP.H,v $
// Revision 1.34  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.33  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.32  2021/11/03 16:47:35  abbott
// Summary: Commented out myFp field in NewSummandPtr (clang says it is unused)
//
// Revision 1.31  2021/10/30 19:03:07  abbott
// Summary: Added some more "override" (clang said they were mising; redmine 1625)
//
// Revision 1.30  2021/10/30 17:38:37  abbott
// Summary: Used keywords override & delete (redmine 1625 & 1627)
//
// Revision 1.29  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.28  2019/10/15 11:55:39  abbott
// Summary: Changed 0 into nullptr (where appropriate)
//
// Revision 1.27  2019/03/19 11:07:06  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.26  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.25  2018/05/17 15:28:08  bigatti
// -- added include SparsePolyIter
// -- removed doxygen initial comment
//
// Revision 1.24  2015/11/04 10:31:22  abbott
// Summary: Consequential changes (after revising SmallFpImpl)
//
// Revision 1.23  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.22  2013/03/25 17:26:31  abbott
// Changed return type of NumTerms (now long, was size_t).
// Cleaned impl of 2 fns.
// There is A LOT of cruft here -- major cleaning needed!!
//
// Revision 1.21  2012/10/24 11:40:56  abbott
// Changed return type of coeff/myCoeff.
//
// Revision 1.20  2012/10/16 10:28:18  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.19  2012/10/05 15:35:43  bigatti
// -- added myAddMonomial
//
// Revision 1.18  2012/09/26 12:25:45  abbott
// Updated to new SmallFpImpl interface.
//
// Revision 1.17  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.16  2011/09/15 11:55:49  abbott
// Added NewSummandPtr::myCoeffIsZero to achieve better encapsulation
// (and avoid getting an error from g++ :-)
//
// Revision 1.15  2011/09/13 15:47:58  abbott
// Moved defn of "grab" to corret place (outside class NewSummandPtr).
//
// Revision 1.14  2011/08/24 10:20:28  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.13  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.12  2009/10/28 16:48:15  bigatti
// -- commented out include <cstddef> (including MemPool.H)
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
// Revision 1.8  2008/04/10 14:53:18  bigatti
// -- added prototype for myPushBack/Front(RingElemConstRawPtr, PPMonoidElemConstRawPtr)
//
// Revision 1.7  2007/12/05 12:11:07  bigatti
// -- cleaning (mostly removing unused code)
//
// Revision 1.6  2007/12/05 10:57:47  bigatti
// -- removed some commented code (abstract definitions in SparsePolyRing)
//
// Revision 1.5  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/09/28 13:15:18  bigatti
// -- minor: modified some comments
//
// Revision 1.3  2007/09/28 13:11:46  bigatti
// -- removed "include "CoCoA/PPMonoid.H"
//
// Revision 1.2  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.18  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.17  2007/03/07 13:42:46  bigatti
// -- removed useless argument and other minor changes
//
// Revision 1.16  2007/01/11 14:04:18  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.15  2006/12/06 17:09:37  cocoa
// -- removed #include "config.H"
//
// Revision 1.14  2006/11/24 17:16:29  cocoa
// -- reorganized includes of header files
//
// Revision 1.13  2006/11/23 18:01:53  cocoa
// -- moved printing functions in unified implementation in SparsePolyRing
// -- simplified "output(f)" for debugging only
//
// Revision 1.12  2006/11/21 18:09:24  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.11  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.10  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.9  2006/10/16 23:19:00  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.8  2006/10/09 21:48:51  cocoa
// Commented out "pointless" template specialization of std::swap.
// Added overloading swap functions in CoCoA namespace.
//
// Revision 1.7  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.6  2006/10/06 10:00:18  cocoa
// Added independent declarations of functions declared as friends -- permits
// compilation with gcc-4.1 series.  Changes some friends into static member
// functions.
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
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.8  2006/05/12 16:10:59  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.7  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.6  2006/03/30 16:59:27  cocoa
// -- changed misleading name: InlineCoeffRing --> InlineCoeffImpl
// -- new: implementation for homomorphisms
// -- rearrangement of code to mimic RingDistrMPolyInlPP
//
// Revision 1.5  2006/03/20 17:27:42  cocoa
// -- changed in DistrMPolyInlFpPP: myMul, myDiv --> myMulByCoeff, myMulByPP, myDivByCoeff
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2006/03/07 10:06:12  cocoa
// -- fixed: PPMonoidElem LPP(f) now returns ConstRefPPMonoidElem
//
// Revision 1.2  2006/02/13 13:17:40  cocoa
// -- fixed: "const PPMonoidElem&" --> "ConstRefPPMonoidElem"
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.7  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
// Revision 1.6  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.5  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.4  2005/07/15 16:34:33  cocoa
// Added iterators for sparse polynomials.
// The code compiles (and the old tests still run).
// It'd Friday evening -- I'm going home before
// getting any ideas about making the iterator code run.
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
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
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
// Revision 1.9  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.8  2004/11/11 13:11:09  cocoa
// -- minor changes for doxygen
//
// Revision 1.7  2004/11/02 15:17:26  cocoa
// -- minor changes for doxygen
//
// Revision 1.6  2004/10/29 15:26:18  cocoa
// -- code fixed for compatibility with OrdvArith
//
// Revision 1.4  2004/07/20 15:37:08  cocoa
// Minor fix for some errors which slipped through the net...
//
// Revision 1.3  2004/07/20 09:20:47  cocoa
// -- myILineCoeffRing is now a reference!
//
// Revision 1.2  2004/07/16 10:11:34  cocoa
// -- now using the new class SmallFpImpl (or SmallFpLogImpl)
// -- updated with "my" coding convenctions
// -- NYI: LC and LCRaw
//
// Revision 1.1  2004/06/25 16:03:58  cocoa
// -- first import
//

#endif
