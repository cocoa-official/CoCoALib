#ifndef CoCoA_DistrMPolyClean_H
#define CoCoA_DistrMPolyClean_H

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


// Header file for the class DistrMPolyClean

#include "CoCoA/assert.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/MemPool.H"

// #include <cstddef> --- already included in MemPool.H
// using std::size_t;
#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;


namespace CoCoA
{

  class DistrMPolyClean
  {
  private:
    // UGLY HACK TO CIRCUMVENT MSVC 2010 bug...
#if defined(_WINDOWS)
  public:
#endif

//???    friend class PolyIter;

/*-- struct summand -------------------------------------------------*/
/**
   \brief summand of a DistrMPolyClean

   Very similar to the DistrMPolyInlPP summand, but the coefficient
   and the power-product are stored "cleanly" as a RingElem and a
   PPMonoidElem.  Being the size of these fileds fixed, we may use a
   static member "ourMemMgr" for summands of any RingDistrMPoly.

 \*-----------------------------------------------------------------*/
    struct summand
    {
    public:
      summand(const ring& R, const PPMonoid& M);
      ~summand();

    public: // data members
      summand* myNext;  ///< pointer to next summand
      RingElem myCoeff; ///< coefficient of the summand (should be non-zero)
      PPMonoidElem myPP; ///< power-product of the summand

    public: // disable default ctor, copy ctor, assignment
      summand() = delete;
      summand(const summand&) = delete;
      summand& operator=(const summand&) = delete;
    };

    // Next class works a bit like auto_ptr<summand>...
    class NewSummandPtr
    {
    public:
      explicit NewSummandPtr(const DistrMPolyClean& f): myPtr(nullptr), myEnv(f) {}
      ~NewSummandPtr() { if (myPtr == nullptr) return;  myPtr->~summand(); myEnv.myMemMgr.free(myPtr); }
      const summand* get() const { return myPtr; }
      summand* operator->() { return myPtr; }
      const summand* operator->() const { return myPtr; }
      summand* relinquish()  { CoCoA_ASSERT(myPtr != nullptr); summand* ans = myPtr; myPtr = nullptr; return ans; }
      void myRenew() { CoCoA_ASSERT(myPtr == nullptr); myPtr = static_cast<summand*>(myEnv.myMemMgr.alloc()); new(myPtr) summand(myEnv.myCoeffRing, myEnv.myPPM); }
      bool myCoeffIsZero() { CoCoA_ASSERT(myPtr != nullptr); return IsZero(myPtr->myCoeff); }

    private: // data members
      summand* myPtr;
      const DistrMPolyClean& myEnv; // for access to myCoeffRing,  myPPM,  myMemMgr.
    };

    static void ourDeleteSummands(summand* ptr/*, const ring& R, const PPMonoid& M*/, MemPool& MemMgr);

  private: // data members (of DistrMPolyClean)
    const ring myCoeffRing;   ///< the ring to which all the coeffs belong
    const PPMonoid myPPM;     ///< the PP monoid to which all the PPs belong
    MemPool& myMemMgr;        ///< the MemPool which manages the memory occupied by the summands
    summand** myEnd;          ///< pointer to facilitate implementing the one-past-the-last iterator
    summand* mySummands;      ///< pointer to the first summand

  public:
    DistrMPolyClean(const ring&, const PPMonoid&, MemPool&);
    DistrMPolyClean(const DistrMPolyClean& copy);
    ~DistrMPolyClean();
    DistrMPolyClean& operator=(const DistrMPolyClean& rhs);
    DistrMPolyClean& operator=(const MachineInt& rhs);
    DistrMPolyClean& operator=(const BigInt& rhs);
    DistrMPolyClean& operator=(const BigRat& rhs);

    friend bool IsCompatible(const DistrMPolyClean& f, const DistrMPolyClean& g);
    static long ourSummandSize(const ring&, const PPMonoid& PPM);
    static void ourSwap(DistrMPolyClean&, DistrMPolyClean&);
    summand* myCopySummand(const summand* original) const;
    // functions depending solely on the structure
    void myAssignZero();

    bool myIsEqual(const summand* const lhs, const summand*const rhs) const;
    friend long NumTerms(const DistrMPolyClean& f);
    friend RingElemRawPtr RawLC(DistrMPolyClean& f);
    friend ConstRefRingElem LC(const DistrMPolyClean& f);
    friend ConstRefPPMonoidElem LPP(const DistrMPolyClean& f);
    friend int CmpLPP(const DistrMPolyClean& f, const DistrMPolyClean& g);
    friend void MoveLMToFront(DistrMPolyClean& f, DistrMPolyClean& g);
    friend void MoveLMToBack(DistrMPolyClean& f, DistrMPolyClean& g);
    void myDeleteLM();
//     friend void wdeg(degree& d, const DistrMPolyClean& f);
//     friend int CmpWDeg(const DistrMPolyClean& f, const DistrMPolyClean& g);
    void myAddMulLM(const DistrMPolyClean& h, const DistrMPolyClean& g, bool SkipLMg);
//    void myWeylAddMul(const DistrMPolyClean& h, const DistrMPolyClean& g);
    void myReductionStep(const DistrMPolyClean& g);
    void myReductionStepGCD(const DistrMPolyClean& g, RingElem& fscale);
    void myAddClear(DistrMPolyClean& g);
    void myAppendClear(DistrMPolyClean& g);
    friend void DivLM(DistrMPolyClean& lhs, const DistrMPolyClean& f, const DistrMPolyClean& g);
    void myMulByCoeff(RingElemConstRawPtr rawc);  ///< EXCEPTION SAFE
    bool myDivByCoeff(RingElemConstRawPtr rawc);  ///< EXCEPTION SAFE
    void myMulByPP(PPMonoidElemConstRawPtr rawpp);
    //    void myWeylMul(PPMonoidElemConstRawPtr rawpp);
    void myPushFront(RingElemConstRawPtr rawc, const std::vector<long>& expv);
    void myPushBack(RingElemConstRawPtr rawc, const std::vector<long>& expv);
    void myPushFront(RingElemConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp);
    void myPushBack(RingElemConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp);
    void myPushFront(summand* t);
    void myPushBack(summand* t);
    void myAddMonomial(const DistrMPolyClean& g); ///< EXCEPTION SAFE
    void myRemoveSummand(summand** prev_link);
    void myInsertSummand(summand* s, summand** prev_link);
    friend bool IsZeroAddLCs(DistrMPolyClean& f, DistrMPolyClean& g);
    void myNegate();

    friend void add(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
    friend void sub(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
    friend bool div(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
    friend void output(std::ostream& out, const DistrMPolyClean& f); ///< for debugging only
    friend bool IsZero(const DistrMPolyClean& f);
    friend bool IsMonomial(const DistrMPolyClean& f);
    friend bool IsEqual(const DistrMPolyClean& f, const DistrMPolyClean& g);
//     friend void WeylMul(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
//     friend void WeylDiv(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
  private:
     void myAddMulSummand(const summand* s, const DistrMPolyClean& g, bool SkipLMg);
    //     void myWeylAddMulSummand(const summand* s, const DistrMPolyClean& g, bool SkipLMg);

  public:
    // iter implements a SparsePolyIter for DistrMPolyClean objects
    class iter: public SparsePolyIterBase
    {
    public:
      explicit iter(const DistrMPolyClean& f);
      iter(const DistrMPolyClean& f, const summand* TermPtr);
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
      const DistrMPolyClean& myPoly;
      const summand* mySummandPtr;
    };

  };


  // Friend functions - must also be declared outside the class (how very tedious).
  bool IsCompatible(const DistrMPolyClean& f, const DistrMPolyClean& g);
  long NumTerms(const DistrMPolyClean& f);
  ConstRefRingElem LC(const DistrMPolyClean& f);
  ConstRefPPMonoidElem LPP(const DistrMPolyClean& f);
  int CmpLPP(const DistrMPolyClean& f, const DistrMPolyClean& g);
  void MoveLMToFront(DistrMPolyClean& f, DistrMPolyClean& g);
  void MoveLMToBack(DistrMPolyClean& f, DistrMPolyClean& g);
  void DivLM(DistrMPolyClean& lhs, const DistrMPolyClean& f, const DistrMPolyClean& g);
  bool IsZeroAddLCs(DistrMPolyClean& f, DistrMPolyClean& g);
  void add(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
  void sub(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
  bool div(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
  void output(std::ostream& out, const DistrMPolyClean& f);
  bool IsZero(const DistrMPolyClean& f);
  bool IsEqual(const DistrMPolyClean& f, const DistrMPolyClean& g);

//   void WeylMul(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);
//   void WeylDiv(DistrMPolyClean& lhs, const DistrMPolyClean& g, const DistrMPolyClean& h);


  //---------------------------------------------------------------------------
  // Inline functions.

  inline DistrMPolyClean::summand::summand(const ring& R, const PPMonoid& M):
      myNext(nullptr),
      myCoeff(R),
      myPP(M)
  {}

  inline DistrMPolyClean::summand::~summand()
  {}


  inline long DistrMPolyClean::ourSummandSize(const ring&, const PPMonoid& /*PPM*/)
  {
    return sizeof(summand);
  }


  inline int CmpLPP(const DistrMPolyClean& f, const DistrMPolyClean& g)
  {
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    return (f.myPPM)->myCmp(raw(f.mySummands->myPP), raw(g.mySummands->myPP));
  }


  inline DistrMPolyClean::iter::iter(const DistrMPolyClean& f):
      myPoly(f),
      mySummandPtr(f.mySummands)
  {}


  inline DistrMPolyClean::iter::iter(const DistrMPolyClean& f, const summand* TermPtr):
      myPoly(f),
      mySummandPtr(TermPtr) // TermPtr *MUST* point to a term in f, or be nullptr.
  {}


  inline DistrMPolyClean::iter* DistrMPolyClean::iter::myClone() const
  {
    return new iter(*this);
  }


  inline bool DistrMPolyClean::iter::IamEnded() const
  {
    return mySummandPtr == nullptr;
  }

  inline bool IsEnded(const DistrMPolyClean::iter& i)
  {
    return i.IamEnded();
  }


  inline void DistrMPolyClean::iter::myForward()
  {
    CoCoA_ASSERT(!IamEnded());
    mySummandPtr = mySummandPtr->myNext;
  }

  inline DistrMPolyClean::iter& DistrMPolyClean::iter::operator++()
  {
////    if (IamEnded()) CoCoA_THROW_ERROR(ERR::PolyIterEnded, "DistrMPolyClean::iter::op++");
    if (IamEnded())  CoCoA_THROW_ERROR2(ERR::PolyIterEnded, "DistrMPolyClean::iter");
    myForward(); //???mySummand = mySummand->myNext;
    return *this;
  }


  inline DistrMPolyClean::iter DistrMPolyClean::iter::operator++(int)
  {
////    if (IamEnded()) CoCoA_THROW_ERROR(ERR::PolyIterEnded, "DistrMPolyClean::iter::op++");
    if (IamEnded())  CoCoA_THROW_ERROR2(ERR::PolyIterEnded, "DistrMPolyClean::iter");
    iter ans(*this);
    myForward(); //???mySummand = mySummand->myNext;
    return ans;
  }


  inline RingElemAlias DistrMPolyClean::iter::myCoeff() const
  {
    CoCoA_ASSERT(!IamEnded());
    return RingElemAlias(mySummandPtr->myCoeff);
  }

  inline RingElemAlias coeff(const DistrMPolyClean::iter& i)
  {
    if (IsEnded(i))  CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myCoeff();
  }


  inline ConstRefPPMonoidElem DistrMPolyClean::iter::myPP() const
  {
    CoCoA_ASSERT(!IamEnded());
    return ConstRefPPMonoidElem(mySummandPtr->myPP);
  }

  inline ConstRefPPMonoidElem PP(const DistrMPolyClean::iter& i)
  {
    if (IsEnded(i))  CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myPP();
  }


//   inline bool DistrMPolyClean::iter::operator==(const iter& i2) const
//   {
//     if (&myPoly != &(i2.myPoly)) CoCoA_THROW_ERROR(ERR::MixedPolyIters, "PolyIter == PolyIter");
//     return mySummandPtr == i2.mySummandPtr;
//   }

//   inline bool DistrMPolyClean::iter::operator!=(const iter& i2) const
//   {
//     if (&myPoly != &(i2.myPoly)) CoCoA_THROW_ERROR(ERR::MixedPolyIters, "PolyIter != PolyIter");
//     return mySummandPtr != i2.mySummandPtr;
//   }


  inline const void* DistrMPolyClean::iter::myPolyPtr() const
  {
    return &myPoly;
  }


  inline const void* DistrMPolyClean::iter::myTermPtr() const
  {
    return mySummandPtr;
  }

  inline void swap(DistrMPolyClean& f, DistrMPolyClean& g)
  {
    DistrMPolyClean::ourSwap(f, g);
  }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/DistrMPolyClean.H,v 1.22 2024/07/31 07:52:46 abbott Exp $
// $Log: DistrMPolyClean.H,v $
// Revision 1.22  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.21  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.20  2021/10/30 19:03:07  abbott
// Summary: Added some more "override" (clang said they were mising; redmine 1625)
//
// Revision 1.19  2021/10/30 17:38:37  abbott
// Summary: Used keywords override & delete (redmine 1625 & 1627)
//
// Revision 1.18  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.17  2019/03/19 11:07:06  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.16  2019/03/18 17:56:49  abbott
// Summary: Changed 0 into nullptr where appropriate
//
// Revision 1.15  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.14  2018/05/17 15:28:08  bigatti
// -- added include SparsePolyIter
// -- removed doxygen initial comment
//
// Revision 1.13  2015/11/04 10:40:19  abbott
// Summary: Added NewSummandPtr (so that summands live in a MemPool)
//
// Revision 1.12  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.11  2014/04/30 16:01:30  abbott
// Summary: Replaced size_t by long; removed some cruft
// Author: JAA
//
// Revision 1.10  2012/10/24 11:40:31  abbott
// Changed return type of coeff/myCoeff.
//
// Revision 1.9  2012/10/16 10:27:34  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.8  2012/10/11 14:27:59  abbott
// Removed "semantically risky" function RefLC/RawLC from DistrMPoly*.
// Reimplemented myRecvTwinFloat in RingDistrMPoly* more cleanly (but
// new impl does make a wasteful copy of the coeff).
//
// Revision 1.7  2012/10/05 15:33:28  bigatti
// -- added myAddMonomial
//
// Revision 1.6  2012/04/03 10:13:13  abbott
// Added ugly hack to circumvent Microsoft MSVC 2010 bug.
//
// Revision 1.5  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.4  2011/08/24 10:20:28  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.3  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.2  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.1  2010/10/08 08:05:28  bigatti
// -- renamed (Ring)DistrMPoly --> (Ring)DistrMPolyClean
//
// Revision 1.11  2009/10/28 16:48:00  bigatti
// -- commented out include <cstddef> (including MemPool.H)
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
// Revision 1.7  2008/04/10 14:53:18  bigatti
// -- added prototype for myPushBack/Front(RingElemConstRawPtr, PPMonoidElemConstRawPtr)
//
// Revision 1.6  2007/12/05 10:57:47  bigatti
// -- removed some commented code (abstract definitions in SparsePolyRing)
//
// Revision 1.5  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/10/05 14:32:50  bigatti
// -- added comment about #include "CoCoA/SparsePolyRing.H"
//
// Revision 1.3  2007/09/27 13:18:54  bigatti
// -- commented out undefined functions (content, gcd, IsOne, IsMinusOne, ..)
// -- added some comments
//
// Revision 1.2  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.18  2007/03/07 13:42:46  bigatti
// -- removed useless argument and other minor changes
//
// Revision 1.17  2007/01/11 14:04:18  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.16  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPolyClean(..)
//
// Revision 1.15  2006/12/06 17:09:37  cocoa
// -- removed #include "config.H"
//
// Revision 1.14  2006/11/23 18:01:53  cocoa
// -- moved printing functions in unified implementation in SparsePolyRing
// -- simplified "output(f)" for debugging only
//
// Revision 1.13  2006/11/21 18:09:24  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPolyClean(..) and RingDistrMPolyClean(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.12  2006/11/14 17:17:13  cocoa
// -- fixed coding convention "our"
//
// Revision 1.11  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPolyClean
//
// Revision 1.10  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.9  2006/10/16 23:19:00  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyCleanInlPP.
//
// Revision 1.8  2006/10/12 20:28:26  cocoa
// Added an assertion.
//
// Revision 1.7  2006/10/09 21:48:51  cocoa
// Commented out "pointless" template specialization of std::swap.
// Added overloading swap functions in CoCoA namespace.
//
// Revision 1.6  2006/10/06 14:04:16  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.5  2006/10/06 10:00:18  cocoa
// Added independent declarations of functions declared as friends -- permits
// compilation with gcc-4.1 series.  Changes some friends into static member
// functions.
//
// Revision 1.4  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.3  2006/07/20 17:06:08  cocoa
// -- moved myStdDeg into SparsePolyRing
//
// Revision 1.2  2006/06/08 16:45:28  cocoa
// -- RingDistrMPolyClean*.H  have been "moved" into RingDistrMPolyClean*.C
// -- some coding conventions fixed in DistrMPolyClean*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.10  2006/05/12 16:10:59  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.9  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPolyClean* and DistrMPolyClean* have been disabled
//
// Revision 1.8  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.7  2006/03/17 18:08:51  cocoa
// -- changed: mul --> myMulByPP
//
// Revision 1.6  2006/03/16 17:52:16  cocoa
// -- changed: mul, div --> myMulByCoeff, myDivByCoeff
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
// Revision 1.2  2006/02/13 13:17:41  cocoa
// -- fixed: "const PPMonoidElem&" --> "ConstRefPPMonoidElem"
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
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
// Main additions are to SparsePolyRing, DistrMPolyClean*.
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
// Revision 1.8  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.7  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
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

#endif
