#ifndef CoCoA_DistrMPolyInlPP_H
#define CoCoA_DistrMPolyInlPP_H

//   Copyright (c)  2005  John Abbott,  Anna M. Bigatti

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


// Header file for the class DistrMPolyInlPP


#include "CoCoA/assert.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/OrdvArith.H"
#include "CoCoA/MemPool.H"

// #include <cstddef> --- already included in MemPool.H
// using std::size_t;
#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;


namespace CoCoA
{

  class DistrMPolyInlPP
  {
  private:
    // Struct used to represent the summand in a DistrMPolyInlPP.
    struct summand
    {
    public: // data members -- deliberately public (to avoid many friend declarations)
      summand* myNext;
      RingElemRawPtr myCoeff;
      OrdvArith::OrdvElem myOrdv[1]; // myOrdv ABSOLUTELY MUST BE THE LAST DATA MEMBER.
                                     // The compiler will allocate this at the end,
                                     // so I can address elements with index greater than
                                     // 1 with impunity if I know the space exists.

    public: // disable default ctor, copy ctor, assignment
      summand() = delete;
      summand(const summand&) = delete;
      summand& operator=(const summand&) = delete;
    };
    summand* NewSummand(const DistrMPolyInlPP& f);
    static void ourDeleteSummands(summand* ptr, const ring& R, MemPool& MemMgr);

    // This class behaves rather like you'd want std::auto_ptr<summand> to behave.
    class AutoPtrSummand
    {
    public:
      explicit AutoPtrSummand(const DistrMPolyInlPP& f);
//???      AutoPtrSummand(const AutoPtrSummand& copy) noexcept;
      ~AutoPtrSummand() noexcept;
//???      AutoPtrSummand& operator=(const AutoPtrSummand& rhs) noexcept;
      summand& operator*() const noexcept;
      summand* operator->() noexcept;
//???      const summand* operator->() const noexcept;
      summand* get() const noexcept;
//???      summand* const get() const noexcept;
      summand* release() noexcept;
//???      void reset(summand* ptr = nullptr) noexcept;
      void myRenew();
    private: // data members
      summand* myPtr;
      MemPool& myMemMgr;
      const ring myR;
      OrdvArith::reference myOrdvArith;
    };


  private: // data members (of DistrMPolyInlPP)
    const ring myCoeffRing;   ///< the ring to which all the coeffs belong
    const PPMonoid myPPM;     ///< the PP monoid to which all the PPs belong
    OrdvArith::reference myOrdvArith;
    MemPool& mySummandMemory; ///< the MemPool which manages the memory occupied by the summands
    summand* mySummands;      ///< pointer to the first summand (or nullptr if poly is zero)
    summand* myLastSummand;   ///< pointer to the first summand (or nullptr if poly is zero)

    friend const ring& CoeffRing(const DistrMPolyInlPP& f);
    friend const PPMonoid& PPM(const DistrMPolyInlPP& f);
//??????    RingElemRawPtr myDenom;  // used only if the coeff ring is a field of fractions.

  public:
    DistrMPolyInlPP(const ring&, const PPMonoid&, const OrdvArith::reference&, MemPool&);
    DistrMPolyInlPP(const DistrMPolyInlPP& copy);
    ~DistrMPolyInlPP();
    DistrMPolyInlPP& operator=(const DistrMPolyInlPP& rhs);
    DistrMPolyInlPP& operator=(const MachineInt& rhs);
    DistrMPolyInlPP& operator=(const BigInt& rhs);
    DistrMPolyInlPP& operator=(const BigRat& rhs);

    friend bool IsCompatible(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
//???    static void CheckCompatible(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g, const ErrorContext& ErrCtx);
    static std::size_t SummandSize(const ring&, const OrdvArith::reference& OA);
    static void ourSwap(DistrMPolyInlPP&, DistrMPolyInlPP&);
    summand* myCopySummand(const summand* original) const;
    // functions depending solely on the structure
    void myAssignZero();

    bool myIsEqual(const summand* const lhs, const summand*const rhs) const;
    friend long NumTerms(const DistrMPolyInlPP& f);
    friend RingElemRawPtr RawLC(DistrMPolyInlPP& f);
    friend RingElemAlias LC(const DistrMPolyInlPP& f);
    //    friend RingElemRawPtr& RawLC(DistrMPolyInlPP& f);
    //    friend const RingElemRawPtr& RawLC(const DistrMPolyInlPP& f);
    //    friend RingElem content(const DistrMPolyInlPP& f);
    friend ConstRefPPMonoidElem LPP(const DistrMPolyInlPP& f);
    friend int CmpLPP(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
    friend void MoveLMToFront(DistrMPolyInlPP& f, DistrMPolyInlPP& g);
    friend void MoveLMToBack(DistrMPolyInlPP& f, DistrMPolyInlPP& g);
    void myDeleteLM();
//     friend void wdeg(degree& d, const DistrMPolyInlPP& f);
//     friend int CmpWDeg(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
    void myAddMulLM(const DistrMPolyInlPP& h, const DistrMPolyInlPP& g, bool SkipLMg);
//    void myWeylAddMul(const DistrMPolyInlPP& h, const DistrMPolyInlPP& g);
    void myReductionStep(const DistrMPolyInlPP& g);
    void myReductionStepGCD(const DistrMPolyInlPP& g, RingElem& fscale);
    void myAddClear(DistrMPolyInlPP& g);
    void myAppendClear(DistrMPolyInlPP& g) noexcept;
    friend void DivLM(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
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
    void myAddMonomial(const DistrMPolyInlPP& g); ///< EXCEPTION SAFE
  private: // these two are impl details
    void myRemoveSummand(summand** prev_link);
    void myInsertSummand(summand* s, summand** prev_link);
  public:
    friend bool IsZeroAddLCs(DistrMPolyInlPP& f, DistrMPolyInlPP& g);
    void myNegate();    ///< MIGHT THROW????

    friend void add(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
    friend void sub(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
    friend bool div(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
    friend void output(std::ostream& out, const DistrMPolyInlPP& f);
    friend bool IsZero(const DistrMPolyInlPP& f);
//     friend bool IsOne(const DistrMPolyInlPP& f);
//     friend bool IsMinusOne(const DistrMPolyInlPP& f);
//     friend bool IsConstant(const DistrMPolyInlPP& f);
//     friend bool IsIndet(std::size_t& index, const DistrMPolyInlPP& f);
    friend bool IsMonomial(const DistrMPolyInlPP& f);
    friend bool IsEqual(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
//     friend void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP&, std::size_t IndetIndex);
//     friend void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, ConstRefPPMonoidElem x);

//     friend void WeylMul(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
//     friend void WeylDiv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
  private:
     void myAddMulSummand(const summand* s, const DistrMPolyInlPP& g, bool SkipLMg);
    //     void myWeylAddMulSummand(const summand* s, const DistrMPolyInlPP& g, bool SkipLMg);

  public:
    // iter implements a SparsePolyIter for DistrMPolyInlPP objects
    class iter: public SparsePolyIterBase
    {
    public:
      explicit iter(const DistrMPolyInlPP& f);
      iter(const DistrMPolyInlPP& f, const summand* TermPtr);
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
      const DistrMPolyInlPP& myPoly;
      const summand* mySummand;
    };

  };


//   void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP&, std::size_t IndetIndex);
//   void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, ConstRefPPMonoidElem x);


  // Friend functions - must also be declared outside the class (how very tedious).
  const ring& CoeffRing(const DistrMPolyInlPP& f);
  const PPMonoid& PPM(const DistrMPolyInlPP& f);

  bool IsCompatible(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
  long NumTerms(const DistrMPolyInlPP& f);
  RingElemAlias LC(const DistrMPolyInlPP& f);
  ConstRefPPMonoidElem LPP(const DistrMPolyInlPP& f);
  int CmpLPP(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
  void MoveLMToFront(DistrMPolyInlPP& f, DistrMPolyInlPP& g);
  void MoveLMToBack(DistrMPolyInlPP& f, DistrMPolyInlPP& g);
  long log(const DistrMPolyInlPP& f, std::size_t IndetIndex);
  void DivLM(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
  bool IsZeroAddLCs(DistrMPolyInlPP& f, DistrMPolyInlPP& g);
  void add(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
  void sub(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
  bool div(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
  void output(std::ostream& out, const DistrMPolyInlPP& f);
  bool IsZero(const DistrMPolyInlPP& f);
  //  bool IsOne(const DistrMPolyInlPP& f);
  //  bool IsMinusOne(const DistrMPolyInlPP& f);
  //  bool IsConstant(const DistrMPolyInlPP& f);
  //  bool IsIndet(std::size_t& index, const DistrMPolyInlPP& f);
  bool IsEqual(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g);
//   void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP&, std::size_t IndetIndex);
//   void deriv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& f, ConstRefPPMonoidElem x);

//   void WeylMul(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);
//   void WeylDiv(DistrMPolyInlPP& lhs, const DistrMPolyInlPP& g, const DistrMPolyInlPP& h);


  //---------------------------------------------------------------------------
  // Inline functions.

  inline std::size_t DistrMPolyInlPP::SummandSize(const ring&, const OrdvArith::reference& OA)
  {
    return sizeof(summand) + sizeof(OrdvArith::OrdvElem)*(OrdvWords(OA)-1);
  }


  inline const ring& CoeffRing(const DistrMPolyInlPP& f)
  {
    return f.myCoeffRing;
  }


  inline const PPMonoid& PPM(const DistrMPolyInlPP& f)
  {
    return f.myPPM;
  }


//   // Just syntactic sugar
//   inline void DistrMPolyInlPP::myPushFront(RingElemConstRawPtr rawc, const std::vector<SmallExponent_t>& expv)
//   {
//     myPushFront(c, &expv[0]);
//   }

//   // Just syntactic sugar
//   inline void DistrMPolyInlPP::myPushBack(RingElemConstRawPtr rawc, const std::vector<SmallExponent_t>& expv)
//   {
//     myPushBack(c, &expv[0]);
//   }


  inline int CmpLPP(const DistrMPolyInlPP& f, const DistrMPolyInlPP& g)
  {
    CoCoA_ASSERT(!IsZero(f) && !IsZero(g));
    return f.myOrdvArith->myCmp(f.mySummands->myOrdv, g.mySummands->myOrdv);
  }


  inline DistrMPolyInlPP::iter::iter(const DistrMPolyInlPP& f):
      myPoly(f),
      mySummand(f.mySummands)
  {}


  inline DistrMPolyInlPP::iter::iter(const DistrMPolyInlPP& f, const summand* TermPtr):
      myPoly(f),
      mySummand(TermPtr) // TermPtr *MUST* point to a term in f, or be nullptr.
  {}


  inline DistrMPolyInlPP::iter* DistrMPolyInlPP::iter::myClone() const
  {
    return new iter(*this);
  }


  inline bool DistrMPolyInlPP::iter::IamEnded() const
  {
    return mySummand == nullptr;
  }

  inline bool IsEnded(const DistrMPolyInlPP::iter& i)
  {
    return i.IamEnded();
  }


  inline void DistrMPolyInlPP::iter::myForward()
  {
    CoCoA_ASSERT(!IamEnded());
    mySummand = mySummand->myNext;
  }

  inline DistrMPolyInlPP::iter& DistrMPolyInlPP::iter::operator++()
  {
    if (IamEnded())
      CoCoA_THROW_ERROR(ERR::PolyIterEnded, "DistrMPolyInlPP::iter::op++");
    myForward();
    return *this;
  }


  inline DistrMPolyInlPP::iter DistrMPolyInlPP::iter::operator++(int)
  {
    if (IamEnded())
      CoCoA_THROW_ERROR(ERR::PolyIterEnded, "DistrMPolyInlPP::iter::op++");
    iter ans(*this);
    myForward();
    return ans;
  }


  inline RingElemAlias DistrMPolyInlPP::iter::myCoeff() const
  {
    CoCoA_ASSERT(!IamEnded());
    return RingElemAlias(CoeffRing(myPoly), mySummand->myCoeff);
  }

  inline RingElemAlias coeff(const DistrMPolyInlPP::iter& i)
  {
    if (IsEnded(i))
      CoCoA_THROW_ERROR(ERR::PolyIterEnded, "coeff(PolyIter)");
    return i.myCoeff();
  }


  inline ConstRefPPMonoidElem DistrMPolyInlPP::iter::myPP() const
  {
    CoCoA_ASSERT(!IamEnded());
    return ConstRefPPMonoidElem(PPM(myPoly), PPMonoidElemConstRawPtr(mySummand->myOrdv));
  }

  inline ConstRefPPMonoidElem PP(const DistrMPolyInlPP::iter& i)
  {
    if (IsEnded(i)) CoCoA_THROW_ERROR(ERR::PolyIterEnded, "PP(PolyIter)");
    return i.myPP();
  }


//   inline bool DistrMPolyInlPP::iter::operator==(const iter& i2) const
//   {
//     if (&myPoly != &(i2.myPoly)) CoCoA_THROW_ERROR(ERR::MixedPolyIters, "PolyIter == PolyIter");
//     return mySummand == i2.mySummand;
//   }

//   inline bool DistrMPolyInlPP::iter::operator!=(const iter& i2) const
//   {
//     if (&myPoly != &(i2.myPoly)) CoCoA_THROW_ERROR(ERR::MixedPolyIters, "PolyIter != PolyIter");
//     return mySummand != i2.mySummand;
//   }


  inline const void* DistrMPolyInlPP::iter::myPolyPtr() const
  {
    return &myPoly;
  }


  inline const void* DistrMPolyInlPP::iter::myTermPtr() const
  {
    return mySummand;
  }


  inline void swap(DistrMPolyInlPP& f, DistrMPolyInlPP& g)
  {
    DistrMPolyInlPP::ourSwap(f, g);
  }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/DistrMPolyInlPP.H,v 1.32 2024/03/16 10:15:11 abbott Exp $
// $Log: DistrMPolyInlPP.H,v $
// Revision 1.32  2024/03/16 10:15:11  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.31  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.30  2021/10/30 19:34:28  abbott
// Summary: Added missing noexcept
//
// Revision 1.29  2021/10/30 19:03:07  abbott
// Summary: Added some more "override" (clang said they were missing; redmine 1625)
//
// Revision 1.28  2021/10/30 17:38:37  abbott
// Summary: Used keywords override & delete (redmine 1625 & 1627)
//
// Revision 1.27  2021/02/10 19:32:00  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.26  2021/01/07 14:11:50  abbott
// Summary: Corrected copyright
//
// Revision 1.25  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.24  2019/03/19 11:07:06  abbott
// Summary: Replaced 0 by nullptr where appropriate
//
// Revision 1.23  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.22  2018/05/17 15:28:08  bigatti
// -- added include SparsePolyIter
// -- removed doxygen initial comment
//
// Revision 1.21  2016/05/18 12:26:50  abbott
// Summary: Replaced data member myEnd by myLastSummand (major change!)
//
// Revision 1.20  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.19  2013/06/27 16:47:12  abbott
// Changed commented out fn.
//
// Revision 1.18  2012/10/24 11:41:39  abbott
// Changed return type of coeff/myCoeff and LC.
//
// Revision 1.17  2012/10/16 10:28:40  abbott
// Replaced  RefRingElem  by  RingElem&
//
// Revision 1.16  2012/10/11 14:27:59  abbott
// Removed "semantically risky" function RefLC/RawLC from DistrMPoly*.
// Reimplemented myRecvTwinFloat in RingDistrMPoly* more cleanly (but
// new impl does make a wasteful copy of the coeff).
//
// Revision 1.15  2012/10/05 15:35:13  bigatti
// -- changed renew --> myRenew
// -- added myAddMonomial
//
// Revision 1.14  2011/11/09 13:38:31  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.13  2011/08/24 10:20:28  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.12  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.11  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.10  2009/10/28 16:48:28  bigatti
// -- commented out include <cstddef> (including MemPool.H)
//
// Revision 1.9  2009/10/02 13:47:07  bigatti
// -- myDivByCoeff now returns bool
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.8  2009/09/28 17:14:41  bigatti
// -- commented out unused functions (div, deriv, *Weyl*)
//
// Revision 1.7  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.6  2008/04/10 14:53:18  bigatti
// -- added prototype for myPushBack/Front(RingElemConstRawPtr, PPMonoidElemConstRawPtr)
//
// Revision 1.5  2007/12/05 10:57:47  bigatti
// -- removed some commented code (abstract definitions in SparsePolyRing)
//
// Revision 1.4  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/28 14:19:49  bigatti
// -- added "my" to some member functions
// -- added some comments
// -- commented out some unused functions (implemented in SparsePolyRing)
//
// Revision 1.2  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.17  2007/03/07 13:42:45  bigatti
// -- removed useless argument and other minor changes
//
// Revision 1.16  2007/01/11 14:04:18  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.15  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.14  2006/12/06 17:09:37  cocoa
// -- removed #include "config.H"
//
// Revision 1.13  2006/11/24 17:16:29  cocoa
// -- reorganized includes of header files
//
// Revision 1.12  2006/11/23 18:01:53  cocoa
// -- moved printing functions in unified implementation in SparsePolyRing
// -- simplified "output(f)" for debugging only
//
// Revision 1.11  2006/11/21 18:09:24  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.10  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.9  2006/10/16 23:18:59  cocoa
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
// Revision 1.11  2006/05/12 16:10:59  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.10  2006/04/27 12:19:26  cocoa
// -- coding conventions (myAddMulSummand)
//
// Revision 1.9  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.8  2006/03/21 17:55:08  cocoa
// -- changed: "my.." coding conventions
//
// Revision 1.7  2006/03/17 18:08:51  cocoa
// -- changed: mul --> myMulByPP
//
// Revision 1.6  2006/03/16 17:56:14  cocoa
// -- changed: mul, div --> myMulByCoeff, myDivByCoeff (+check IsOne(c))
// -- changed: myOutputSelf for homomorphism
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
// Revision 1.8  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
//
// Revision 1.7  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.6  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
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
// Revision 1.15  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.14  2004/11/11 13:35:48  cocoa
// -- minor changes for doxygen
//
// Revision 1.13  2004/11/08 14:03:32  cocoa
// -- small changes for doxygen support
//
// Revision 1.12  2004/11/02 15:17:26  cocoa
// -- minor changes for doxygen
//
// Revision 1.11  2004/10/29 15:26:17  cocoa
// -- code fixed for compatibility with OrdvArith
//
// Revision 1.9  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
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
// Revision 1.3  2003/10/09 12:48:17  cocoa
// New coding convention for rings.
//
// Revision 1.2  2003/10/01 10:35:32  cocoa
// - applied "my" coding convention to PPMonoid and PPOrdering
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.5  2003/09/17 15:17:16  abbott
// Added missing "reference" to third arg of deriv.
//
// Revision 1.4  2003/06/23 17:03:06  abbott
// Minor cleaning prior to public release.
// Name change, and consequences.
//
// Revision 1.3  2003/05/08 10:40:35  abbott
// Checking in prior to a probable major overhaul (this code is a mess).
// Numerous changes, mostly due to name changes for rings and PPMonoids.
// Added function deriv to compute derivatives.
//
// Revision 1.2  2002/11/18 17:53:26  bigatti
// - AddMulSummand is now private
//
// Revision 1.1  2002/11/15 15:50:31  abbott
// Initial revision
//
//

#endif
