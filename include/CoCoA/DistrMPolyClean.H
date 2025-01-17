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
   PPMonoidElem.  Being the size of these fields fixed, we may use a
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

#endif
