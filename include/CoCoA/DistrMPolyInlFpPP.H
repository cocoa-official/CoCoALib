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
    if (IsEnded(i))  CoCoA_THROW_ERROR1(ERR::PolyIterEnded);
    return i.myPP();
  }


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

#endif
