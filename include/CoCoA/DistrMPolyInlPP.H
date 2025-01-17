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

#endif
