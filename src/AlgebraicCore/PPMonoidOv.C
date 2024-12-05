//   Copyright (c)  2001-2017,2021  John Abbott and Anna M. Bigatti
//   Author:  2005-2010  John Abbott

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


// Implementation of classes PPMonoidOvImpl

#include "CoCoA/PPMonoidOv.H"

#include "CoCoA/BigInt.H"
#include "CoCoA/DivMask.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OrdvArith.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"
#include "CoCoA/verbose.H"
//#include "CoCoA/VectorOps.H" // for debugging only


#include <algorithm>
using std::min;
using std::max;
#include <iostream>
using std::ostream;
#include<limits>
using std::numeric_limits;
#include <memory>
using std::unique_ptr;
#include <vector>
using std::vector;


namespace CoCoA
{

  class PPMonoidOvImpl: public PPMonoidBase
  {
    typedef PPMonoidElemRawPtr RawPtr;           // just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; // just to save typing
    typedef OrdvArith::OrdvElem OrdvElem;        // just to save typing

    static const unsigned long ourMaxExp;        // defined below, it is just numeric_limits<SmallExponent_t>::max()

  public:
    PPMonoidOvImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    ~PPMonoidOvImpl();
  public: // disable copy ctor and assignment
    PPMonoidOvImpl(const PPMonoidOvImpl& copy) = delete;
    PPMonoidOvImpl& operator=(const PPMonoidOvImpl& rhs) = delete;

  public:
    void contents() const; // FOR DEBUGGING ONLY

    const std::vector<PPMonoidElem>& myIndets() const override;                  ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

    // The functions below are operations on power products owned by PPMonoidOvImpl
    const PPMonoidElem& myOne() const override;
    PPMonoidElemRawPtr myNew() const override;                                   ///< ctor from nothing
    PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const override;      ///< ctor by assuming ownership
    PPMonoidElemRawPtr myNew(const std::vector<long>& expv) const override;      ///< ctor from exp vector
//NYI    PPMonoidElemRawPtr myNew(const std::vector<BigInt>& EXPV) const;///< ctor from exp vector
    void myDelete(RawPtr rawpp) const override;                                  ///< dtor, frees pp

    void mySwap(RawPtr rawpp1, RawPtr rawpp2) const override;                    ///< swap(pp1, pp2);
    void myAssignOne(RawPtr rawpp) const override;                               ///< pp = 1
    void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const override;              ///< p = pp1
    void myAssign(RawPtr rawpp, const std::vector<long>& expv) const override;   ///< pp = expv (assign from exp vector)

    void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = pp1*pp2
    void myMulIndetPower(RawPtr rawpp, long indet, long exp) const override;           ///< pp *= indet^exp, assumes exp >= 0
    void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = pp1/pp2
    void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override; ///< pp = pp1/gcd(pp1,pp2)
    void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = gcd(pp1,pp2)
    void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = lcm(pp1,pp2)
    void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const override;                  ///< pp = radical(pp1)
    void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const override;  ///< pp = pp1^exp, assumes exp >=  0
    void myPowerOverflowCheck(ConstRawPtr rawpp1, long exp) const override;            ///< throw if pp1^exp would overflow, assumes exp >= 0

    bool myIsOne(ConstRawPtr rawpp) const override;                                ///< is pp = 1?
    bool myIsIndet(long& index, ConstRawPtr rawpp) const override;                 ///< true iff pp is an indet
    bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;       ///< are pp1 & pp2 coprime?
    bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;         ///< is pp1 equal to pp2?
    bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;     ///< does pp2 divide pp1?
    bool myIsSqFree(ConstRawPtr rawpp) const override;                             ///< is pp equal to its radical?

    int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;              ///< -1,0,1 as pp1 < = > pp2
    long myStdDeg(ConstRawPtr rawpp) const override;                               ///< standard degree of pp
    void myWDeg(degree& d, ConstRawPtr rawpp) const override;                      ///< d = grading(pp)
    int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;          ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long) const override; ///< as myCmpWDeg wrt the first weights
    long myExponent(ConstRawPtr rawpp, long indet) const override;                 ///< exponent of indet in pp
    void myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const override; ///< EXP = exponent of indet in pp
    void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const override;   ///< expv[i] = exponent(pp,i)
    void myBigExponents(std::vector<BigInt>& v, ConstRawPtr rawpp) const override; ///< get exponents, SHOULD BE myExponents ???
    void myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const override;       ///< v[i] = true if i-th indet has exponent != 0
    void myOutputSelf(std::ostream& out) const override;                           ///< out << PPM
    // INHERITED DEFINITION of virtual  void myOutput(std::ostream& out, ConstRawPtr rawpp) const;
    void myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const;        ///< print pp in debugging format???


  private: // auxiliary functions
    OrdvElem* myOrdv(RawPtr) const;
    const OrdvElem* myOrdv(ConstRawPtr) const;

    void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const override; ///< used by PPWithMask
    void myComputeExpv(std::vector<long>& expv, RawPtr rawpp) const;
    bool myCheckExponents(const std::vector<long>& expv) const;
//???    void mySetExpv(RawPtr, const std::vector<long>& expv) const;

  private: // data members
    ///@name Class members
    //@{
    OrdvArith::reference myOrdvArith;  //??? should be const
    const long myEntrySize;
    mutable vector<long> myExpv1;  // buffer space
    mutable vector<long> myExpv2;  // buffer space
    mutable MemPool myMemMgr;     // IMPORTANT: this must come *before* myIndetVector and myOnePtr.
//???    std::vector<long> myDelta;
    vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
    unique_ptr<PPMonoidElem> myOnePtr;
    //@}
  };

  // static constant value
  const unsigned long PPMonoidOvImpl::ourMaxExp = numeric_limits<SmallExponent_t>::max();

  // File local inline functions

  inline PPMonoidOvImpl::OrdvElem* PPMonoidOvImpl::myOrdv(RawPtr rawpp) const
  {
    return static_cast<OrdvElem*>(rawpp.myRawPtr());

  }

  inline const PPMonoidOvImpl::OrdvElem* PPMonoidOvImpl::myOrdv(ConstRawPtr rawpp) const
  {
    return static_cast<const OrdvElem*>(rawpp.myRawPtr());
  }


  bool PPMonoidOvImpl::myCheckExponents(const std::vector<long>& expv) const
  {
    // Check len(expv) == myNumIndets.
    // Check exps are non-neg and not too big.
    if (len(expv) != myNumIndets) return false;
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] < 0 || static_cast<unsigned long>(expv[i]) > numeric_limits<SmallExponent_t>::max()) return false;
    return true;
  }


  //----   Constructors & destructor   ----//

  PPMonoidOvImpl::PPMonoidOvImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord):
      PPMonoidBase(ord, IndetNames),
      myOrdvArith(NewOrdvArith(ord)),
      myEntrySize(sizeof(OrdvElem)*OrdvWords(myOrdvArith)),
      myExpv1(myNumIndets),
      myExpv2(myNumIndets),
      myMemMgr(myEntrySize, "PPMonoidOvImpl.myMemMgr"),
      myIndetVector()
  {
    VerboseLog VERBOSE("PPMonoidOvImpl(syms,ord)");
    VERBOSE(99) << "-- called --" << std::endl;
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myOnePtr.reset(new PPMonoidElem(PPMonoid(this)));
    myIndetVector.reserve(myNumIndets);
    {
      // IMPORTANT: this block destroys pp *before* the call to myRefCountZero.
      PPMonoidElem pp(PPMonoid(this));
      vector<long> expv(myNumIndets);
      for (long i=0; i < myNumIndets; ++i)
      {
        expv[i] = 1;
        myAssign(raw(pp), expv);
        myIndetVector.push_back(pp);
        expv[i] = 0;
      }
    }
    myRefCountZero();
  }


  PPMonoidOvImpl::~PPMonoidOvImpl()
  {}


  // mySmartPtr is private
//   const OrdvArith::reference& OrdvA(const PPMonoid& PPM)
//   {
//     if (!IsPPMonoidOv(PPM))
//       CoCoA_THROW_ERROR("PPM must be PPMonoidOvImpl", OrdvA);
//     return dynamic_cast<const PPMonoidOvImpl*>(PPM.mySmartPtr()).myOrdvArith;
//   }


/////////////////////////////////////////////////////////////////////////////


  const std::vector<PPMonoidElem>& PPMonoidOvImpl::myIndets() const
  {
    return myIndetVector;
  }


  const PPMonoidElem& PPMonoidOvImpl::myOne() const
  {
    return *myOnePtr;
  }


  PPMonoidElemRawPtr PPMonoidOvImpl::myNew() const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssignOne(rawpp); // cannot throw
    return rawpp;
  }

  PPMonoidElemRawPtr PPMonoidOvImpl::myNew(PPMonoidElemConstRawPtr rawcopypp) const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssign(rawpp, rawcopypp); // cannot throw
    return rawpp;
  }


  PPMonoidElemRawPtr PPMonoidOvImpl::myNew(const std::vector<long>& expv) const
  {
    CoCoA_ASSERT(myCheckExponents(expv));
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssign(rawpp, expv); // cannot throw
    return rawpp;
  }


  void PPMonoidOvImpl::myAssignOne(RawPtr rawpp) const
  {
    myOrdvArith->myAssignZero(myOrdv(rawpp));
  }


  void PPMonoidOvImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    if (rawpp == rawpp1) return;
    myOrdvArith->myAssign(myOrdv(rawpp), myOrdv(rawpp1));
  }

  void PPMonoidOvImpl::myAssign(RawPtr rawpp, const vector<long>& expv) const
  {
    CoCoA_ASSERT(myCheckExponents(expv));
    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), expv);
  }


  void PPMonoidOvImpl::myDelete(RawPtr rawpp) const
  {
    myMemMgr.free(rawpp.myRawPtr());
  }


  void PPMonoidOvImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
  {
    if (rawpp1 == rawpp2) return;
    myOrdvArith->mySwap(myOrdv(rawpp1), myOrdv(rawpp2));
  }


  void PPMonoidOvImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    myOrdvArith->myMul(myOrdv(rawpp), myOrdv(rawpp1), myOrdv(rawpp2));
  }


  void PPMonoidOvImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    CoCoA_ASSERT("Exponent Overflow" && ourMaxExp - myExponent(rawpp, indet) >= static_cast<unsigned long>(exp));
    myOrdvArith->myMulIndetPower(myOrdv(rawpp), indet, exp);
  }


  void PPMonoidOvImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    myOrdvArith->myDiv(myOrdv(rawpp), myOrdv(rawpp1), myOrdv(rawpp2));
  }


  void PPMonoidOvImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
    vector<long> myExpv2(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp1));
    myOrdvArith->myComputeExpv(myExpv2, myOrdv(rawpp2));

    for (long i = 0; i < myNumIndets; ++i)
    {
      if (myExpv1[i] > myExpv2[i])
        myExpv1[i] -= myExpv2[i];
      else
        myExpv1[i] = 0;
    }
    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), myExpv1);
  }


  void PPMonoidOvImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
    vector<long> myExpv2(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp1));
    myOrdvArith->myComputeExpv(myExpv2, myOrdv(rawpp2));

    for (long i = 0; i < myNumIndets; ++i)
      myExpv1[i] = min(myExpv1[i], myExpv2[i]);

    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), myExpv1);
  }


  void PPMonoidOvImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
    vector<long> myExpv2(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp1));
    myOrdvArith->myComputeExpv(myExpv2, myOrdv(rawpp2));

    for (long i = 0; i < myNumIndets; ++i)
      myExpv1[i] = max(myExpv1[i], myExpv2[i]);

    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), myExpv1);
  }


  void PPMonoidOvImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp1));

    for (long i = 0; i < myNumIndets; ++i)
      myExpv1[i] = (myExpv1[i] > 0);

    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), myExpv1);
  }


  void PPMonoidOvImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(exp >= 0);
#ifdef CoCoA_DEBUG
    myPowerOverflowCheck(rawpp1, exp);
#endif
    myOrdvArith->myPower(myOrdv(rawpp), myOrdv(rawpp1), exp);
  }


  void PPMonoidOvImpl::myPowerOverflowCheck(ConstRawPtr rawpp, long exp) const
  {
    CoCoA_ASSERT(exp >= 0);
    myOrdvArith->myPowerOverflowCheck(myOrdv(rawpp), exp);
  }


  bool PPMonoidOvImpl::myIsOne(ConstRawPtr rawpp) const
  {
    return myOrdvArith->myIsZero(myOrdv(rawpp));
  }


  bool PPMonoidOvImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
  {
    return myOrdvArith->myIsIndet(index, myOrdv(rawpp));
  }


  bool PPMonoidOvImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
    vector<long> myExpv2(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp1));
    myOrdvArith->myComputeExpv(myExpv2, myOrdv(rawpp2));

    for (long i = 0; i < myNumIndets; ++i)
      if (myExpv1[i] != 0 && myExpv2[i] != 0) return false;

    return true;
  }


  bool PPMonoidOvImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdvArith->myCmp(myOrdv(rawpp1), myOrdv(rawpp2))==0;
  }


  bool PPMonoidOvImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
    vector<long> myExpv2(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp1));
    myOrdvArith->myComputeExpv(myExpv2, myOrdv(rawpp2));

    for (long i = 0; i < myNumIndets; ++i)
      if (myExpv1[i] < myExpv2[i]) return false;

    return true;
  }


  bool PPMonoidOvImpl::myIsSqFree(ConstRawPtr rawpp) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp));
    for (long i = 0; i < myNumIndets; ++i)
      if (myExpv1[i] > 1) return false;
    return true;
  }


  int PPMonoidOvImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdvArith->myCmp(myOrdv(rawpp1), myOrdv(rawpp2));
  }


// // should potentially skip the first few packed ordv entries???
// int PPMonoidOvImpl::myHomogCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
// {
//   return myOrdvArith->myCmp(myOrdv(rawpp1), myOrdv(rawpp2));
// }


  long PPMonoidOvImpl::myStdDeg(ConstRawPtr rawpp) const
  {
    return myOrdvArith->myStdDeg(myOrdv(rawpp));
  }


  void PPMonoidOvImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
  {
    myOrdvArith->myWDeg(d, myOrdv(rawpp));
  }


  int PPMonoidOvImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdvArith->myCmpWDeg(myOrdv(rawpp1), myOrdv(rawpp2));
  }


  int PPMonoidOvImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long i) const
  {
    return myOrdvArith->myCmpWDegPartial(myOrdv(rawpp1), myOrdv(rawpp2), i);
  }


  long PPMonoidOvImpl::myExponent(ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    return myOrdvArith->myExponent(myOrdv(rawpp), indet);
  }

  void PPMonoidOvImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    EXP = myExponent(rawpp, indet);
  }


  void PPMonoidOvImpl::myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    myOrdvArith->myComputeExpv(expv, myOrdv(rawpp));
  }


  void PPMonoidOvImpl::myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    std::vector<long> v(myNumIndets);
    myOrdvArith->myComputeExpv(v, myOrdv(rawpp));
    for (long i=0; i < myNumIndets; ++i)  expv[i] = v[i];
  }


  void PPMonoidOvImpl::myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    // SLUG SLUG SLUG
    vector<long> expv(myNumIndets);
    myOrdvArith->myComputeExpv(expv, myOrdv(rawpp));
    for (int i=0; i < myNumIndets; ++i)
      if (expv[i] != 0) v[i] = true;
  }


  void PPMonoidOvImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
  {
#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp));
    vector<SmallExponent_t> expv(myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
      expv[i] = static_cast<SmallExponent_t>(myExpv1[i]); // no problem as exponent must be non-neg.
    DivMaskImpl->myAssignFromExpv(dm, &expv[0], myNumIndets);
  }


  void PPMonoidOvImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PPMonoidOv(" << myNumIndets << ", " << myOrd << ")";
  }


  void PPMonoidOvImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const
  {
    if (!out) return;  // short-cut for bad ostreams

#ifdef CoCoA_THREADSAFE_HACK
    vector<long> myExpv1(myNumIndets);
#endif
    myOrdvArith->myComputeExpv(myExpv1, myOrdv(rawpp));
    out << "DEBUG PP: myNumIndets=" << myNumIndets << ", exps=[";
    for (long i=0; i < myNumIndets; ++i)
      out << myExpv1[i] << " ";
    out << "]" << std::endl;
  }


  PPMonoid NewPPMonoidOv(const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
    // Sanity check on the indet names given.
    const long nvars = NumIndets(ord);

    if (len(IndetNames) != nvars)
      CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPPMonoidOv(IndetNames,ord)");
    if (!AreDistinct(IndetNames))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidOv(IndetNames,ord)");
    if (!AreArityConsistent(IndetNames))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidOv(IndetNames,ord)");

    return PPMonoid(new PPMonoidOvImpl(IndetNames, ord));
  }


  PPMonoid NewPPMonoidOv(const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  { return NewPPMonoidOv(IndetNames, OrdCtor(len(IndetNames))); }


  bool IsPPMonoidOv(const PPMonoid& PPM)
  {
    return dynamic_cast<const PPMonoidOvImpl*>(PPM.operator->()) != nullptr;
  }

} // end of namespace CoCoA
