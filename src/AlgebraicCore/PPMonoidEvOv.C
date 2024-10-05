//   Copyright (c)  2005,2007,2010,2021  John Abbott and Anna M. Bigatti

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


// Implementation of class PPMonoidEvOvImpl

#include "CoCoA/PPMonoidEvOv.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DivMask.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OrdvArith.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"

#include <algorithm>
using std::min;
using std::max;
//using std::swap;
#include <cstring>
using std::memcpy;
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

  /*-- class PPMonoidEvOvImpl ---------------------------------------*/
  /**

  \brief Implementation of power product monoid for fast generic use

  PPMonoidEvOv implements a power product monoid for generic use as
  it stores exp vector and order vector.
  Compared with PPMonoidSafe, PPMonoidEvOv is:
  - slower for gcd/lcm because it has to update the order vector;
  - faster for comparisons expecially with matrix defined orderings.

  So this type is good for you if
  (1) you do not perform many gcd/lcm
  (2) you need efficiency in ordering test

  */
  /*-----------------------------------------------------------------*/

  class PPMonoidEvOvImpl: public PPMonoidBase
  {
    typedef PPMonoidElemRawPtr RawPtr;           // just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; // just to save typing
    typedef OrdvArith::OrdvElem OrdvElem;        // just to save typing

    static const unsigned long ourMaxExp;        // defined below; value is just numeric_limits<SmallExponent_t>::max()

  public:
    PPMonoidEvOvImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    ~PPMonoidEvOvImpl();
  public: // disable copy ctor and assignment
    PPMonoidEvOvImpl(const PPMonoidEvOvImpl& copy) = delete;
    PPMonoidEvOvImpl& operator=(const PPMonoidEvOvImpl& rhs) = delete;

  public:
    void contents() const; // FOR DEBUGGING ONLY

    const std::vector<PPMonoidElem>& myIndets() const override;                  ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

    // The functions below are operations on power products owned by PPMonoidEvOvImpl
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
    void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const override;                   ///< pp = radical(pp1)
    void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const override;   ///< pp = pp1^exp, assumes exp >= 0
    void myPowerOverflowCheck(ConstRawPtr rawpp1, long exp) const override;            ///< throw if pp1^exp would overflow, assumes exp >= 0

    bool myIsOne(ConstRawPtr rawpp) const override;                              ///< is pp = 1?
    bool myIsIndet(long& index, ConstRawPtr rawpp) const override;               ///< true iff pp is an indet
    bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;     ///< are pp1 & pp2 coprime?
    bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;       ///< is pp1 equal to pp2?
    bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< does pp2 divide pp1?
    bool myIsSqFree(ConstRawPtr rawpp) const override;                           ///< is pp equal to its radical?

    int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;               ///< -1,0,1 as pp1 < = > pp2
    long myStdDeg(ConstRawPtr rawpp) const override;                                ///< standard degree of pp
    void myWDeg(degree& d, ConstRawPtr rawpp) const override;                       ///< d = grading(pp)
    int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;           ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long) const override; ///< as myCmpWDeg wrt the first weights
    long myExponent(ConstRawPtr rawpp, long indet) const override;                  ///< exponent of indet in pp
    void myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const override;  ///< EXP = exponent of indet in pp
    void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const override;    ///< expv[i] = exponent(pp,i)
    void myBigExponents(std::vector<BigInt>& v, ConstRawPtr rawpp) const override;  ///< get exponents, SHOULD BE myExponents ???
    void myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const override;        ///< v[i] = true if i-th indet has exponent != 0
    void myOutputSelf(std::ostream& out) const override;                            ///< out << PPM
    // INHERITED DEFINITION of virtual  void myOutput(std::ostream& out, ConstRawPtr rawpp) const;
    void myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const;   ///< print pp in debugging format???


  private: // auxiliary functions
    SmallExponent_t* myExpv(RawPtr) const;
    const SmallExponent_t* myExpv(ConstRawPtr) const;
    OrdvElem* myOrdv(RawPtr) const;
    const OrdvElem* myOrdv(ConstRawPtr) const;

    void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const override; ///< used by PPWithMask
    void myComputeOrdv(RawPtr) const;
    bool myCheckExponents(const std::vector<long>& expv) const;

  private: // data members
    ///@name Class members
    //@{
    OrdvArith::reference myOrdvArith;  //??? should be const
    const long myOrdvSize;        ///< used only in myExpv
    const long myEntrySize;       ///< size in bytes
    mutable MemPool myMemMgr;     // IMPORTANT: this must come *before* myIndetVector and myOnePtr.
//???    vector<SmallExponent_t> myDelta;
    vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
    unique_ptr<PPMonoidElem> myOnePtr;
    //@}
  };


  // static variable
  const unsigned long PPMonoidEvOvImpl::ourMaxExp = numeric_limits<SmallExponent_t>::max();


  // File local inline functions

  inline SmallExponent_t* PPMonoidEvOvImpl::myExpv(RawPtr rawpp) const
  {
    return reinterpret_cast<SmallExponent_t*>(static_cast<char*>(rawpp.myRawPtr()) + myOrdvSize);
  }


  inline const SmallExponent_t* PPMonoidEvOvImpl::myExpv(ConstRawPtr rawpp) const
  {
    return reinterpret_cast<const SmallExponent_t*>(static_cast<const char*>(rawpp.myRawPtr()) + myOrdvSize);
  }


  inline PPMonoidEvOvImpl::OrdvElem* PPMonoidEvOvImpl::myOrdv(RawPtr rawpp) const
  {
    return static_cast<OrdvElem*>(rawpp.myRawPtr());

  }

  inline const PPMonoidEvOvImpl::OrdvElem* PPMonoidEvOvImpl::myOrdv(ConstRawPtr rawpp) const
  {
    return static_cast<const OrdvElem*>(rawpp.myRawPtr());
  }


  bool PPMonoidEvOvImpl::myCheckExponents(const std::vector<long>& expv) const
  {
    // Check expv.size == myNumIndets.
    // Check exps are non-neg and not too big.
    if (len(expv) != myNumIndets) return false;
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] < 0 || static_cast<unsigned long>(expv[i]) > numeric_limits<SmallExponent_t>::max()) return false;
    return true;
  }


  //----   Constructors & destructor   ----//

  PPMonoidEvOvImpl::PPMonoidEvOvImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord):
      PPMonoidBase(ord, IndetNames),
      myOrdvArith(NewOrdvArith(ord)),
      myOrdvSize(sizeof(OrdvElem)*OrdvWords(myOrdvArith)),
      myEntrySize(myOrdvSize + sizeof(SmallExponent_t)*myNumIndets),
      myMemMgr(myEntrySize, "PPMonoidEvOvImpl.myMemMgr"),
      myIndetVector()
  {
    // std::cout << "------PPMonoidEvOvImpl:NewOrdvArith-called" << std::endl;
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myOnePtr.reset(new PPMonoidElem(PPMonoid(this)));
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


  PPMonoidEvOvImpl::~PPMonoidEvOvImpl()
  {}

/////////////////////////////////////////////////////////////////////////////


  inline void PPMonoidEvOvImpl::myComputeOrdv(RawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);
    vector<long> ExpvCopy(myNumIndets);
    for (long i=0; i < myNumIndets; ++i)
      ExpvCopy[i] = IntegerCast<long>(expv[i]);
    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), ExpvCopy);
  }


  const std::vector<PPMonoidElem>& PPMonoidEvOvImpl::myIndets() const
  {
    return myIndetVector;
  }


  const PPMonoidElem& PPMonoidEvOvImpl::myOne() const
  {
    return *myOnePtr;
  }


  PPMonoidElemRawPtr PPMonoidEvOvImpl::myNew() const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssignOne(rawpp); // cannot throw
    return rawpp;
  }

  PPMonoidElemRawPtr PPMonoidEvOvImpl::myNew(PPMonoidElemConstRawPtr rawcopypp) const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssign(rawpp, rawcopypp); // cannot throw
    return rawpp;
  }


  PPMonoidElemRawPtr PPMonoidEvOvImpl::myNew(const std::vector<long>& expv) const
  {
    CoCoA_ASSERT(myCheckExponents(expv));
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssign(rawpp, expv); // cannot throw
    return rawpp;
  }


  void PPMonoidEvOvImpl::myAssignOne(RawPtr rawpp) const
  {
    SmallExponent_t* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = 0;
    myOrdvArith->myAssignZero(myOrdv(rawpp));
  }


  void PPMonoidEvOvImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    if (rawpp == rawpp1) return;
//     // This code assumes that myEntrySize is an exact multiple of sizeof(int).
//     int* const expv = static_cast<int*>(rawpp.myRawPtr());
//     const int* const expv1 = static_cast<const int*>(rawpp1.myRawPtr());
//     const long NumWords = myEntrySize/sizeof(int);
//     std::copy(expv1, expv1+NumWords, expv);  // does this work???

    memcpy(myOrdv(rawpp), myOrdv(rawpp1), myEntrySize);

// This would be a cleaner (but slower) way of achieving the same result...
//     SmallExponent_t* const exp = myExpv(rawpp);
//     const SmallExponent_t* const exp_src = myExpv(src);
//     for (long i=0 ; i<myNumIndets ; ++i )  exp[i] = exp_src[i];
//     myOrdvArith->myAssign(myOrdv(rawpp), myOrdv(src));
  }

  void PPMonoidEvOvImpl::myAssign(RawPtr rawpp, const vector<long>& expv) const
  {
    CoCoA_ASSERT(myCheckExponents(expv));

    SmallExponent_t* const expv2 = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
      expv2[i] = IntegerCast<SmallExponent_t>(expv[i]);

    myOrdvArith->myAssignFromExpv(myOrdv(rawpp), expv);
  }


  void PPMonoidEvOvImpl::myDelete(RawPtr rawpp) const
  {
    myMemMgr.free(rawpp.myRawPtr());
  }


  void PPMonoidEvOvImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
  {
    if (rawpp1 == rawpp2) return;
    // This code assumes that myEntrySize is an exact multiple of sizeof(int)
    int* v1 = static_cast<int*>(rawpp1.myRawPtr());
    int* v2 = static_cast<int*>(rawpp2.myRawPtr());
    const long NumWords = myEntrySize/sizeof(int);
    for (long i=0; i < NumWords; ++i)
      std::swap(v1[i], v2[i]);
  }


  void PPMonoidEvOvImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);
    for (long i=0; i < myNumIndets; ++i)
    {
      CoCoA_ASSERT("Exponent Overflow" && expv1[i] <= std::numeric_limits<SmallExponent_t>::max()-expv2[i]);
      expv[i] = expv1[i] + expv2[i];
    }
    myOrdvArith->myMul(myOrdv(rawpp), myOrdv(rawpp1), myOrdv(rawpp2));
  }


  void PPMonoidEvOvImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    SmallExponent_t* const expv = myExpv(rawpp);
    // If CoCoA_DEBUG active, check for exponent overflow...
    CoCoA_ASSERT("Exponent Overflow" && ourMaxExp - expv[indet] >= static_cast<unsigned long>(exp));
    expv[indet] += static_cast<SmallExponent_t>(exp);  // cast to keep M$ compiler quiet
    myOrdvArith->myMulIndetPower(myOrdv(rawpp), indet, exp);
  }


  void PPMonoidEvOvImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);
    for (long i=0; i < myNumIndets; ++i)
    {
      CoCoA_ASSERT("Exponent Underflow" && expv1[i] >= expv2[i]);
      expv[i] = expv1[i] - expv2[i];
    }
    myOrdvArith->myDiv(myOrdv(rawpp), myOrdv(rawpp1), myOrdv(rawpp2));
  }


  void PPMonoidEvOvImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
    {
      if (expv1[i] > expv2[i])
        expv[i] = expv1[i] - expv2[i];
      else
        expv[i] = 0;
    }
    myComputeOrdv(rawpp);
  }


  void PPMonoidEvOvImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = min(expv1[i], expv2[i]);

    myComputeOrdv(rawpp);
  }


  void PPMonoidEvOvImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = max(expv1[i], expv2[i]);

    myComputeOrdv(rawpp);
  }


  void PPMonoidEvOvImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = (expv1[i] > 0);

    myComputeOrdv(rawpp);
  }


  void PPMonoidEvOvImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long LongExp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(LongExp >= 0);
#ifdef CoCoA_DEBUG
    myPowerOverflowCheck(rawpp1, LongExp);
#endif
    if (static_cast<unsigned long>(LongExp) > ourMaxExp)
      CoCoA_THROW_ERROR(ERR::ExpTooBig, "PPMonoidEvOvImpl::myPowerSmallExp");
    const SmallExponent_t exp = static_cast<SmallExponent_t>(LongExp);

    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = exp * expv1[i];
    myOrdvArith->myPower(myOrdv(rawpp), myOrdv(rawpp1), exp);
  }


  void PPMonoidEvOvImpl::myPowerOverflowCheck(ConstRawPtr rawpp, long LongExp) const
  {
    if (LongExp == 0 || LongExp == 1) return;
    CoCoA_ASSERT(LongExp >= 0);
    const char* const FnName = "PPMonoidEvOvImpl::myPowerOverflowCheck";
    if (static_cast<unsigned long>(LongExp) > ourMaxExp)
      CoCoA_THROW_ERROR(ERR::ExpTooBig, FnName);
    const SmallExponent_t exp = static_cast<SmallExponent_t>(LongExp);
    const SmallExponent_t limit = ourMaxExp/exp;

    const SmallExponent_t* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (expv[i] > limit)
        CoCoA_THROW_ERROR(ERR::ExpTooBig, FnName);
    }
    // Check separately for overflow in ordv
    myOrdvArith->myPowerOverflowCheck(myOrdv(rawpp), exp);
  }


  bool PPMonoidEvOvImpl::myIsOne(ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv[i] != 0) return false;

    return true;
  }


  bool PPMonoidEvOvImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);
    long j = myNumIndets;
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (expv[i] == 0) continue;
      if (j != myNumIndets || expv[i] != 1) return false;
      j = i;
    }
    if (j == myNumIndets) return false;
    index = j;
    return true;
  }


  bool PPMonoidEvOvImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] != 0 && expv2[i] != 0) return false;

    return true;
  }


  bool PPMonoidEvOvImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdvArith->myCmp(myOrdv(rawpp1), myOrdv(rawpp2)) == 0;
  }


  bool PPMonoidEvOvImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] < expv2[i]) return false;

    return true;
  }


  bool PPMonoidEvOvImpl::myIsSqFree(ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv[i] > 1) return false;

    return true;
  }


  int PPMonoidEvOvImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdvArith->myCmp(myOrdv(rawpp1), myOrdv(rawpp2));
  }


// // should potentially skip the first few packed ordv entries???
// int PPMonoidEvOvImpl::myHomogCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
// {
//   return myOrdvArith->myCmp(myOrdv(rawpp1), myOrdv(rawpp2));
// }


  long PPMonoidEvOvImpl::myStdDeg(ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);
    long d=0;
    for (long i=0; i < myNumIndets; ++i)
      d += expv[i];
    return d;
  }


  void PPMonoidEvOvImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
  {
    myOrdvArith->myWDeg(d, myOrdv(rawpp));
  }


  int PPMonoidEvOvImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdvArith->myCmpWDeg(myOrdv(rawpp1), myOrdv(rawpp2));
  }


  int PPMonoidEvOvImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long i) const
  { 
    return myOrdvArith->myCmpWDegPartial(myOrdv(rawpp1), myOrdv(rawpp2), i);
  }


  long PPMonoidEvOvImpl::myExponent(ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    return IntegerCast<long>(myExpv(rawpp)[indet]);
  }

  void PPMonoidEvOvImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    EXP = myExpv(rawpp)[indet];
  }


  void PPMonoidEvOvImpl::myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const SmallExponent_t* const v = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)
      expv[i] = IntegerCast<long>(v[i]);
  }


  void PPMonoidEvOvImpl::myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const SmallExponent_t* const v = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)  expv[i] = v[i];
  }


  void PPMonoidEvOvImpl::myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    const SmallExponent_t* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] != 0) v[i] = true;
  }


  void PPMonoidEvOvImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
  {
    DivMaskImpl->myAssignFromExpv(dm, myExpv(rawpp), myNumIndets);
  }


  void PPMonoidEvOvImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PPMonoidEvOv(" << myNumIndets << ", " << myOrd << ")";
  }


  void PPMonoidEvOvImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "DEBUG PP: myNumIndets=" << myNumIndets << ", exps=[";
    for (long i=0; i < myNumIndets; ++i)
      out << myExponent(rawpp, i) << " ";
    out << "]" << std::endl;
  }


  PPMonoid NewPPMonoidEvOv(const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
    // Sanity check on the indet names given.
    const long nvars = NumIndets(ord);

    if (len(IndetNames) != nvars)
      CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPPMonoidEvOv(IndetNames,ord)");
    if (!AreDistinct(IndetNames))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidEvOv(IndetNames,ord)");
    if (!AreArityConsistent(IndetNames))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidEvOv(IndetNames,ord)");

    // Inefficient quadratic loop -- speed is probably not important.
    for (long i=0; i < nvars; ++i)
      for (long j=i+1; j < nvars; ++j)
        if (IndetNames[i] == IndetNames[j])
          CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidEvOv(IndetNames,ord)");

    return PPMonoid(new PPMonoidEvOvImpl(IndetNames, ord));
  }

  PPMonoid NewPPMonoidEvOv(const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  {
    return NewPPMonoidEvOv(IndetNames, OrdCtor(len(IndetNames)));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPMonoidEvOv.C,v 1.45 2022/02/18 14:11:56 abbott Exp $
// $Log: PPMonoidEvOv.C,v $
// Revision 1.45  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.44  2021/10/30 19:36:35  abbott
// Summary: Added in some missing "override" (according to clang)
//
// Revision 1.43  2021/10/30 13:47:39  abbott
// Summary: Used keyword override where appropriate (redmine 1625)
//
// Revision 1.42  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.41  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.40  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.39  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.38  2019/03/04 10:31:45  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.37  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.36  2017/12/01 17:29:21  bigatti
// // -- updated Copyright line
// // -- removed doxygen initial comment
// // -- some commented out debugging info
//
// Revision 1.35  2017/04/18 12:50:06  abbott
// Summary: Corrected ifdef use of CoCoA_THREADSAFE_HACK and CoCoA_DEBUG
//
// Revision 1.34  2016/11/03 12:25:25  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.33  2015/12/01 13:11:01  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.32  2015/06/30 12:54:35  abbott
// Summary: Added new fn myIndetsIn
// Author: JAA
//
// Revision 1.31  2015/04/16 20:20:22  abbott
// Summary: Fixed incorrectly edited code
// Author: JAA
//
// Revision 1.30  2015/04/16 16:36:33  abbott
// Summary: Cleaned impls of myPowerOverflowCheck
// Author: JAA
//
// Revision 1.29  2015/04/13 14:42:08  abbott
// Summary: Added myPowerOverflowCheck (1st version)
// Author: JAA
//
// Revision 1.28  2014/07/03 15:36:35  abbott
// Summary: Cleaned up impl of PPMonoids: moved myIndetSymbols & myNumIndets to base class
// Author: JAA
//
// Revision 1.27  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.26  2013/03/15 11:00:50  abbott
// Added check for exponent overflow when powering a PP.
// Merged PPMonoidEv and PPMonoidEvZZ implementations into a single file.
// Implemented new interface for pseudo-ctors for PPMonoidEv which uses a "flag"
// to say whether exponents are big or not.
//
// Revision 1.25  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.24  2012/01/26 16:50:55  bigatti
// -- changed back_inserter into insert
//
// Revision 1.23  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.22  2011/05/03 12:13:12  abbott
// Added static const data member ourMaxExp.
// Code is more readable, and compiler doesn't grumble any more.
//
// Revision 1.21  2011/03/22 22:38:15  abbott
// Fixed some wrong static_casts inside some CoCoA_ASSERTs.
//
// Revision 1.20  2011/03/10 17:27:11  bigatti
// -- changed unsigned long into long in some CoCoA_ASSERT
// -- removed assert in myCmpWDegPartial (done in OrdvArith)
//
// Revision 1.19  2011/03/10 16:39:34  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.18  2010/12/17 16:09:51  abbott
// Corrected used of myIndetSymbols in some assertions.
//
// Revision 1.17  2010/11/30 11:18:11  bigatti
// -- renamed IndetName --> IndetSymbol
//
// Revision 1.16  2010/11/05 16:21:08  bigatti
// -- added ZZExponents
//
// Revision 1.15  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.14  2010/02/03 16:13:52  abbott
// Added new single word tags for specifying the ordering in PPMonoid
// pseudo-ctors.
//
// Revision 1.13  2010/02/02 16:44:31  abbott
// Added radical & IsRadical (via mem fns myRadical & myIsRadical)
// for PPMonoidElems.
//
// Revision 1.12  2009/12/23 18:53:52  abbott
// Major change to conversion functions:
//   convert(..) is now a procedure instead of a function
//   IsConvertible(..) replaces the former convert(..) function
//   Added new NumericCast conversion function (placeholder for BOOST feature)
//   Consequent changes in code which uses these features.
//
// Revision 1.11  2009/09/22 14:01:33  bigatti
// -- added myCmpWDegPartial (ugly name, I know....)
// -- cleaned up and realigned code in PPMonoid*.C files
//
// Revision 1.10  2008/03/26 16:52:04  abbott
// Added exponent overflow checks (also for ordvs) when CoCoA_DEBUG is active.
//
// Revision 1.9  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.8  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.7  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.6  2007/09/25 16:32:30  abbott
// Several minor changes to silence gcc-4.3:
//    more #includes,
//    and fixed a template problemm in RegisterServerOps.C
//
// Revision 1.5  2007/05/31 14:54:31  bigatti
// -- now using AreDistinct and AreArityConsistent for sanity check on
//    indet names
//
// Revision 1.3  2007/05/03 10:35:23  abbott
// Added new PPMonoidEvZZ with (virtually) unlimited exponents.
// Modified test-PPMonoid1.C accordingly.
// Added warning in doc about silent/unchecked exponent overflow in other
// PPMonoids.
//
// Revision 1.2  2007/03/23 18:38:42  abbott
// Separated the "convert" functions (and CheckedCast) into their own files.
// Many consequential changes.  Also corrected conversion to doubles.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.14  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.13  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.12  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.11  2006/12/06 17:35:58  cocoa
// -- style: RawPtr args are now called "raw.."
//
// Revision 1.10  2006/11/27 13:06:23  cocoa
// Anna and Michael made me check without writing a proper message.
//
// Revision 1.9  2006/11/24 17:04:32  cocoa
// -- reorganized includes of header files
//
// Revision 1.8  2006/11/23 17:39:11  cocoa
// -- added #include
//
// Revision 1.7  2006/11/16 11:27:20  cocoa
// -- reinserted myRefCountZero(): sometimes really necessary, in general safe
//
// Revision 1.6  2006/11/14 17:29:20  cocoa
// -- commented out myRefCountZero() (not necessary???)
//
// Revision 1.5  2006/10/16 23:18:59  cocoa
// Corrected use of std::swap and various special swap functions.
// Improved myApply memfn for homs of RingDistrMPolyInlPP.
//
// Revision 1.4  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.2  2006/06/21 17:07:10  cocoa
// Fixed IsIndet bug -- why are there three almost identical copies of code?
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.6  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.5  2006/03/14 17:21:18  cocoa
// Moved concrete PPMonoid impls entirely into their respective .C files.
// Now the corresponding .H files are very compact.
//
// Revision 1.4  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.3  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.2  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.8  2005/10/11 16:37:30  cocoa
// Added new small prime finite field class (see RingFpDouble).
//
// Cleaned makefiles and configuration script.
//
// Tidied PPMonoid code (to eliminate compiler warnings).
//
// Fixed bug in RingFloat::myIsInteger.
//
// Revision 1.7  2005/08/08 16:36:32  cocoa
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
// Revision 1.5  2005/07/08 15:09:28  cocoa
// Added new symbol class (to represent names of indets).
// Integrated the new class into concrete polynomial rings
// and PPMonoid -- many consequential changes.
// Change ctors for the "inline" sparse poly rings: they no
// longer expect a PPMonoid, but build their own instead
// (has to be a PPMonoidOv).
//
// Revision 1.4  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.3  2005/06/23 15:42:41  cocoa
// Fixed typo in GNU fdl -- all doc/*.txt files affected.
// Minor corrections to PPMonoid (discovered while writing doc).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.3  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.4  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.3  2004/11/11 13:41:48  cocoa
// -- change: cout --> GlobalLogput()
//
// Revision 1.2  2004/11/02 14:56:33  cocoa
// -- changed *Print* into *Output* (myPrint --> myOutput)
// -- changed *Var* into *Indet* (myPrintVarName --> myOutputIndetName)
// -- removed suffix "IgnoreDivMask"
// -- added myComputeDivMask
// -- improved storing of IndetNames
// -- changed ExpvElem into SmallExponent_t
//
// Revision 1.1  2004/10/29 15:31:25  cocoa
// -- new PPMonoid for compatibility with OrdvArith (without DivMask)
//
