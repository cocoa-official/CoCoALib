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
  - faster for comparisons especially with matrix defined orderings.

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
