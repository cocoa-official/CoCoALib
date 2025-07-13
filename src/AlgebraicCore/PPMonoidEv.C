//   Copyright (c)  2004-2011,2021  John Abbott and Anna M. Bigatti

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


// Implementation of class PPMonoidEvImpl

#include "CoCoA/PPMonoidEv.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DivMask.H"
#include "CoCoA/MemPool.H"
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

  /*-- class PPMonoidEvSmallExpImpl ---------------------------------------*/
  /*-- class PPMonoidEvImpl ---------------------------------------*/
  /**

  \brief Implementation power product monoid with exponent vector

  PPMonoidEvSmallExpImpl implements a power product monoid for safe testing.
  It stores the exponents as a C vector of SmallExponent_t.
  It is not designed to be efficient (but it might be ;-)

  */
  /*-----------------------------------------------------------------*/

  class PPMonoidEvSmallExpImpl: public PPMonoidBase
  {
    typedef PPMonoidElemRawPtr RawPtr;           // just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; // just to save typing

    static const unsigned long ourMaxExp;        // defined below; value is just numeric_limits<SmallExponent_t>::max()

    class CmpBase
    {
    public:
      CmpBase(long NumIndets, long GradingDim);
      virtual ~CmpBase();
    public:
      virtual int myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const =0;
      virtual void myWDeg(degree& d, const SmallExponent_t* v) const =0;
      //      virtual int myCmpWDeg(const SmallExponent_t* v1, const SmallExponent_t* v2) const =0;
      virtual int myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long) const =0;
      int myCmpWDeg(const SmallExponent_t* v1, const SmallExponent_t* v2) const { return myCmpWDegPartial(v1, v2, myGradingDim); }
      
    protected:
      long myNumIndets;
      long myGradingDim;
    };

    class LexImpl: public CmpBase
    {
    public:
      LexImpl(long NumIndets);
      int myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const override;
      void myWDeg(degree& d, const SmallExponent_t* v) const override;
      int myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long) const override;
    };


    class StdDegLexImpl: public CmpBase
    {
    public:
      StdDegLexImpl(long NumIndets);
      int myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const override;
      void myWDeg(degree& d, const SmallExponent_t* v) const override;
      int myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long) const override;
    };


    class StdDegRevLexImpl: public CmpBase
    {
    public:
      StdDegRevLexImpl(long NumIndets);
      int myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const override;
      void myWDeg(degree& d, const SmallExponent_t* v) const override;
      int myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long) const override;
    };


    class MatrixOrderingImpl: public CmpBase
    {
    public:
      MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      int myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const override;
      void myWDeg(degree& d, const SmallExponent_t* v) const override;
      int myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long) const override;
    private:
      std::vector< std::vector<BigInt> > myOrderMatrix;
    };


  public:
    PPMonoidEvSmallExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    ~PPMonoidEvSmallExpImpl();
  public: // disable copy ctor and assignment
    PPMonoidEvSmallExpImpl(const PPMonoidEvSmallExpImpl& copy) = delete;
    PPMonoidEvSmallExpImpl& operator=(const PPMonoidEvSmallExpImpl& rhs) = delete;

  public:
    void contents() const; // FOR DEBUGGING ONLY

    const std::vector<PPMonoidElem>& myIndets() const override;               ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

    // The functions below are operations on power products owned by PPMonoidEvSmallExpImpl
    const PPMonoidElem& myOne() const override;
    PPMonoidElemRawPtr myNew() const override;                                ///< ctor from nothing
    PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const override;   ///< ctor by assuming ownership
    PPMonoidElemRawPtr myNew(const std::vector<long>& expv) const override;   ///< ctor from exp vector
    void myDelete(RawPtr rawpp) const override;                               ///< dtor, frees pp

    void mySwap(RawPtr rawpp1, RawPtr rawpp2) const override;                 ///< swap(pp1, pp2);
    void myAssignOne(RawPtr rawpp) const override;                            ///< pp = 1
    void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const override;           ///< p = pp1
    void myAssign(RawPtr rawpp, const std::vector<long>& expv) const override;///< pp = expv (assign from exp vector, all exps >= 0)

    void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = pp1*pp2
    void myMulIndetPower(RawPtr rawpp, long indet, long exp) const override;           ///< pp *= indet^exp, assumes exp >= 0
    void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = pp1/pp2
    void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override; ///< pp = pp1/gcd(pp1,pp2)
    void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = gcd(pp1,pp2)
    void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = lcm(pp1,pp2)
    void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const override;                   ///< pp = radical(pp1)
    void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const override;   ///< pp = pp1^exp, assumes exp >= 0
    void myPowerOverflowCheck(ConstRawPtr rawpp1, long exp) const override;            ///< throw if pp1^exp would overflow, assumes exp >= 0

    bool myIsOne(ConstRawPtr rawpp) const override;                                   ///< is pp = 1?
    bool myIsIndet(long& index, ConstRawPtr rawpp) const override;                    ///< true iff pp is an indet
    bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;          ///< are pp1 & pp2 coprime?
    bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;            ///< is pp1 equal to pp2?
    bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;        ///< does pp2 divide pp1?
    bool myIsSqFree(ConstRawPtr rawpp) const override;                                ///< is pp equal to its radical?

    int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;                 ///< -1,0,1 as pp1 < = > pp2
    long myStdDeg(ConstRawPtr rawpp) const override;                                  ///< standard degree of pp
    void myWDeg(degree& d, ConstRawPtr rawpp) const override;                         ///< d = grading(pp)
    int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;             ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long) const override; ///< as myCmpWDeg wrt the first weights
    long myExponent(ConstRawPtr rawpp, long indet) const override;                    ///< exponent of indet in pp
    void myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const override;    ///< EXP = exponent of indet in pp
    void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const override;      ///< expv[i] = exponent(pp,i)
    void myBigExponents(std::vector<BigInt>& v, ConstRawPtr rawpp) const override;    ///< get exponents, SHOULD BE myExponents ???
    void myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const override;          ///< v[i] = true if i-th indet has exponent != 0
    void myOutputSelf(std::ostream& out) const override;                              ///< out << PPM
    // INHERITED DEFINITION of virtual  void myOutput(std::ostream& out, ConstRawPtr rawpp) const;
    void myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const;           ///< print pp in debugging format???


  private: // auxiliary functions
    SmallExponent_t* myExpv(RawPtr) const;
    const SmallExponent_t* myExpv(ConstRawPtr) const;

    void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const override; ///< used by PPWithMask
    bool myCheckExponents(const std::vector<long>& expv) const;

  private: // data members
    ///@name Class members
    //@{
    const long myEntrySize;       ///< size in ?bytes???
    mutable MemPool myMemMgr;     // IMPORTANT: this must come *before* myIndetVector and myOnePtr.
    vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
    unique_ptr<CmpBase> myOrdPtr;         ///< actual implementation of the ordering [should be const???]
    unique_ptr<PPMonoidElem> myOnePtr;
    //@}
  };

  // static variable
  const unsigned long PPMonoidEvSmallExpImpl::ourMaxExp = numeric_limits<SmallExponent_t>::max();


  // dtor not inline for some compiler (?) on PPC
  PPMonoidEvSmallExpImpl::CmpBase::~CmpBase() {}


  // File local inline functions

  inline SmallExponent_t* PPMonoidEvSmallExpImpl::myExpv(RawPtr rawpp) const
  {
    return static_cast<SmallExponent_t*>(rawpp.myRawPtr());
  }


  inline const SmallExponent_t* PPMonoidEvSmallExpImpl::myExpv(ConstRawPtr rawpp) const
  {
    return static_cast<const SmallExponent_t*>(rawpp.myRawPtr());
  }


  bool PPMonoidEvSmallExpImpl::myCheckExponents(const std::vector<long>& expv) const
  {
    // Check expv.size == myNumIndets.
    // Check exps are non-neg and not too big.
    if (len(expv) != myNumIndets) return false;
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] < 0 || static_cast<unsigned long>(expv[i]) > numeric_limits<SmallExponent_t>::max()) return false;
    return true;
  }


  //----   Constructors & destructor   ----//

  PPMonoidEvSmallExpImpl::PPMonoidEvSmallExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord):
      PPMonoidBase(ord, IndetNames),
      myEntrySize(sizeof(SmallExponent_t)*myNumIndets),
      myMemMgr(myEntrySize, "PPMonoidEvSmallExpImpl.myMemMgr"),
      myIndetVector()
  {
    // Put the correct implementation in myOrdPtr...  [VERY UGLY CODE!!!]
    if (IsLex(ord))
      myOrdPtr.reset(new LexImpl(myNumIndets));
    else if (IsStdDegLex(ord))
      myOrdPtr.reset(new StdDegLexImpl(myNumIndets));
    else if (IsStdDegRevLex(ord))
      myOrdPtr.reset(new StdDegRevLexImpl(myNumIndets));
    else
      myOrdPtr.reset(new MatrixOrderingImpl(myNumIndets, GradingDim(ord), OrdMat(ord)));

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


  PPMonoidEvSmallExpImpl::~PPMonoidEvSmallExpImpl()
  {}

/////////////////////////////////////////////////////////////////////////////



  const std::vector<PPMonoidElem>& PPMonoidEvSmallExpImpl::myIndets() const
  {
    return myIndetVector;
  }


  const PPMonoidElem& PPMonoidEvSmallExpImpl::myOne() const
  {
    return *myOnePtr;
  }


  PPMonoidElemRawPtr PPMonoidEvSmallExpImpl::myNew() const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssignOne(rawpp); // cannot throw
    return rawpp;
  }

  PPMonoidElemRawPtr PPMonoidEvSmallExpImpl::myNew(PPMonoidElemConstRawPtr rawcopypp) const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssign(rawpp, rawcopypp); // cannot throw
    return rawpp;
  }


  PPMonoidElemRawPtr PPMonoidEvSmallExpImpl::myNew(const std::vector<long>& expv) const
  {
    CoCoA_ASSERT(myCheckExponents(expv));
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    myAssign(rawpp, expv); // cannot throw
    return rawpp;
  }


  void PPMonoidEvSmallExpImpl::myAssignOne(RawPtr rawpp) const
  {
    SmallExponent_t* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = 0;
  }


  void PPMonoidEvSmallExpImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    if (rawpp == rawpp1) return;

//     SmallExponent_t* const expv = myExpv(rawpp);
//     const SmallExponent_t* const expv1 = myExpv(rawpp1);
//     std::copy(expv1, expv1+myNumIndets, expv);
    memcpy(myExpv(rawpp), myExpv(rawpp1), myNumIndets*sizeof(SmallExponent_t));
  }

  void PPMonoidEvSmallExpImpl::myAssign(RawPtr rawpp, const vector<long>& expv1) const
  {
    CoCoA_ASSERT(myCheckExponents(expv1));

    SmallExponent_t* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = IntegerCast<SmallExponent_t>(expv1[i]);
  }


  void PPMonoidEvSmallExpImpl::myDelete(RawPtr rawpp) const
  {
    myMemMgr.free(rawpp.myRawPtr());
  }


  void PPMonoidEvSmallExpImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
  {
    if (rawpp1 == rawpp2) return;
    SmallExponent_t* const expv1 = myExpv(rawpp1);
    SmallExponent_t* const expv2 = myExpv(rawpp2);
    for (long i = 0; i < myNumIndets; ++i)
      std::swap(expv1[i], expv2[i]);
  }


  void PPMonoidEvSmallExpImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
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
  }


  void PPMonoidEvSmallExpImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    SmallExponent_t* const expv = myExpv(rawpp);
    // If CoCoA_DEBUG active, check for exponent overflow...
    CoCoA_ASSERT("Exponent Overflow" && ourMaxExp - expv[indet] >= static_cast<unsigned long>(exp));
    expv[indet] += static_cast<SmallExponent_t>(exp);  // cast to keep M$ compiler quiet
  }


  void PPMonoidEvSmallExpImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
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
  }


  void PPMonoidEvSmallExpImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] > expv2[i])
        expv[i] = expv1[i] - expv2[i];
      else
        expv[i] = 0;
  }


  void PPMonoidEvSmallExpImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = min(expv1[i], expv2[i]);
  }


  void PPMonoidEvSmallExpImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = max(expv1[i], expv2[i]);
  }


  void PPMonoidEvSmallExpImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = (expv1[i] > 0);
  }


  void PPMonoidEvSmallExpImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long LongExp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(LongExp >= 0);
#ifdef CoCoA_DEBUG
    myPowerOverflowCheck(rawpp1, LongExp);
#endif
    if (static_cast<unsigned long>(LongExp) > ourMaxExp)
      CoCoA_THROW_ERROR1(ERR::ExpTooBig);
    const SmallExponent_t exp = static_cast<SmallExponent_t>(LongExp);

    SmallExponent_t* const expv = myExpv(rawpp);
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = exp * expv1[i];
  }


  void PPMonoidEvSmallExpImpl::myPowerOverflowCheck(ConstRawPtr rawpp, long LongExp) const
  {
    if (LongExp == 0 || LongExp == 1) return;
    CoCoA_ASSERT(LongExp >= 0);
    if (static_cast<unsigned long>(LongExp) > ourMaxExp)
      CoCoA_THROW_ERROR1(ERR::ExpTooBig);
    const SmallExponent_t exp = static_cast<SmallExponent_t>(LongExp);
    const SmallExponent_t limit = ourMaxExp/exp;

    const SmallExponent_t* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (expv[i] > limit)
        CoCoA_THROW_ERROR1(ERR::ExpTooBig);
    }
  }


  bool PPMonoidEvSmallExpImpl::myIsOne(ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv[i] != 0) return false;

    return true;
  }


  bool PPMonoidEvSmallExpImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
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


  bool PPMonoidEvSmallExpImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] != 0 && expv2[i] != 0) return false;

    return true;
  }


  bool PPMonoidEvSmallExpImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] != expv2[i]) return false;
    return true;
  }


  bool PPMonoidEvSmallExpImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const SmallExponent_t* const expv1 = myExpv(rawpp1);
    const SmallExponent_t* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] < expv2[i]) return false;

    return true;
  }


  bool PPMonoidEvSmallExpImpl::myIsSqFree(ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv[i] > 1) return false;

    return true;
  }


  int PPMonoidEvSmallExpImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdPtr->myCmpExpvs(myExpv(rawpp1), myExpv(rawpp2));
  }


  long PPMonoidEvSmallExpImpl::myStdDeg(ConstRawPtr rawpp) const
  {
    const SmallExponent_t* const expv = myExpv(rawpp);
    long d=0;
    for (long i=0; i < myNumIndets; ++i)
      d += expv[i];
    return d;
  }


  void PPMonoidEvSmallExpImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
  {
    myOrdPtr->myWDeg(d, myExpv(rawpp));
  }


  int PPMonoidEvSmallExpImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdPtr->myCmpWDeg(myExpv(rawpp1), myExpv(rawpp2));
  }


  int PPMonoidEvSmallExpImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long i) const
  {
    return myOrdPtr->myCmpWDegPartial(myExpv(rawpp1), myExpv(rawpp2), i);
  }


  long PPMonoidEvSmallExpImpl::myExponent(ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(indet < myNumIndets);
    return IntegerCast<long>(myExpv(rawpp)[indet]);
  }

  void PPMonoidEvSmallExpImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(indet < myNumIndets);
    EXP = myExpv(rawpp)[indet];
  }


  void PPMonoidEvSmallExpImpl::myExponents(std::vector<long>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    const SmallExponent_t* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)
      v[i] = IntegerCast<long>(expv[i]);
  }


  void PPMonoidEvSmallExpImpl::myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const SmallExponent_t* const v = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)  expv[i] = v[i];
  }


  void PPMonoidEvSmallExpImpl::myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    const SmallExponent_t* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] != 0)
        v[i] = true;
  }

  void PPMonoidEvSmallExpImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
  {
    DivMaskImpl->myAssignFromExpv(dm, myExpv(rawpp), myNumIndets);
  }


  void PPMonoidEvSmallExpImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PPMonoidEv(" << myNumIndets << ", " << myOrd <<")";
  }


  void PPMonoidEvSmallExpImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "DEBUG PP: myNumIndets=" << myNumIndets << ", exps=[";
    for (long i=0; i < myNumIndets; ++i)
      out << myExponent(rawpp, i) << " ";
    out << "]" << std::endl;
  }


//--- orderings ----------------------------------------

  PPMonoidEvSmallExpImpl::CmpBase::CmpBase(long NumIndets, long GradingDim):
    myNumIndets(NumIndets),
    myGradingDim(GradingDim)
  {
    CoCoA_ASSERT(NumIndets > 0);
    CoCoA_ASSERT(NumIndets < 1000000); // complain about ridiculously large number of indets
    CoCoA_ASSERT(GradingDim >= 0);
    CoCoA_ASSERT(GradingDim <= NumIndets);
  }


//--- LexImpl

  PPMonoidEvSmallExpImpl::LexImpl::LexImpl(long NumIndets):
    CmpBase(NumIndets, 0)
  {
  }


  int PPMonoidEvSmallExpImpl::LexImpl::myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const
  {
    for (long i=0; i<myNumIndets; ++i)
      if (v1[i] != v2[i])
      {
        if (v1[i] > v2[i]) return 1; else return -1;
      }
    return 0;
  }


  void PPMonoidEvSmallExpImpl::LexImpl::myWDeg(degree& /*d*/, const SmallExponent_t* /*v*/) const
  {
    // deliberately does nothing because GradingDim=0
    //??? should it assign:  d = []
  }


  int PPMonoidEvSmallExpImpl::LexImpl::myCmpWDegPartial(const SmallExponent_t* /*v1*/, const SmallExponent_t* /*v2*/, long /*i*/) const
  {
    return 0; // GradingDim=0, and all degrees in Z^0 are equal
  }

//--- StdDegLexImpl

  PPMonoidEvSmallExpImpl::StdDegLexImpl::StdDegLexImpl(long NumIndets):
    CmpBase(NumIndets, 1)
  {
  }


  void PPMonoidEvSmallExpImpl::StdDegLexImpl::myWDeg(degree& d, const SmallExponent_t* v) const
  {
    long deg=0;
    for (long i=0; i < myNumIndets; ++i)
      deg += v[i];
    SetComponent(d, 0, deg);
  }


  int PPMonoidEvSmallExpImpl::StdDegLexImpl::myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const
  {
    long deg1=0, deg2=0;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2)
    {
      if (deg1 > deg2) return 1; else return -1;
    }

    for (long i=0; i < myNumIndets; ++i)
      if (v1[i] != v2[i])
      {
        if (v1[i] > v2[i]) return 1; else return -1;
      }
    return 0;
  }


  int PPMonoidEvSmallExpImpl::StdDegLexImpl::myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long i) const
  {
    if (i==0) return 0; // ie: wrt to 0-grading
    long deg1=0, deg2=0;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2)
    {
      if (deg1 > deg2) return 1; else return -1;
    }
    return 0;
  }


//--- StdDegRevLexImpl

  PPMonoidEvSmallExpImpl::StdDegRevLexImpl::StdDegRevLexImpl(long NumIndets):
    CmpBase(NumIndets, 1)
  {
  }


  void PPMonoidEvSmallExpImpl::StdDegRevLexImpl::myWDeg(degree& d, const SmallExponent_t* v) const
  {
    long deg=0;
    for (long i=0; i < myNumIndets; ++i)
      deg += v[i];
    SetComponent(d, 0, deg);
  }


  int PPMonoidEvSmallExpImpl::StdDegRevLexImpl::myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const
  {
    long deg1=0, deg2=0;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2)
    {
      if (deg1 > deg2) return 1; else return -1;
    }

    for (long i = myNumIndets-1; i>0; --i)
      if (v1[i] != v2[i])
      {
        if (v1[i] > v2[i]) return -1; else return 1;
      }
    return 0;
  }


  int PPMonoidEvSmallExpImpl::StdDegRevLexImpl::myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long PartialGrDim) const
  {
    CoCoA_ASSERT(PartialGrDim==0 || PartialGrDim==1);
    if (PartialGrDim==0) return 0; // ie: wrt to 0-grading
    long deg1=0, deg2=0;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2)
    {
      if (deg1 > deg2) return 1; else return -1;
    }
    return 0;
  }


//--- MatrixOrderingImpl

  PPMonoidEvSmallExpImpl::MatrixOrderingImpl::MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix):
    CmpBase(NumIndets, GradingDim)
  {
    CoCoA_ASSERT(NumIndets > 0);
    CoCoA_ASSERT(0 <= GradingDim && GradingDim <= NumIndets);
    CoCoA_ASSERT(NumRows(OrderMatrix) == NumIndets);
    CoCoA_ASSERT(NumCols(OrderMatrix) == NumIndets);

    myOrderMatrix.resize(NumIndets, vector<BigInt>(NumIndets));
    for (long i=0; i < NumIndets; ++i)
      for (long j=0; j < NumIndets; ++j)
        myOrderMatrix[i][j] = ConvertTo<BigInt>(OrderMatrix(i,j));
  }


  void PPMonoidEvSmallExpImpl::MatrixOrderingImpl::myWDeg(degree& d, const SmallExponent_t* v) const
  {
    BigInt deg;
    for (long g=0; g < myGradingDim; ++g)
    {
      deg = 0;
      for (long i=0; i < myNumIndets; ++i)
        deg += v[i] * myOrderMatrix[g][i];
      SetComponent(d, g, deg);
    }
  }


  static int MatrixOrderingCmpArrays(const vector<vector<BigInt> >& M, const SmallExponent_t* v1, const SmallExponent_t* v2, long UpTo)
  {
    CoCoA_ASSERT(0 <= UpTo && UpTo <= len(M));
    BigInt s1, s2; // automatically set to 0
    const long ncols = len(M[0]);
    for (long r=0; r < UpTo; ++r)
    {
      for (long i=0; i < ncols; ++i)
      {
        s1 += v1[i]*M[r][i];
        s2 += v2[i]*M[r][i];
      }
      if (s1 != s2)
      {
        if (s1 > s2) return 1; else return -1;
      }
      s1 = 0;
      s2 = 0;
    }
    return 0;
  }


  int PPMonoidEvSmallExpImpl::MatrixOrderingImpl::myCmpExpvs(const SmallExponent_t* v1, const SmallExponent_t* v2) const
  {
    return MatrixOrderingCmpArrays(myOrderMatrix, v1, v2, myNumIndets);
  }


//   int PPMonoidEvSmallExpImpl::MatrixOrderingImpl::myCmpWDeg(const SmallExponent_t* v1, const SmallExponent_t* v2) const
//   {
//     return MatrixOrderingCmpArrays(myOrderMatrix, v1, v2, myGradingDim);
//   }


  int PPMonoidEvSmallExpImpl::MatrixOrderingImpl::myCmpWDegPartial(const SmallExponent_t* v1, const SmallExponent_t* v2, long PartialGrDim) const
  {
    return MatrixOrderingCmpArrays(myOrderMatrix, v1, v2, PartialGrDim);
  }


  //////////////////////////////////////////////////////////////////
  // BIG EXPONENTS

  class PPMonoidEvBigExpImpl: public PPMonoidBase
  {
    typedef PPMonoidElemRawPtr RawPtr;           // just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; // just to save typing

    class CmpBase
    {
    public:
      CmpBase(long NumIndets, long GradingDim);
      virtual ~CmpBase() {};
    public:
      virtual int myCmpExpvs(const BigInt* v1, const BigInt* v2) const =0;
      virtual void myWDeg(degree& d, const BigInt* v) const =0;
      virtual int myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const =0;
      //      virtual int myCmpWDeg(const BigInt* v1, const BigInt* v2) const =0;
      int myCmpWDeg(const BigInt* v1, const BigInt* v2) const { return myCmpWDegPartial(v1, v2, myGradingDim); }
    protected:
      long myNumIndets;
      long myGradingDim;
    };

    class LexImpl: public CmpBase
    {
    public:
      LexImpl(long NumIndets);
      int myCmpExpvs(const BigInt* v1, const BigInt* v2) const override;
      void myWDeg(degree& d, const BigInt* v) const override;
      int myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const override;
    };


    class StdDegLexImpl: public CmpBase
    {
    public:
      StdDegLexImpl(long NumIndets);
      int myCmpExpvs(const BigInt* v1, const BigInt* v2) const override;
      void myWDeg(degree& d, const BigInt* v) const override;
      int myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const override;
    };


    class StdDegRevLexImpl: public CmpBase
    {
    public:
      StdDegRevLexImpl(long NumIndets);
      int myCmpExpvs(const BigInt* v1, const BigInt* v2) const override;
      void myWDeg(degree& d, const BigInt* v) const override;
      int myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const override;
    };


    class MatrixOrderingImpl: public CmpBase
    {
    public:
      MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix);
      int myCmpExpvs(const BigInt* v1, const BigInt* v2) const override;
      void myWDeg(degree& d, const BigInt* v) const override;
      int myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const override;
    private:
      std::vector< std::vector<BigInt> > myOrderMatrix;
    };


  public:
    PPMonoidEvBigExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    ~PPMonoidEvBigExpImpl();
  public: // disable copy ctor and assignment
    PPMonoidEvBigExpImpl(const PPMonoidEvBigExpImpl& copy) = delete;
    PPMonoidEvBigExpImpl& operator=(const PPMonoidEvBigExpImpl& rhs) = delete;

  public:
    void contents() const; // FOR DEBUGGING ONLY

    const std::vector<PPMonoidElem>& myIndets() const override;               ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

    // The functions below are operations on power products owned by PPMonoidEvBigExpImpl
    const PPMonoidElem& myOne() const override;
    PPMonoidElemRawPtr myNew() const override;                                ///< ctor from nothing
    PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const override;   ///< ctor by assuming ownership
    PPMonoidElemRawPtr myNew(const std::vector<long>& expv) const override;   ///< ctor from exp vector
    PPMonoidElemRawPtr myNew(const std::vector<BigInt>& EXPV) const override; ///< ctor from exp vector
    void myDelete(RawPtr rawpp) const override;                               ///< dtor, frees pp

    void mySwap(RawPtr rawpp1, RawPtr rawpp2) const override;                 ///< swap(pp1, pp2);
    void myAssignOne(RawPtr rawpp) const override;                            ///< pp = 1
    void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const override;           ///< p = pp1
    void myAssign(RawPtr rawpp, const std::vector<long>& expv) const override;///< pp = expv (assign from exp vector)

    void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = pp1*pp2
    void myMulIndetPower(RawPtr rawpp, long indet, long exp) const override;           ///< pp *= indet^exp, assumes exp >= 0
    void myMulIndetPower(RawPtr rawpp, long indet, const BigInt& EXP) const override;  ///< pp *= indet^EXP
    void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = pp1/pp2
    void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override; ///< pp = pp1/gcd(pp1,pp2)
    void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = gcd(pp1,pp2)
    void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< pp = lcm(pp1,pp2)
    void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const override;                   ///< pp = radical(pp1)
    void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const override;   ///< pp = pp1^exp, assumes exp >= 0
    void myPowerOverflowCheck(ConstRawPtr rawpp1, long exp) const override;            ///< throw if pp1^exp would overflow, assumes exp >= 0

    bool myIsOne(ConstRawPtr rawpp) const override;                                ///< is pp = 1?
    bool myIsIndet(long& index, ConstRawPtr rawpp) const override;                 ///< true iff pp is an indet
    bool myIsIndetPosPower(long& indet, BigInt& EXP, ConstRawPtr rawpp) const override;
    bool myIsIndetPosPower(long& indet, long& pow, ConstRawPtr rawpp) const override;
    bool myIsIndetPosPower(ConstRawPtr rawpp) const override;    
    bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;       ///< are pp1 & pp2 coprime?
    bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;         ///< is pp1 equal to pp2?
    bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;     ///< does pp2 divide pp1?
    bool myIsSqFree(ConstRawPtr rawpp) const override;                             ///< is pp equal to its radical?

    int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;              ///< -1,0,1 as pp1 < = > pp2
    long myStdDeg(ConstRawPtr rawpp) const override;                               ///< standard degree of pp
    void myWDeg(degree& d, ConstRawPtr rawpp) const override;                      ///< d = grading(pp)
    int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;          ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long GrDim ) const override; ///< as myCmpWDeg wrt the first GrDim weights
    long myExponent(ConstRawPtr rawpp, long indet) const override;                 ///< exponent of indet in pp
    void myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const override; ///< EXP = exponent of indet in pp
    void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const override;   ///< expv[i] = exponent(pp,i)
    void myBigExponents(std::vector<BigInt>& v, ConstRawPtr rawpp) const override; ///< get exponents, SHOULD BE myExponents ???
    void myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const override;       ///< v[i] = true if i-th indet has exponent != 0
    void myOutputSelf(std::ostream& out) const override;                           ///< out << PPM
    // INHERITED DEFINITION of virtual  void myOutput(std::ostream& out, ConstRawPtr rawpp) const;
    void myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const;        ///< print pp in debugging format???


  private: // auxiliary functions
    BigInt* myExpv(RawPtr) const;
    const BigInt* myExpv(ConstRawPtr) const;

    void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const override; ///< used by PPWithMask
    bool myCheckExponents(const std::vector<long>& expv) const;
    bool myCheckExponents(const std::vector<BigInt>& EXPV) const;

  private: // data members
    ///@name Class members
    //@{
    const long myEntrySize;     ///< size in ????
    mutable MemPool myMemMgr;     // IMPORTANT: this must come *before* myIndetVector and myOnePtr.
    std::vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
    std::unique_ptr<CmpBase> myOrdPtr;         ///< actual implementation of the ordering [should be const???]
    std::unique_ptr<PPMonoidElem> myOnePtr;
    //@}
  };


  // File local inline functions

  inline BigInt* PPMonoidEvBigExpImpl::myExpv(RawPtr rawpp) const
  {
    return static_cast<BigInt*>(rawpp.myRawPtr());
  }


  inline const BigInt* PPMonoidEvBigExpImpl::myExpv(ConstRawPtr rawpp) const
  {
    return static_cast<const BigInt*>(rawpp.myRawPtr());
  }


  bool PPMonoidEvBigExpImpl::myCheckExponents(const std::vector<long>& expv) const
  {
    // Check len(expv) == myNumIndets.
    // Check exps are non-neg and not too big.
    if (len(expv) != myNumIndets) return false;
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] < 0) return false;
    return true;
  }

  bool PPMonoidEvBigExpImpl::myCheckExponents(const std::vector<BigInt>& expv) const
  {
    // Check len(expv) == myNumIndets.
    // Check exps are non-neg and not too big.
    if (len(expv) != myNumIndets) return false;
    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] < 0) return false;
    return true;
  }


  //----   Constructors & destructor   ----//

  PPMonoidEvBigExpImpl::PPMonoidEvBigExpImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord):
      PPMonoidBase(ord, IndetNames),
      myEntrySize(sizeof(BigInt)*NumIndets(ord)),
      myMemMgr(myEntrySize, "PPMonoidEvBigExpImpl.myMemMgr"),
      myIndetVector()
  {
    // Put the correct implementation in myOrdPtr...  [VERY UGLY CODE!!!]
    if (IsLex(ord))
      myOrdPtr.reset(new LexImpl(myNumIndets));
    else if (IsStdDegLex(ord))
      myOrdPtr.reset(new StdDegLexImpl(myNumIndets));
    else if (IsStdDegRevLex(ord))
      myOrdPtr.reset(new StdDegRevLexImpl(myNumIndets));
    else
      myOrdPtr.reset(new MatrixOrderingImpl(myNumIndets, GradingDim(ord), OrdMat(ord)));

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


  PPMonoidEvBigExpImpl::~PPMonoidEvBigExpImpl()
  {}

/////////////////////////////////////////////////////////////////////////////



  const std::vector<PPMonoidElem>& PPMonoidEvBigExpImpl::myIndets() const
  {
    return myIndetVector;
  }


  const PPMonoidElem& PPMonoidEvBigExpImpl::myOne() const
  {
    return *myOnePtr;
  }


  PPMonoidElemRawPtr PPMonoidEvBigExpImpl::myNew() const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    BigInt* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i) // BUG BUG: if this throws, rawpp is leaked
      new(&expv[i]) BigInt;
    myAssignOne(rawpp); // cannot throw
    return rawpp;
  }

  PPMonoidElemRawPtr PPMonoidEvBigExpImpl::myNew(PPMonoidElemConstRawPtr rawcopypp) const
  {
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    BigInt* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i) // BUG BUG: if this throws, rawpp is leaked
      new(&expv[i]) BigInt;
    myAssign(rawpp, rawcopypp); // cannot throw
    return rawpp;
  }


  PPMonoidElemRawPtr PPMonoidEvBigExpImpl::myNew(const std::vector<long>& expv) const
  {
    CoCoA_ASSERT(myCheckExponents(expv));
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    BigInt* const expv1 = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i) // BUG BUG: if this throws, rawpp is leaked
      new(&expv1[i]) BigInt;
    myAssign(rawpp, expv); // cannot throw
    return rawpp;
  }


  PPMonoidElemRawPtr PPMonoidEvBigExpImpl::myNew(const std::vector<BigInt>& EXPV) const
  {
    CoCoA_ASSERT(myCheckExponents(EXPV));
    PPMonoidElemRawPtr rawpp(myMemMgr.alloc());
    BigInt* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i) // BUG BUG: if this throws, rawpp is leaked
      new(&expv[i]) BigInt;
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = EXPV[i];
    return rawpp;
  }


  void PPMonoidEvBigExpImpl::myAssignOne(RawPtr rawpp) const
  {
    BigInt* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = 0;
  }


  void PPMonoidEvBigExpImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    if (rawpp == rawpp1) return;

    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);
    std::copy(&expv1[0], &expv1[myNumIndets], &expv[0]);
  }

  void PPMonoidEvBigExpImpl::myAssign(RawPtr rawpp, const vector<long>& expv1) const
  {
    CoCoA_ASSERT(myCheckExponents(expv1));

    BigInt* const expv = myExpv(rawpp);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = expv1[i];
  }


  void PPMonoidEvBigExpImpl::myDelete(RawPtr rawpp) const
  {
    BigInt* const expv = myExpv(rawpp);;
    for (long i=0; i < myNumIndets; ++i)
      expv[i].~BigInt();
    myMemMgr.free(rawpp.myRawPtr());
  }


  void PPMonoidEvBigExpImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
  {
    if (rawpp1 == rawpp2) return;
    BigInt* const expv1 = myExpv(rawpp1);
    BigInt* expv2 = myExpv(rawpp2);
    for (long i = 0; i < myNumIndets; ++i)
      std::swap(expv1[i], expv2[i]);
  }


  void PPMonoidEvBigExpImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = expv1[i] + expv2[i];
  }


  void PPMonoidEvBigExpImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(exp >= 0);
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    BigInt* const expv = myExpv(rawpp);
    expv[indet] += exp;
  }


  // This fn overrides the default defn in PPMonoid.C (which throws for big exps)
  void PPMonoidEvBigExpImpl::myMulIndetPower(RawPtr rawpp, long indet, const BigInt& EXP) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    CoCoA_ASSERT(EXP >= 0);
    BigInt* const expv = myExpv(rawpp);
    expv[indet] += EXP;
  }


  void PPMonoidEvBigExpImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);
    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = expv1[i] - expv2[i]; // assert that result is positive?
  }


  void PPMonoidEvBigExpImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] > expv2[i])
        expv[i] = expv1[i] - expv2[i];
      else
        expv[i] = 0;
  }


  void PPMonoidEvBigExpImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = min(expv1[i], expv2[i]);
  }


  void PPMonoidEvBigExpImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    // No worries about aliasing.
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = max(expv1[i], expv2[i]);
  }


  void PPMonoidEvBigExpImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = (expv1[i] > 0);
  }


  void PPMonoidEvBigExpImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const  // assumes exp >= 0
  {
    CoCoA_ASSERT(exp >= 0);
    BigInt* const expv = myExpv(rawpp);
    const BigInt* const expv1 = myExpv(rawpp1);

    for (long i = 0; i < myNumIndets; ++i)
      expv[i] = exp * expv1[i];
  }


  void PPMonoidEvBigExpImpl::myPowerOverflowCheck(ConstRawPtr /*rawpp*/, long exp) const
  {
    (void)(exp); // just to avoid compiler warning
    CoCoA_ASSERT(exp >= 0);
    // nothing to check -- exps here cannot overflow :-)
  }


  bool PPMonoidEvBigExpImpl::myIsOne(ConstRawPtr rawpp) const
  {
    const BigInt* const expv = myExpv(rawpp);

    for (long i = 0; i < myNumIndets; ++i)
      if (!IsZero(expv[i])) return false;

    return true;
  }


  bool PPMonoidEvBigExpImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
  {
    const BigInt* const expv = myExpv(rawpp);
    long j = myNumIndets;
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (IsZero(expv[i])) continue;
      if (j != myNumIndets || !IsOne(expv[i])) return false;
      j = i;
    }
    if (j == myNumIndets) return false;
    index = j;
    return true;
  }


  bool PPMonoidEvBigExpImpl::myIsIndetPosPower(long& index, BigInt& EXP, ConstRawPtr rawpp) const
  {
    const BigInt* const expv = myExpv(rawpp);
    long j = myNumIndets;
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (IsZero(expv[i])) continue;
      if (j != myNumIndets) return false;
      j = i;
    }
    if (j == myNumIndets) return false;
    index = j;
    EXP = expv[index];
    return true;
  }

  bool PPMonoidEvBigExpImpl::myIsIndetPosPower(long& index, long& pow, ConstRawPtr rawpp) const
  {
    BigInt POW;
    long TmpIndex, TmpPow;
    if (!myIsIndetPosPower(TmpIndex, POW, rawpp)) return false;
    if (!IsConvertible(TmpPow, POW))
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    index = TmpIndex;
    pow = TmpPow;
    return true;
  }

  bool PPMonoidEvBigExpImpl::myIsIndetPosPower(ConstRawPtr rawpp) const
  {
    const BigInt* const expv = myExpv(rawpp);
    long j = myNumIndets;
    for (long i = 0; i < myNumIndets; ++i)
    {
      if (IsZero(expv[i])) continue;
      if (j != myNumIndets) return false;
      j = i;
    }
    return j != myNumIndets;
  }


  bool PPMonoidEvBigExpImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (!IsZero(expv1[i]) && !IsZero(expv2[i])) return false;

    return true;
  }


  bool PPMonoidEvBigExpImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] != expv2[i]) return false;
    return true;
  }


  bool PPMonoidEvBigExpImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const BigInt* const expv1 = myExpv(rawpp1);
    const BigInt* const expv2 = myExpv(rawpp2);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv1[i] < expv2[i]) return false;

    return true;
  }


  bool PPMonoidEvBigExpImpl::myIsSqFree(ConstRawPtr rawpp) const
  {
    const BigInt* const expv = myExpv(rawpp);

    for (long i = 0; i < myNumIndets; ++i)
      if (expv[i] > 1) return false;

    return true;
  }


  int PPMonoidEvBigExpImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdPtr->myCmpExpvs(myExpv(rawpp1), myExpv(rawpp2));
  }


  long PPMonoidEvBigExpImpl::myStdDeg(ConstRawPtr rawpp) const
  {
    const BigInt* const expv = myExpv(rawpp);
    BigInt d;
    for (long i=0; i < myNumIndets; ++i)
      d += expv[i];
    long ans;
    if (!IsConvertible(ans, d))
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return ans;
  }


  void PPMonoidEvBigExpImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
  {
    myOrdPtr->myWDeg(d, myExpv(rawpp));
  }


  int PPMonoidEvBigExpImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myOrdPtr->myCmpWDeg(myExpv(rawpp1), myExpv(rawpp2));
  }


  int PPMonoidEvBigExpImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long GrDim) const
  {
    return myOrdPtr->myCmpWDegPartial(myExpv(rawpp1), myExpv(rawpp2), GrDim);
  }


  long PPMonoidEvBigExpImpl::myExponent(ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    long ans;
    if (!IsConvertible(ans, myExpv(rawpp)[indet]))
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return ans;
  }

  void PPMonoidEvBigExpImpl::myBigExponent(BigInt& EXP, ConstRawPtr rawpp, long indet) const
  {
    CoCoA_ASSERT(0 <= indet && indet < myNumIndets);
    EXP = myExpv(rawpp)[indet];
  }


  void PPMonoidEvBigExpImpl::myExponents(vector<long>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    const BigInt* const expv = myExpv(rawpp);
    bool OK = true;
    for (long i=0; i < myNumIndets; ++i)
    {
      if (!IsConvertible(v[i], expv[i]))
      {
        OK = false;
        if (expv[i] > 0)
          v[i] = numeric_limits<long>::max();
        else
          v[i] = numeric_limits<long>::min();
      }
      if (!OK)
        CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    }
  }


  void PPMonoidEvBigExpImpl::myBigExponents(vector<BigInt>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const BigInt* const v = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)  expv[i] = v[i];
  }
  

  void PPMonoidEvBigExpImpl::myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    const BigInt* const expv = myExpv(rawpp);
    for (long i=0; i < myNumIndets; ++i)
      if (!IsZero(expv[i]))
        v[i] = true;
  }


  void PPMonoidEvBigExpImpl::myComputeDivMask(DivMask& /*dm*/, const DivMaskRule& /*DivMaskImpl*/, ConstRawPtr /*rawpp*/) const
  {
    CoCoA_THROW_ERROR1(ERR::NYI);
// BUG BUG can't do this (yet)    DivMaskImpl->myAssignFromExpv(dm, myExpv(rawpp), myNumIndets);
  }


  void PPMonoidEvBigExpImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PPMonoidBigEv(" << myNumIndets << ", " << myOrd <<")";
  }


  void PPMonoidEvBigExpImpl::myDebugPrint(std::ostream& out, ConstRawPtr rawpp) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "DEBUG PP: myNumIndets=" << myNumIndets << ", exps=[";
    for (long i=0; i < myNumIndets; ++i)
      out << myExponent(rawpp, i) << " ";
    out << "]" << std::endl;
  }


//--- orderings ----------------------------------------

  PPMonoidEvBigExpImpl::CmpBase::CmpBase(long NumIndets, long GradingDim):
    myNumIndets(NumIndets),
    myGradingDim(GradingDim)
  {
    CoCoA_ASSERT(NumIndets > 0);
    CoCoA_ASSERT(NumIndets < 1000000); // complain about ridiculously large number of indets
    CoCoA_ASSERT(0 <= GradingDim && GradingDim <= NumIndets);
  }


//--- LexImpl

  PPMonoidEvBigExpImpl::LexImpl::LexImpl(long NumIndets):
    CmpBase(NumIndets, 0)
  {
  }


  int PPMonoidEvBigExpImpl::LexImpl::myCmpExpvs(const BigInt* v1, const BigInt* v2) const
  {
    for (long i=0; i<myNumIndets; ++i)
      if (v1[i] != v2[i]) return (v1[i]>v2[i] ? 1 : -1 );
    return 0;
  }


  void PPMonoidEvBigExpImpl::LexImpl::myWDeg(degree& /*d*/, const BigInt* /*v*/) const
  {
    // deliberately does nothing because GradingDim=0
    //??? should it assign:  d = []
  }


  int PPMonoidEvBigExpImpl::LexImpl::myCmpWDegPartial(const BigInt* /*v1*/, const BigInt* /*v2*/, long /*GrDim*/) const
  {
    return 0; // GradingDim=0, and all degrees in Z^0 are equal
  }

//--- StdDegLexImpl

  PPMonoidEvBigExpImpl::StdDegLexImpl::StdDegLexImpl(long NumIndets):
    CmpBase(NumIndets, 1)
  {
  }


  void PPMonoidEvBigExpImpl::StdDegLexImpl::myWDeg(degree& d, const BigInt* v) const
  {
    BigInt deg;
    for (long i=0; i < myNumIndets; ++i)
      deg += v[i];
    SetComponent(d, 0, deg);
  }


  int PPMonoidEvBigExpImpl::StdDegLexImpl::myCmpExpvs(const BigInt* v1, const BigInt* v2) const
  {
    BigInt deg1;
    BigInt deg2;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2) return (deg1>deg2 ? 1 : -1 );

    for (long i=0; i < myNumIndets; ++i)
      if (v1[i] != v2[i]) return (v1[i]>v2[i] ? 1 : -1 );
    return 0;
  }


  int PPMonoidEvBigExpImpl::StdDegLexImpl::myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const
  {
    CoCoA_ASSERT(0 <= GrDim && GrDim <= 1);
    if (GrDim == 0) return 0; // ie: wrt to 0-grading
    BigInt deg1;
    BigInt deg2;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2) return (deg1>deg2 ? 1 : -1 );
    return 0;
  }


//--- StdDegRevLexImpl

  PPMonoidEvBigExpImpl::StdDegRevLexImpl::StdDegRevLexImpl(long NumIndets):
    CmpBase(NumIndets, 1)
  {
  }


  void PPMonoidEvBigExpImpl::StdDegRevLexImpl::myWDeg(degree& d, const BigInt* v) const
  {
    BigInt deg;
    for (long i=0; i < myNumIndets; ++i)
      deg += v[i];
    SetComponent(d, 0, deg);
  }


  int PPMonoidEvBigExpImpl::StdDegRevLexImpl::myCmpExpvs(const BigInt* v1, const BigInt* v2) const
  {
    BigInt deg1;
    BigInt deg2;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2) return (deg1>deg2 ? 1 : -1 );

    for (long i = myNumIndets-1; i>0; --i)
      if (v1[i] != v2[i]) return (v1[i]>v2[i] ? -1 : 1 );
    return 0;
  }


  int PPMonoidEvBigExpImpl::StdDegRevLexImpl::myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const
  {
    CoCoA_ASSERT(0 <= GrDim && GrDim <= 1);
    if (GrDim == 0) return 0; // ie: wrt to 0-grading
    BigInt deg1;
    BigInt deg2;
    for (long i=0; i < myNumIndets; ++i)
    {
      deg1 += v1[i];
      deg2 += v2[i];
    }
    if (deg1 != deg2) return (deg1>deg2 ? 1 : -1 );
    return 0;
  }


//--- MatrixOrderingImpl

  PPMonoidEvBigExpImpl::MatrixOrderingImpl::MatrixOrderingImpl(long NumIndets, long GradingDim, const ConstMatrixView& OrderMatrix):
    CmpBase(NumIndets, GradingDim)
  {
    CoCoA_ASSERT(0 <= NumIndets);
    CoCoA_ASSERT(0 <= GradingDim && GradingDim <= NumIndets);
    CoCoA_ASSERT(NumRows(OrderMatrix) == NumIndets);
    CoCoA_ASSERT(NumCols(OrderMatrix) == NumIndets);

    myOrderMatrix.resize(NumIndets, vector<BigInt>(NumIndets));
    for (long i=0; i < NumIndets; ++i)
      for (long j=0; j < NumIndets; ++j)
        myOrderMatrix[i][j] = ConvertTo<BigInt>(OrderMatrix(i,j));
  }


  void PPMonoidEvBigExpImpl::MatrixOrderingImpl::myWDeg(degree& d, const BigInt* v) const
  {
    BigInt deg;
    for (long g=0; g < myGradingDim; ++g)
    {
      deg = 0;
      for (long i=0; i < myNumIndets; ++i)
        deg += v[i] * myOrderMatrix[g][i];
      SetComponent(d, g, deg);
    }
  }


  static int MatrixOrderingCmpArrays(const vector<vector<BigInt> >& M, const BigInt* v1, const BigInt* v2, long UpTo)
  {
    CoCoA_ASSERT(0 <= UpTo && UpTo <= len(M));
    BigInt s1, s2; // automatically set to 0
    const long ncols = len(M[0]);
    for (long r=0; r < UpTo; ++r)
    {
      for (long i=0; i < ncols; ++i)
      {
        s1 += v1[i]*M[r][i];
        s2 += v2[i]*M[r][i];
      }
      if (s1 != s2) return (s1>s2 ? 1 : -1 );
      s1 = 0;
      s2 = 0;
    }
    return 0;
  }


  int PPMonoidEvBigExpImpl::MatrixOrderingImpl::myCmpExpvs(const BigInt* v1, const BigInt* v2) const
  {
    return MatrixOrderingCmpArrays(myOrderMatrix, v1, v2, myNumIndets);
  }


  int PPMonoidEvBigExpImpl::MatrixOrderingImpl::myCmpWDegPartial(const BigInt* v1, const BigInt* v2, long GrDim) const
  {
    CoCoA_ASSERT(0 <= GrDim && GrDim <= myGradingDim);
    return MatrixOrderingCmpArrays(myOrderMatrix, v1, v2, GrDim);
  }


  //////////////////////////////////////////////////////////////////
  // Pseudo-ctors

  PPMonoid NewPPMonoidEv(const std::vector<symbol>& IndetNames, const PPOrdering& ord, PPExpSize ExpSize)
  {
    // Sanity check on the indet names given.
    const long nvars = NumIndets(ord);

    if (len(IndetNames) != nvars ||
        !AreDistinct(IndetNames) ||
        !AreArityConsistent(IndetNames))
      CoCoA_THROW_ERROR1(ERR::BadIndetNames);

    if (ExpSize == PPExpSize::small)
      return PPMonoid(new PPMonoidEvSmallExpImpl(IndetNames, ord));
    return PPMonoid(new PPMonoidEvBigExpImpl(IndetNames, ord));
  }

  PPMonoid NewPPMonoidEv(const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor, PPExpSize ExpSize)
  {
    return NewPPMonoidEv(IndetNames, OrdCtor(len(IndetNames)), ExpSize);
  }



} // end of namespace CoCoA
