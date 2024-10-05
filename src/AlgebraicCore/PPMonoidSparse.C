//   Copyright (c)  2005,2007,2010,2015,2021  John Abbott and Anna M. Bigatti

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


// Implementation of class PPMonoidSparseImpl

#include "CoCoA/PPMonoidSparse.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DivMask.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/assert.H"
#include "CoCoA/degree.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <algorithm>
using std::swap;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
#include <list>
using std::list;
#include <memory>
using std::unique_ptr;
#include <new>
//for placement new
#include <string>
using std::string;
#include <utility>
using std::pair;
#include <vector>
using std::vector;


namespace CoCoA
{

  // Oh great!  There's no doc.  2022-02-09
  // A PP in a PPMonoidSparse is represented as a std::list<IndexExp>
  // where the struct IndexExp is defined immediately below:
  //  it is just a pair of longs being indeterminate index, and exponent
  //  with the rule that the exponent must be positive (surely non-zero).
  // The entries in a list are ordered so that the indet indexes are
  // in (strictly) increasing order.
  // IDEA (NYI): include a first IndexExp with index -1 and exp being
  // the total degree (this would be helpful for deg-compat orderings).

  namespace // anonymous namespace for file local definitions
  {
    struct IndexExp
    {
    public:
      IndexExp(long IndetIndex, long exp): myIndetIndex(IndetIndex), myExp(exp) {}
    public:
      long myIndetIndex;
      long myExp;
    };


    std::ostream& operator<<(std::ostream& out, const IndexExp& ie)
    {
    if (!out) return out;  // short-cut for bad ostreams
      out << "IndexExp(" << ie.myIndetIndex << ", " << ie.myExp << ")";
      return out;
    }

  } // end of anon namespace


/////////////////////////////////////////////////////////////////////////////

  class PPMonoidSparseImpl: public PPMonoidBase
  {
  protected:
    // This pseudo-ctor is the only fn which calls the ctor:
    friend PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames,
                                      const PPOrdering& ord);
    PPMonoidSparseImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
    virtual ~PPMonoidSparseImpl() {};
  public: // disable copy ctor and assignment
    PPMonoidSparseImpl(const PPMonoidSparseImpl&) = delete;
    PPMonoidSparseImpl& operator=(const PPMonoidSparseImpl&) = delete;

  public:
    typedef PPMonoidElemRawPtr RawPtr;           ///< just to save typing
    typedef PPMonoidElemConstRawPtr ConstRawPtr; ///< just to save typing
  private:
    typedef std::list<IndexExp> value_t; // ***ACTUAL REPRESENTATION TYPE***
    static value_t& import(RawPtr rawpp);
    static const value_t& import(ConstRawPtr rawpp);

    // functions every PPMonoid must implement
///???    const PPOrdering& myOrdering() const  { return myPPO; }
    const std::vector<PPMonoidElem>& myIndets() const override;            ///< std::vector whose n-th entry is n-th indet as PPMonoidElem

    const PPMonoidElem& myOne() const override;
    PPMonoidElemRawPtr myNew() const override;                              ///< ctor from nothing
    PPMonoidElemRawPtr myNew(PPMonoidElemConstRawPtr rawpp) const override; ///< ctor from another pp
    PPMonoidElemRawPtr myNew(const std::vector<long>& expv) const override; ///< ctor from exp vector
//NYI    virtual PPMonoidElemRawPtr myNew(const std::vector<BigInt>& v) const;     ///< ctor from exp vector
    PPMonoidElemRawPtr myNewIndet(long k) const;  // just for creating the n-th indet
    void myDelete(RawPtr rawpp) const override;                                ///< dtor, frees pp

    void mySwap(RawPtr rawpp1, RawPtr rawpp2) const override;                  ///< swap(pp1, pp2)
    void myAssignOne(RawPtr rawpp) const override;                             ///< pp = 1
    void myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const override;            ///< pp = pp1
    void myAssign(RawPtr rawpp, const std::vector<long>& expv) const override; ///< pp = v (assign from exp vector)

    void myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;  ///< pp = pp1*pp2
    void myMulIndetPower(RawPtr rawpp, long indet, long exp) const override;          ///< pp *= indet^exp
    void myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;  ///< pp = pp1/pp2
    void myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;///< pp = pp1/gcd(pp1,pp2)
    void myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;  ///< pp = gcd(pp1,pp2)
    void myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;  ///< pp = lcm(pp1,pp2)
    void myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const override;                  ///< pp = radical(pp1)
    void myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const override;  // pp = pp1^exp (non-trivial), assumes exp >= 0
///NYI    virtual void myPowerBigExp(RawPtr rawpp, ConstRawPtr rawpp1, const BigInt& EXP) const;// pp = pp1^EXP (non-trivial), assumes EXP >= 0
    void myPowerOverflowCheck(ConstRawPtr rawpp1, long exp) const override;           ///< throw if pp1^exp would overflow, assumes exp >= 0

    bool myIsOne(ConstRawPtr rawpp) const override;                            ///< true iff pp = 1?
    bool myIsIndet(long& index, ConstRawPtr rawpp) const override;             ///< true iff pp is an indet
    bool myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;   ///< are pp1, pp2 coprime?
    bool myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;     ///< is pp1 equal to pp2?
    bool myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override; ///< is pp1 is divisible by pp2?
    bool myIsSqFree(ConstRawPtr rawpp) const override;                         ///< is pp equal to its radical?

    int myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;          ///< <0, =0, >0 as pp1 < = > pp2
    int myCmpOrdvs(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const;      ///< <0, =0, >0 as pp1 < = > pp2 (compare just first n compts)
    //    virtual int myHomogCmp(ConstRawPtr rawt1, ConstRawPtr rawt2) const;   ///< <0, =0, >0 as t1 < = > t2 assuming t1 and t2 have same multi-degree
    //    virtual int myHomogDegRevLex(ConstRawPtr rawt1, ConstRawPtr rawt2) const; ///< <0, =0, >0 as t1 < = > t2 ??? degrevlex assuming t1 and t2 have same multi-degree TO BE REMOVED

    long myStdDeg(ConstRawPtr rawpp) const override;                           ///< standard degree of pp
    void myWDeg(degree& d, ConstRawPtr rawpp) const override;                  ///< d = grading(pp)
    int myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const override;      ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2)
    int myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const override;      ///< <0, =0, >0 as wdeg(pp1) < = > wdeg(pp2) wrt the first n rows of weights
    long myExponent(ConstRawPtr rawpp, long indet) const override;             ///< degree of pp in indet
    void myBigExponent(BigInt& E, ConstRawPtr rawpp, long indet) const override; ///< degree of pp in indet

    void myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const override; ///< get exponents, SHOULD BE vector<BigInt> ????
    void myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const override; ///< get exponents, SHOULD BE vector<BigInt> ????
    void myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const override;    ///< v[i] = true if i-th indet has exponent != 0
    void myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const override; ///< computes the DivMask for pp according to DivMaskImpl
    void myOutputSelf(std::ostream& out) const override;                        ///< print value of PPMonoid
///USE DEFAULT IMPL IN PPMONOIDBASE    virtual void myOutput(std::ostream& out, ConstRawPtr rawpp) const;   ///< NOT PURE!!
///???    virtual void myOutput_OM(OpenMath::OutputChannel& OMOut, ConstRawPtr rawpp) const;///< NOT PURE!!

  protected: // Data members
    std::vector<PPMonoidElem> myIndetVector; ///< the indets as PPMonoidElems
    std::unique_ptr<PPMonoidElem> myOnePtr;
///???    PPOrdering myPPO;
  };


  inline PPMonoidSparseImpl::value_t& PPMonoidSparseImpl::import(RawPtr rawpp)
  {
    return *static_cast<value_t*>(rawpp.myRawPtr());
  }

  inline const PPMonoidSparseImpl::value_t& PPMonoidSparseImpl::import(ConstRawPtr rawpp)
  {
    return *static_cast<const value_t*>(rawpp.myRawPtr());
  }


  PPMonoidSparseImpl::PPMonoidSparseImpl(const std::vector<symbol>& IndetNames, const PPOrdering& ord):
      PPMonoidBase(ord, IndetNames)
///???      myPPO(ord)
  {
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myOnePtr.reset(new PPMonoidElem(PPMonoid(this)));
    myIndetVector.reserve(myNumIndets);
    // Fill vector of indets
    for (long i=0; i < myNumIndets; ++i)
      myIndetVector.push_back(PPMonoidElem(PPMonoid(this), myNewIndet(i)));

    myRefCountZero();  // ignore the "internal" references created above
  }


  const std::vector<PPMonoidElem>& PPMonoidSparseImpl::myIndets() const
  {
    return myIndetVector;
  }


  const PPMonoidElem& PPMonoidSparseImpl::myOne() const
  {
    return *myOnePtr;
  }


  PPMonoidElemRawPtr PPMonoidSparseImpl::myNew() const
  {
    unique_ptr<value_t> tmp;
    tmp.reset(new value_t());
    return PPMonoidElemRawPtr(tmp.release());
  }


  PPMonoidElemRawPtr PPMonoidSparseImpl::myNew(PPMonoidElemConstRawPtr rawpp) const
  {
    unique_ptr<value_t> tmp;
    tmp.reset(new value_t(import(rawpp)));
    return PPMonoidElemRawPtr(tmp.release());
  }


  PPMonoidElemRawPtr PPMonoidSparseImpl::myNew(const std::vector<long>& expv) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    unique_ptr<value_t> tmp;
    tmp.reset(new value_t());

    for (long i=0; i < myNumIndets; ++i)
      if (expv[i] != 0)
        tmp->push_back(IndexExp(i, expv[i]));

    return PPMonoidElemRawPtr(tmp.release());
  }

//NYI    virtual PPMonoidElemRawPtr myNew(const std::vector<BigInt>& v) const;     ///< ctor from exp vector

  PPMonoidElemRawPtr PPMonoidSparseImpl::myNewIndet(long k) const
  {
    CoCoA_ASSERT(0 <= k && k < myNumIndets);
    unique_ptr<value_t> tmp;
    tmp.reset(new value_t());

    tmp->push_back(IndexExp(k, 1));

    return PPMonoidElemRawPtr(tmp.release());
  }


  void PPMonoidSparseImpl::myDelete(RawPtr rawpp) const
  {
    unique_ptr<value_t> tmp;
    tmp.reset(&import(rawpp));
  }


  void PPMonoidSparseImpl::mySwap(RawPtr rawpp1, RawPtr rawpp2) const
  {
    std::swap(import(rawpp1), import(rawpp2));
  }


  void PPMonoidSparseImpl::myAssignOne(RawPtr rawpp) const
  {
    import(rawpp).clear();
  }

  
  void PPMonoidSparseImpl::myAssign(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    import(rawpp) = import(rawpp1);
  }

  
  void PPMonoidSparseImpl::myAssign(RawPtr rawpp, const std::vector<long>& expv) const
  {
    RawPtr rawrhs = myNew(expv);
    mySwap(rawpp, rawrhs);
    myDelete(rawrhs);
  }

  
  void PPMonoidSparseImpl::myMul(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator i = import(rawpp1).begin();
    const value_t::const_iterator endi = import(rawpp1).end();
    value_t::const_iterator j = import(rawpp2).begin();
    const value_t::const_iterator endj = import(rawpp2).end();
    while (i != endi && j != endj)
    {
      if (i->myIndetIndex < j->myIndetIndex)
      {
        ans->push_back(*i++);
        continue;
      }
      if (i->myIndetIndex > j->myIndetIndex)
      {
        ans->push_back(*j++);
        continue;
      }
      ans->push_back(IndexExp(i->myIndetIndex, i->myExp + j->myExp)); // assume combined exp is non-zero
      ++i;
      ++j;
    }
    // At most only one of the two lines below will actually do anything
    ans->insert(ans->end(), i, endi);
    ans->insert(ans->end(), j, endj);

    std::swap(import(rawpp), *ans);
  }

  void PPMonoidSparseImpl::myMulIndetPower(RawPtr rawpp, long indet, long exp) const
  {
    CoCoA_ASSERT(exp > 0);
    value_t::iterator it = import(rawpp).begin();
    const value_t::iterator endit = import(rawpp).end();
    while (it != endit && it->myIndetIndex < indet) ++it;
    if (it != endit && it->myIndetIndex == indet)
    {
      it->myExp += exp; // assume the sum is non-zero
      return;
    }
    import(rawpp).insert(it, IndexExp(indet, exp));
  }


  void PPMonoidSparseImpl::myDiv(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator it1 = import(rawpp1).begin();
    const value_t::const_iterator end1 = import(rawpp1).end();
    value_t::const_iterator it2 = import(rawpp2).begin();
    const value_t::const_iterator end2 = import(rawpp2).end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex < it2->myIndetIndex)
      {
        ans->push_back(*it1);
        ++it1;
        continue;
      }
      CoCoA_ASSERT("Exponent Underflow" && it1->myIndetIndex == it2->myIndetIndex);
      CoCoA_ASSERT("Exponent Underflow" && it1->myExp >= it2->myExp);
      if (it1->myExp > it2->myExp)
        ans->push_back(IndexExp(it1->myIndetIndex, it1->myExp - it2->myExp));
      ++it1;
      ++it2;
    }
    CoCoA_ASSERT(it2 == end2);
    ans->insert(ans->end(), it1, end1);

    std::swap(import(rawpp), *ans);
  }


  void PPMonoidSparseImpl::myColon(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator it1 = import(rawpp1).begin();
    const value_t::const_iterator end1 = import(rawpp1).end();
    value_t::const_iterator it2 = import(rawpp2).begin();
    const value_t::const_iterator end2 = import(rawpp2).end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex < it2->myIndetIndex)
      {
        ans->push_back(*it1);
        ++it1;
        continue;
      }
      if (it1->myIndetIndex > it2->myIndetIndex) { ++it2; continue; }
      if (it1->myExp > it2->myExp)
        ans->push_back(IndexExp(it1->myIndetIndex, it1->myExp - it2->myExp));
      ++it1;
      ++it2;
    }
    ans->insert(ans->end(), it1, end1);
    // Ignore any remaining factors in it2

    std::swap(import(rawpp), *ans);
  }


  void PPMonoidSparseImpl::myGcd(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator it1 = import(rawpp1).begin();
    const value_t::const_iterator end1 = import(rawpp1).end();
    value_t::const_iterator it2 = import(rawpp2).begin();
    const value_t::const_iterator end2 = import(rawpp2).end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex < it2->myIndetIndex) { ++it1; continue; }
      if (it1->myIndetIndex > it2->myIndetIndex) { ++it2; continue; }
      ans->push_back(IndexExp(it1->myIndetIndex, std::min(it1->myExp, it2->myExp)));
      ++it1;
      ++it2;
    }
    // Ignore any remaining factors in it1 or it2

    std::swap(import(rawpp), *ans);
  }

  
  void PPMonoidSparseImpl::myLcm(RawPtr rawpp, ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator it1 = import(rawpp1).begin();
    const value_t::const_iterator end1 = import(rawpp1).end();
    value_t::const_iterator it2 = import(rawpp2).begin();
    const value_t::const_iterator end2 = import(rawpp2).end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex < it2->myIndetIndex)
      {
        ans->push_back(*it1);
        ++it1;
        continue;
      }
      if (it1->myIndetIndex > it2->myIndetIndex)
      {
        ans->push_back(*it2);
        ++it2;
        continue;
      }
      ans->push_back(IndexExp(it1->myIndetIndex, std::max(it1->myExp, it2->myExp)));
      ++it1;
      ++it2;
    }
    // At most only one of the two lines below will actually do anything
    ans->insert(ans->end(), it1, end1);
    ans->insert(ans->end(), it2, end2);

    std::swap(import(rawpp), *ans);
  }


  void PPMonoidSparseImpl::myRadical(RawPtr rawpp, ConstRawPtr rawpp1) const
  {
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator it = import(rawpp1).begin();
    const value_t::const_iterator endit = import(rawpp1).end();
    while (it != endit)
    {
      ans->push_back(IndexExp(it->myIndetIndex, 1));
      ++it;
    }

    std::swap(import(rawpp), *ans);
  }


  void PPMonoidSparseImpl::myPowerSmallExp(RawPtr rawpp, ConstRawPtr rawpp1, long exp) const
  {
    CoCoA_ASSERT(exp >= 0);
#ifdef CoCoA_DEBUG
    myPowerOverflowCheck(rawpp1, exp);
#endif
    unique_ptr<value_t> ans;
    ans.reset(new value_t);
    value_t::const_iterator it = import(rawpp1).begin();
    const value_t::const_iterator endit = import(rawpp1).end();
    while (it != endit)
    {
      ans->push_back(IndexExp(it->myIndetIndex, exp*it->myExp));
      ++it;
    }

    std::swap(import(rawpp), *ans);
  }


  void PPMonoidSparseImpl::myPowerOverflowCheck(ConstRawPtr rawpp, long exp) const
  {
    CoCoA_ASSERT(exp >= 0);
    // nothing to check as exps are unlimited

    if (exp == 0 || exp == 1) return;
    // *** BUG *** in the next line "long" should be the type used inside "IndexExp"
    const long limit = numeric_limits<long>::max()/exp;
    value_t::const_iterator it = import(rawpp).begin();
    const value_t::const_iterator endit = import(rawpp).end();
    while (it != endit)
    {
      if (it->myExp > limit)
        CoCoA_THROW_ERROR(ERR::ExpTooBig, "PPMonoidSparseImpl::myPowerOverflowCheck");
      ++it;
    }
  }


  bool PPMonoidSparseImpl::myIsOne(ConstRawPtr rawpp) const
  {
    return import(rawpp).empty();
  }


  bool PPMonoidSparseImpl::myIsIndet(long& index, ConstRawPtr rawpp) const
  {
    value_t::const_iterator it = import(rawpp).begin();
    const value_t::const_iterator endit = import(rawpp).end();
    if (it == endit) return false;
    if (it->myExp != 1) return false;
    const long IndetIndex = it->myIndetIndex;
    ++it;
    if (it != endit) return false;
    index = IndetIndex;
    return true;
  }


  bool PPMonoidSparseImpl::myIsCoprime(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator it1 = T1.begin();
    const value_t::const_iterator end1 = T1.end();
    value_t::const_iterator it2 = T2.begin();
    const value_t::const_iterator end2 = T2.end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex == it2->myIndetIndex) return false;
      if (it1->myIndetIndex < it2->myIndetIndex) ++it1;
      else ++it2;
    }
    return true;
  }


bool PPMonoidSparseImpl::myIsEqual(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator it1 = T1.begin();
    const value_t::const_iterator end1 = T1.end();
    value_t::const_iterator it2 = T2.begin();
    const value_t::const_iterator end2 = T2.end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex != it2->myIndetIndex || it1->myExp != it2->myExp) return false;
      ++it1;
      ++it2;
    }
    return (it1 == end1) && (it2 == end2);
  }


  bool PPMonoidSparseImpl::myIsDivisible(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    value_t::const_iterator it1 = T1.begin();
    const value_t::const_iterator end1 = T1.end();
    value_t::const_iterator it2 = T2.begin();
    const value_t::const_iterator end2 = T2.end();
    while (it1 != end1 && it2 != end2)
    {
      if (it1->myIndetIndex > it2->myIndetIndex) return false;
      if (it1->myIndetIndex < it2->myIndetIndex) { ++it1; continue; }
      if (it1->myExp < it2->myExp) return false;
      ++it1;
      ++it2;
    }
    return (it2 == end2);
  }


  bool PPMonoidSparseImpl::myIsSqFree(ConstRawPtr rawpp) const
  {
    const value_t& T = import(rawpp);
    value_t::const_iterator it = T.begin();
    const value_t::const_iterator endit = T.end();
    while (it != endit)
    {
      if (it->myExp > 1) return false;
      ++it;
    }
    return true;
  }


  int PPMonoidSparseImpl::myCmp(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myCmpOrdvs(rawpp1, rawpp2, myNumIndets);
  }

    //    virtual int myHomogCmp(ConstRawPtr rawt1, ConstRawPtr rawt2) const;   ///< <0, =0, >0 as t1 < = > t2 assuming t1 and t2 have same multi-degree
    //    virtual int myHomogDegRevLex(ConstRawPtr rawt1, ConstRawPtr rawt2) const; ///< <0, =0, >0 as t1 < = > t2 ??? degrevlex assuming t1 and t2 have same multi-degree TO BE REMOVED

  long PPMonoidSparseImpl::myStdDeg(ConstRawPtr rawpp) const
  {
    const value_t& T = import(rawpp);
    value_t::const_iterator it = T.begin();
    const value_t::const_iterator endit = T.end();
    long deg = 0;
    while (it != endit)
    {
      deg += it->myExp;  //??? BUG check for overflow????
      ++it;
    }
    return deg;
  }


  int PPMonoidSparseImpl::myCmpOrdvs(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const
  {
    // This fn computes the ordvs of pp1 & pp2 one compt at a time, and compares them lexicographically.
    if (n < 0 || n > NumIndets(myOrdering())) CoCoA_THROW_ERROR(ERR::BadIndex, "myCmpOrdvs");
    const ConstMatrixView& M = OrdMat(myOrdering());
    const value_t& T1 = import(rawpp1);
    const value_t& T2 = import(rawpp2);
    typedef value_t::const_iterator iter;
    const iter end1 = T1.end();
    const iter end2 = T2.end();
    for (int i=0; i < n; ++i)
    {
      BigInt comp1;
      for (iter it1 = T1.begin(); it1 != end1; ++it1)
      {
        const BigInt Mij = ConvertTo<BigInt>(M(i, it1->myIndetIndex));
        comp1 += it1->myExp * Mij;
      }
      BigInt comp2;
      for (iter it2 = T2.begin(); it2 != end2; ++it2)
      {
        const BigInt Mij = ConvertTo<BigInt>(M(i, it2->myIndetIndex));
        comp2 += it2->myExp * Mij;
      }
      if (comp1 == comp2) continue;
      return sign(comp1-comp2);
    }
    return 0;
  }


  void PPMonoidSparseImpl::myWDeg(degree& d, ConstRawPtr rawpp) const
  {
    const int dim = GradingDim(d);
    if (dim != GradingDim(myOrdering())) CoCoA_THROW_ERROR(ERR::MixedDegrees, "PPMonoidSparseImpl::myWDeg");
    const ConstMatrixView& M = OrdMat(myOrdering());
    const value_t& T = import(rawpp);
    typedef value_t::const_iterator iter;
    const iter endit = T.end();
    for (int i=0; i < dim; ++i)
    {
      BigInt comp;
      for (iter it = T.begin(); it != endit; ++it)
      {
        BigInt Mij = ConvertTo<BigInt>(M(i, it->myIndetIndex));
        comp += it->myExp * Mij;
      }
      SetComponent(d, i, comp);
    }
  }


  int PPMonoidSparseImpl::myCmpWDeg(ConstRawPtr rawpp1, ConstRawPtr rawpp2) const
  {
    return myCmpOrdvs(rawpp1, rawpp2, GradingDim(myOrdering()));
  }


  int PPMonoidSparseImpl::myCmpWDegPartial(ConstRawPtr rawpp1, ConstRawPtr rawpp2, long n) const
  {
    if (n < 0 || n > GradingDim(myOrdering())) CoCoA_THROW_ERROR(ERR::BadIndex, "myCmpWDegPartial");
    return myCmpOrdvs(rawpp1, rawpp2, n);
  }


  long PPMonoidSparseImpl::myExponent(ConstRawPtr rawpp, long var) const
  {
    CoCoA_ASSERT(0 <= var && var <= myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator it = T.begin();
    const value_t::const_iterator endit = T.end();
    while (it != endit)
    {
      if (it->myIndetIndex == var) return it->myExp;
      if (it->myIndetIndex > var) return 0;
      ++it;
    }
    return 0;
  }

  void PPMonoidSparseImpl::myBigExponent(BigInt& E, ConstRawPtr rawpp, long indet) const
  {
    E = myExponent(rawpp, indet);
  }

  
  void PPMonoidSparseImpl::myExponents(std::vector<long>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator it = T.begin();
    const value_t::const_iterator endit = T.end();
    for (long j=0; j < myNumIndets; ++j) expv[j] = 0;
    while (it != endit)
    {
      expv[it->myIndetIndex] = it->myExp;
      ++it;
    }
  }

  void PPMonoidSparseImpl::myBigExponents(std::vector<BigInt>& expv, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(expv) == myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator it = T.begin();
    const value_t::const_iterator endit = T.end();
    for (long j=0; j < myNumIndets; ++j) expv[j] = 0;
    while (it != endit)
    {
      expv[it->myIndetIndex] = it->myExp;
      ++it;
    }
  }

  
  void PPMonoidSparseImpl::myIndetsIn(std::vector<bool>& v, ConstRawPtr rawpp) const
  {
    CoCoA_ASSERT(len(v) == myNumIndets);
    const value_t& T = import(rawpp);
    value_t::const_iterator it = T.begin();
    const value_t::const_iterator endit = T.end();
    while (it != endit)
    {
      v[it->myIndetIndex] = true; // ASSUMES it->myExp != 0
      ++it;
    }
  }


  void PPMonoidSparseImpl::myComputeDivMask(DivMask& dm, const DivMaskRule& DivMaskImpl, ConstRawPtr rawpp) const
  {
    // !!!!!SLUG!!!!! appallingly inefficient :-(
//    const value_t& T = import(rawpp);
    vector<SmallExponent_t> expv(myNumIndets);
//    const value_t::const_iterator endit = T.end();
//    for (value_t::const_iterator it = T.begin(); it != endit; ++it)
    for (const auto& IndetPwr: import(rawpp))
      expv[IndetPwr.myIndetIndex] = IndetPwr.myExp;
//      expv[it->myIndetIndex] = it->myExp;
    DivMaskImpl->myAssignFromExpv(dm, &expv[0], myNumIndets);
  }

  
  void PPMonoidSparseImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "PPMonoidSparse(" << myNumIndets << ", " << myOrd << ")";
  }

/// USE GENERIC IMPL in PPMonoid.C  std::ostream& PPMonoidSparse::myOutput(std::ostream& out, ConstRawPtr rawpp) const;   ///< NOT PURE!!
//???    virtual OpenMath::OutputChannel& myOutput_OM(OpenMath::OutputChannel& OMOut, ConstRawPtr rawpp) const;///< NOT PURE!!



  //----------------------------------------------------------------------
  // Pseudo-ctor for sparse PPMonoids

  PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrdering& ord)
  {
  // Sanity check on the indet names given.
    const long nvars = NumIndets(ord);
  
    if (len(IndetNames) != nvars)
      CoCoA_THROW_ERROR(ERR::BadNumIndets, "NewPPMonoidSparse(IndetNames,ord)");
    if (!AreDistinct(IndetNames))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidSparse(IndetNames,ord)");
    if (!AreArityConsistent(IndetNames))
      CoCoA_THROW_ERROR(ERR::BadIndetNames, "NewPPMonoidSparse(IndetNames,ord)");
  
    return PPMonoid(new PPMonoidSparseImpl(IndetNames, ord));
  }

  PPMonoid NewPPMonoidSparse(const std::vector<symbol>& IndetNames, const PPOrderingCtor& OrdCtor)
  {
    return NewPPMonoidSparse(IndetNames, OrdCtor(len(IndetNames)));
  }



}// end of namespace CoCoA




// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/PPMonoidSparse.C,v 1.41 2024/07/03 09:35:15 bigatti Exp $
// $Log: PPMonoidSparse.C,v $
// Revision 1.41  2024/07/03 09:35:15  bigatti
// Summary: error codes (redmine 92):
// removed BadIntedIndex, BadDegIndex, BadComptIndex
//
// Revision 1.40  2023/03/13 19:55:05  abbott
// Summary: Switched to "new" iterator for loop
//
// Revision 1.39  2022/02/18 14:11:56  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.38  2022/02/11 09:49:35  abbott
// Summary: Updated copyright notices (redmine 855)
//
// Revision 1.37  2022/02/08 20:18:54  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.36  2021/10/30 13:47:39  abbott
// Summary: Used keyword override where appropriate (redmine 1625)
//
// Revision 1.35  2021/01/07 15:16:52  abbott
// Summary: Corrected copyright
//
// Revision 1.34  2020/06/17 15:49:25  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.33  2020/02/11 16:56:41  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.32  2020/02/11 16:12:18  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.31  2019/03/04 10:33:12  abbott
// Summary: Changed auto_ptr into unqiue_ptr
//
// Revision 1.30  2018/05/22 14:16:40  abbott
// Summary: Split BigRat into BigRat (class defn + ctors) and BigRatOps
//
// Revision 1.29  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.28  2018/05/17 15:38:23  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.27  2017/04/18 12:50:07  abbott
// Summary: Corrected ifdef use of CoCoA_THREADSAFE_HACK and CoCoA_DEBUG
//
// Revision 1.26  2017/02/22 12:32:45  abbott
// Summary: Changed exponents to long (from int)
//
// Revision 1.25  2016/11/11 14:15:33  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.24  2016/11/03 12:25:25  abbott
// Summary: Changed IsRadical (for PPMonoidElem) into IsSqFree
//
// Revision 1.23  2015/12/01 13:11:01  abbott
// Summary: Changed mem fn PPOrderingCtor::myCtor into operator(); also for ModuleOrderingCtor; see issue 829
//
// Revision 1.22  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.21  2015/06/30 12:55:21  abbott
// Summary: Added new fn myIndetsIn
// Author: JAA
//
// Revision 1.20  2015/04/16 20:20:58  abbott
// Summary: Silly change to avoid a compiler warning about unused param
// Author: JAA
//
// Revision 1.19  2015/04/16 16:36:33  abbott
// Summary: Cleaned impls of myPowerOverflowCheck
// Author: JAA
//
// Revision 1.18  2015/04/13 14:42:08  abbott
// Summary: Added myPowerOverflowCheck (1st version)
// Author: JAA
//
// Revision 1.17  2014/07/31 14:45:17  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.16  2014/07/03 15:36:35  abbott
// Summary: Cleaned up impl of PPMonoids: moved myIndetSymbols & myNumIndets to base class
// Author: JAA
//
// Revision 1.15  2014/05/14 15:57:15  bigatti
// -- added "using" for clang with superpedantic flag
//
// Revision 1.14  2014/04/30 16:10:31  abbott
// Summary: Removed pointless include
// Author: JAA
//
// Revision 1.13  2012/04/02 15:03:02  bigatti
// -- using insert instead of copy+back_inserter
//
// Revision 1.12  2012/02/08 16:13:18  bigatti
// -- using  insert  instead of  back_inserter
//
// Revision 1.11  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.10  2011/06/23 16:07:13  abbott
// Added incomplete but compiling PPMonoidSparse: first prototype,
// simple rather than fast!
//
// Revision 1.9  2010/11/30 11:18:11  bigatti
// -- renamed IndetName --> IndetSymbol
//
// Revision 1.8  2010/10/06 14:10:24  abbott
// Added increments to the ref count in ring and PPMonoid ctors to make
// them exception safe.
//
// Revision 1.7  2008/06/05 14:57:58  bigatti
// -- added: include for MemPool and DivMask
//
// Revision 1.6  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.5  2007/12/04 14:27:06  bigatti
// -- changed "log(pp, i)" into "exponent(pp, i)"
//
// Revision 1.4  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/05/31 14:54:31  bigatti
// -- now using AreDistinct and AreArityConsistent for sanity check on
//    indet names
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/03/08 17:43:11  cocoa
// Swapped order of args to the NewPPMonoid pseudo ctors.
//
// Revision 1.8  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.7  2007/01/18 14:22:44  cocoa
// -- added "raw" to RawPtr arguments
//
// Revision 1.6  2007/01/15 13:39:54  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.5  2006/11/24 17:04:32  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/10/19 13:56:24  cocoa
// Added #include<new> whenever needed (i.e. for placement new).
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.5  2006/04/28 16:33:51  cocoa
// Used SmartPtrIRC for PPOrderings.
//
// Revision 1.4  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.3  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.2  2005/10/24 15:22:12  cocoa
// Fixed a couple of buglets
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1  2005/09/22 18:04:17  cocoa
// It compiles; the tests run OK.  The examples compile.
// No documentation -- the mindless eurocrats have rendered
// me mindless too.
//
