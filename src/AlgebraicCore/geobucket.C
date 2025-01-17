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


#include "CoCoA/geobucket.H"

#include "CoCoA/assert.H"
#include "CoCoA/utils.H"

#include <iostream>
using std::ostream;
using std::endl;

//???#include <algorithm>


namespace CoCoA
{

  //  static const long gbk_minlen = 4*4*4;   // suggested
  //  static const long gbk_factor = 4;   // suggested
  static const long gbk_minlen = 128;
  static const long gbk_factor = 4;
  static const long gbk_numbuckets = 20;  // guarantees no realloc for < 2^(5+2*20) terms

  //----------------------------------------------------------------------//
  // inline functions
  //----------------------------------------------------------------------//

  //----------  bucket functions  ----------//


  geobucket::bucket::bucket(const SparsePolyRing& P, long MaxLen):
      myCoeff(CoeffRing(P),1), myPoly(P)
  {
    myMaxLen = MaxLen;
    myApproxLen = 0;
    //  IamNormalized = true;
  }


  geobucket::bucket::bucket(const bucket& b):
      myCoeff(b.myCoeff), myPoly(b.myPoly)
  {
    myMaxLen = b.myMaxLen;
    myApproxLen = b.myApproxLen;
    //  IamNormalized = b.IamNormalized;
  }


  geobucket::bucket::~bucket()
  {
  }


  inline void geobucket::bucket::myNormalize(void)
  {
    if (IsOne(myCoeff)) return;
    SparsePolyRingPtr(owner(myPoly))->myMulByCoeff(raw(myPoly), raw(myCoeff));
    myCoeff = 1;
  }


  void geobucket::bucket::myAddClear(RingElem& f, long FLen)
  {
    myNormalize();
    SparsePolyRingPtr(owner(myPoly))->myAddClear(raw(myPoly), raw(f));
    myApproxLen += FLen;
  }


  void geobucket::bucket::myAddClear(bucket& b)
  {
    b.myNormalize();
    myAddClear(b.myPoly, b.myApproxLen);
    b.myApproxLen = 0;
  }


  inline bool geobucket::bucket::myIsZeroAddLCs(const SparsePolyRing& P, geobucket::bucket& b1, geobucket::bucket& b2)
  {
    CoCoA_ASSERT(owner(b1.myPoly) == P);
    CoCoA_ASSERT(owner(b2.myPoly) == P);
    b1.myNormalize(); //  I think I may assume b1 is normalized (???)
    b2.myNormalize();
    --b2.myApproxLen;
    if ( !P->myIsZeroAddLCs(raw(b1.myPoly), raw(b2.myPoly))) return false;
    --b1.myApproxLen;
    return true;
  }


  inline void MoveLMToFront(geobucket::bucket& b1, geobucket::bucket& b2)
  {
    b1.myNormalize();
    b2.myNormalize();
    ++b1.myApproxLen;
    --b2.myApproxLen;
    SparsePolyRingPtr(owner(b1.myPoly))->myMoveLMToFront(raw(b1.myPoly), raw(b2.myPoly));
  }


  inline void MoveLMToFront(const SparsePolyRing& P, geobucket::bucket& b1, geobucket::bucket& b2)
  {
    b1.myNormalize();
    b2.myNormalize();
    ++b1.myApproxLen;
    --b2.myApproxLen;
    P->myMoveLMToFront(raw(b1.myPoly), raw(b2.myPoly));
  }


  void geobucket::bucket::myMulByCoeff(ConstRefRingElem coeff)
  {
    myCoeff *= coeff;
    //  IamNormalized = (IsOne(myCoeff));
  }


  void geobucket::bucket::myDivByCoeff(ConstRefRingElem coeff)
  {
    myNormalize();
    SparsePolyRingPtr(owner(myPoly))->myDivByCoeff(raw(myPoly), raw(coeff));
  }


  RingElem content(const geobucket::bucket& b)
  {
    CoCoA_ASSERT(!IsZero(b));
    RingElem c(content(b.myPoly));
    c *= b.myCoeff;
    return c;
  }


  ConstRefRingElem poly(geobucket::bucket& b)
  {
    b.myNormalize();
    return b.myPoly;
  }

//----------  geobucket constructors & destructors  ----------//

  geobucket::geobucket(const SparsePolyRing& P):
      myPolyRing(P)
  {
    IhaveLM = true;
    myBuckets.reserve(gbk_numbuckets);
    myPushBackZeroBucket(gbk_minlen);
  }


  geobucket::~geobucket()
  {
  }

//---------- friends ----------//

  ostream& operator<<(ostream& out, const geobucket& g)
  {
    if (!out) return out;  // short-cut for bad ostreams

    for (long i=0; i<len(g) ; ++i)
      out << "    /*myBuckets[" << i << "]*/ +"
          << g.myBuckets[i].myCoeff << " * ("
          << g.myBuckets[i].myPoly  << ")"<< endl;
    out << endl;
    return out;
  }


  void PrintLengths(std::ostream& out, const geobucket& g)
  {
    if (!out) return;  // short-cut for bad ostreams

    for (long i=0; i<len(g) ; ++i)
      out << "    -- len(myBuckets[" << i << "])"
          << " MaxLen = "    << g.myBuckets[i].myMaxLen
          << " ApproxLen = " << g.myBuckets[i].myApproxLen
          << " Len = "       << NumTerms(g.myBuckets[i].myPoly)
          << endl;
    out << endl;
  }


  void AddClear(RingElem& f, geobucket& gbk)
  {
    for ( long i=0 ; i<len(gbk) ; ++i )
    {
      gbk.myBuckets[i].myNormalize();
      SparsePolyRingPtr(owner(f))->myAddClear(raw(f), raw(gbk.myBuckets[i].myPoly));
      gbk.myBuckets[i].myApproxLen = 0;
    }
    gbk.IhaveLM = true;
  }


  long len(const geobucket& g)
  { return g.myLen(); }
  

  RingElem content(const geobucket& f)
  {
    CoCoA_ASSERT(!IsZero(f));
    RingElem cnt(CoeffRing(f));

    for ( long i=0 ; i < len(f) ; ++i )
      if (!IsZero(f.myBuckets[i]))
        if (IsOne(cnt = gcd(cnt, content(f.myBuckets[i])))) return cnt;
    return cnt;
  }


  void RemoveBigContent(geobucket& gbk)
  {
    RingElem cnt(content(gbk));
//  if ( IsBig(cnt) )
    gbk.myDivByCoeff(cnt);
  }


  //----------  geobucket functions  ----------//


  long geobucket::myLen() const
  { return len(myBuckets); }


  void geobucket::myAddClear(RingElem& f, long len)
  {
    IhaveLM = false;
    myBuckets[myBucketIndex(len)].myAddClear(f, len);
  }


  void geobucket::myDeleteLM(void)
  {
    CoCoA_ASSERT(IhaveLM);
    IhaveLM = false;
    --myBuckets[0].myApproxLen;
    myPolyRing->myDeleteLM(raw(myBuckets[0].myPoly));
  }


  void geobucket::mySetLM() const
  {
    myBuckets[0].myNormalize();
    if (myLen()==1)
    {
      IhaveLM = true;
      return;
    }
    int CMP;
    long gbk_len=myLen(), i;
    const SparsePolyRing& P = myPolyRing;
    bucket& b0 = myBuckets[0];

    while (!IhaveLM)
    {
      i=1;
      IhaveLM = true;
      if ( IsZero(b0) )
      {
        b0.myApproxLen = 0;
        while (i<gbk_len && IsZero(myBuckets[i]))
          ++i;
        if (i < gbk_len)
        {
          if (myBuckets[i].myApproxLen < gbk_minlen)
            b0.myAddClear(myBuckets[i]);
          else
            MoveLMToFront(b0, myBuckets[i]);
        }
        ++i;
      }
      for ( ; i < gbk_len; ++i)
        if ( !IsZero(myBuckets[i]) )
        {
          CMP = P->myCmpLPP(raw(b0.myPoly), raw(myBuckets[i].myPoly));
          if (CMP < 0)
          {
            if (myBuckets[i].myApproxLen < gbk_minlen)
              b0.myAddClear(myBuckets[i]);
            else
              MoveLMToFront(b0, myBuckets[i]);
          }
          if (CMP == 0)
            if (bucket::myIsZeroAddLCs(P, b0, myBuckets[i]))
            {
              IhaveLM = false;
              if (myBuckets[i].myApproxLen < gbk_minlen)
                b0.myAddClear(myBuckets[i]);
              break;
            }
        }
    }
  }


  void geobucket::myDivByCoeff(ConstRefRingElem coeff)  //to be tested: ANNA ???
  {
    if (!IsOne(coeff))
      for ( long l = myLen() ; l-- > 0 ; )  myBuckets[l].myDivByCoeff(coeff);
  }


  void geobucket::myMulByCoeff(ConstRefRingElem coeff)
  {
    if (!IsOne(coeff))
      for ( long l = myLen() ; l-- > 0 ; ) myBuckets[l].myMulByCoeff(coeff);
  }


  void geobucket::myPushBackZeroBucket(long MaxLen)
  {
    bucket b(myPolyRing, MaxLen);
    b.myApproxLen = 0;
    myBuckets.push_back(b);
  }


  void geobucket::myCascadeFrom(long i)
  {
    long MaxLen;

    while ( myBuckets[i].myApproxLen > (MaxLen=myBuckets[i].myMaxLen)*3 )
    {
      myBuckets[i].myApproxLen = NumTerms(myBuckets[i].myPoly);
      if ( myBuckets[i].myApproxLen <= MaxLen*2 ) break;
      IhaveLM = false;
      ++i;
      if (i==myLen()) myPushBackZeroBucket(MaxLen*gbk_factor);
      myBuckets[i].myAddClear(myBuckets[i-1]);
    }
  }


  long geobucket::myBucketIndex(long len)
  {
    long i=0;
    long bkt_len=gbk_minlen;
    while ( len>bkt_len )
    {
      ++i;
      bkt_len *= gbk_factor;
      if (i==myLen()) myPushBackZeroBucket(bkt_len);
    }
    return i;
  }


  void geobucket::myAddMulLM(ConstRefRingElem monom,
                           ConstRefRingElem g,
                           long gLen)
  { myAddMulLM(monom, g, gLen, SparsePolyRingBase::DontSkipLMg); }


  void geobucket::myAddMulLM(ConstRefRingElem monom,
                           ConstRefRingElem g,
                           long gLen,
                           SparsePolyRingBase::SkipLMFlag skip)
  {
    IhaveLM = false;
    if (gLen==0) gLen = NumTerms(g);
    long  i = myBucketIndex(gLen);
    if (skip == SparsePolyRingBase::SkipLMg) --gLen;
    myBuckets[i].myApproxLen += gLen;
    myBuckets[i].myNormalize();
    myPolyRing->myAddMulLM(raw(myBuckets[i].myPoly), raw(monom), raw(g), skip);
    myCascadeFrom(i);
  }


  void MoveLMToFront(RingElem& f, geobucket& gbk)
  {
    if (!gbk.IhaveLM) gbk.mySetLM();
    gbk.myPolyRing->myMoveLMToFront(raw(f),raw(gbk.myBuckets[0].myPoly));
    //  --gbk.myBuckets[0].myApproxLen;
    gbk.IhaveLM = false;
  }


  void MoveLMToBack(RingElem& f, geobucket& gbk)
  {
    if (!gbk.IhaveLM) gbk.mySetLM();
    gbk.myPolyRing->myMoveLMToBack(raw(f),raw(gbk.myBuckets[0].myPoly));
    //  --gbk.myBuckets[0].myApproxLen;
    gbk.IhaveLM = false;
  }


  RingElemAlias LC(const geobucket& gbk)
  {
    if (!gbk.IhaveLM) gbk.mySetLM();
    return LC(poly(gbk.myBuckets[0]));
  }


  void ReductionStep(geobucket& gbk, ConstRefRingElem g, long RedLen)
  {
    CoCoA_ASSERT(!IsZero(g));
    CoCoA_ASSERT(gbk.IhaveLM);
    CoCoA_ASSERT(!IsZero(poly(gbk.myBuckets[0])));

    const SparsePolyRing& P = gbk.myPolyRing;
    RingElem tmp_poly(P);

    P->myDivLM(raw(tmp_poly), raw(poly(gbk.myBuckets[0])), raw(g));
    P->myNegate(raw(tmp_poly), raw(tmp_poly));
    gbk.myDeleteLM();
    gbk.myAddMulLM(tmp_poly, g, RedLen, SparsePolyRingBase::SkipLMg);
  }


  void ReductionStepGCD(geobucket& gbk, ConstRefRingElem g, RingElem& gbkScale, long RedLen)
  {
    CoCoA_ASSERT(!IsZero(g));
    CoCoA_ASSERT(gbk.IhaveLM);
    CoCoA_ASSERT(!IsZero(gbk.myBuckets[0]));

    const SparsePolyRing& P = gbk.myPolyRing;
    RingElem tmp_poly(P);
    RingElem gScale(CoeffRing(P));
    {
      RingElem junk(CoeffRing(P));
      ///JAA  gbkScale = 1;
      ///JAA  ComputeFScaleAndGScale(CoeffRing(P), RawLC(gbk), P->myRawLC(raw(g)), raw(gbkScale), raw(gScale));
      GcdQuot(junk, gScale, gbkScale, LC(gbk), LC(g));
      if (IsMinusOne(gbkScale)) gbkScale = -gbkScale; else gScale = -gScale;
      if (IsInvertible(gbkScale))
      {
        gScale /= gbkScale;
        gbkScale = 1;
      }
    }
    gbk.myMulByCoeff(gbkScale);
    P->myDivLM(raw(tmp_poly), raw(poly(gbk.myBuckets[0])), raw(g));
    P->myNegate(raw(tmp_poly), raw(tmp_poly));
    gbk.myDeleteLM();
    gbk.myAddMulLM(tmp_poly, g, RedLen, SparsePolyRingBase::SkipLMg);
  }



} // end of namespace CoCoA


//-------------------------------------------------------
// Prototype template class for geobuckets
// Orig Author: Vibhash Singh (2022); copyright handed over to J Abbott & A M Bigatti
// 2023-12-10 Not yet fully tested.

// namespace CoCoA
// {
//   template<class T, size_t (*bucketIndex)(const T&), T (*binaryFunc)(const T&, const T&)> 
//   class GeoBucket 
//   {
//   public:
//     GeoBucket(const T &identity) : identity(identity) {};  
//     T myTotal() const { return myNormalize(); }
    
//     void compute(const T& N)
//     {
//       const size_t i = bucketIndex(N);
//       if (buckets.size() <= i)  { buckets.resize(i + 1, identity); }
//       buckets[i] = binaryFunc(buckets[i], N);
          
//       // cascade if necessary
//       size_t prevIndex = i;
//       while (true)
//       {
//         const size_t newIndex = bucketIndex(buckets[prevIndex]);

//         if (newIndex == prevIndex)
//           return;

//         if (buckets.size() <= newIndex) { buckets.resize(newIndex + 1, identity); }
//         buckets[newIndex] = binaryFunc(buckets[newIndex], buckets[prevIndex]);
//         buckets[prevIndex] = identity;
//         prevIndex = newIndex;
//       }
//     }

//   private:  
//     mutable vector<T> buckets;
//     T identity;

//     T myNormalize() const
//     {
//       T total = identity;
//       for (T& q: buckets)
//       {
//         total = binaryFunc(total, q);

//         // reset bucket initial value to identity
//         q = identity;
//       }
      
//       // Put total into correct bucket
//       const size_t i = bucketIndex(total);
//       if (i >= buckets.size()) buckets.resize(i + 1, identity);
//       buckets[i] = total;
//       return total;
//     }
//   };
// }
