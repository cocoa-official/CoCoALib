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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/geobucket.C,v 1.23 2023/12/10 20:27:29 abbott Exp $
// $Log: geobucket.C,v $
// Revision 1.23  2023/12/10 20:27:29  abbott
// Summary: Added commented out template code for geobuckets
//
// Revision 1.22  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.21  2020/02/11 16:56:42  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.20  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.19  2019/11/14 17:57:33  abbott
// Summary: Removed cruft
//
// Revision 1.18  2018/04/10 14:30:49  bigatti
// -- fixed includes
//
// Revision 1.17  2018/04/09 16:34:25  bigatti
// --
// back to original name
//
// Revision 1.2  2018/04/06 15:35:53  bigatti
// -- fixed includes
//
// Revision 1.1  2018/04/06 15:19:24  bigatti
// -- renaming geobucket.C
//
// Revision 1.15  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.14  2015/05/04 13:17:30  bigatti
// -- added MoveLMToBack
//
// Revision 1.13  2015/04/24 15:40:58  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.12  2014/07/07 13:27:05  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.11  2014/04/30 16:10:38  bigatti
// -- changed size/size_t into len/long
//
// Revision 1.10  2014/01/28 13:09:40  bigatti
// -- simplified function content
//
// Revision 1.9  2012/10/24 12:24:22  abbott
// Changed return type of LC.
//
// Revision 1.8  2012/10/16 09:55:27  abbott
// Replaced  RefRingElem  by  RingElem&  (several times)
//
// Revision 1.7  2010/12/26 13:04:37  abbott
// Changed "GlobalXXXput" into corresponding std C++ stream
// (even in commented out code).
//
// Revision 1.6  2009/05/14 09:21:47  abbott
// Changed empty for loop into a while loop -- makes code more readable.
//
// Revision 1.5  2008/05/30 12:45:22  abbott
// Corrected indentation.  Added std:: prefix to size_t.
//
// Revision 1.4  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.3  2007/09/24 14:26:08  abbott
// Added parens to shut up gcc-4.3; might actually make code more readable.
// Added a few spaces here and there (definitely makes code more readable :-)
//
// Revision 1.2  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.5  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.4  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.3  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/06/20 17:25:27  cocoa
// -- added function geobucket::myAddMul(monom, g, gLen);   [without SkipLMFlag]
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.11  2006/04/27 15:58:09  cocoa
// -- minor tidying
//
// Revision 1.10  2006/04/27 14:15:07  cocoa
// -- added constant [ gbk_numbuckets = 20 ] to avoid realloc of geobuckets
//
// Revision 1.9  2006/04/10 10:24:12  cocoa
// -- changed all unsigned int/short into size_t (ex-bug for 4x4)
//
// Revision 1.8  2006/03/21 13:47:03  cocoa
// -- removed wrong CoCoA_ASSERT from MoveLM
//
// Revision 1.7  2006/03/16 12:43:01  cocoa
// -- changed: myMul, myDiv --> myMulByCoeff, myDivByCoeff
//
// Revision 1.6  2006/03/01 14:24:17  cocoa
// -- removed DivMask from geobuckets (not generally useful)
//
// Revision 1.5  2006/01/17 18:18:50  cocoa
// -- changed: fscale, gscale --> gbkScale, gScale
// -- minor cleaning
//
// Revision 1.4  2006/01/17 15:44:56  cocoa
// -- chamges by Max for operations with modules
//
// Revision 1.3  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.2  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.2  2005/07/01 16:08:15  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.3  2005/03/11 18:22:41  cocoa
// -- removed: ComputeFScaleAndGScale (now we use GcdQuot)
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.10  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.9  2004/11/11 13:24:01  cocoa
// -- change: PrintLengths now takes an ostream as first argument
//
// Revision 1.8  2004/10/29 16:17:11  cocoa
// -- new function LPPwMask  (a bit dirty for allowing geobuckets without DivMask
// -- new constructor with DivMask::base
//
// Revision 1.7  2004/09/29 09:21:25  cocoa
// -- fixed the fscale/gscale bug introduced by John
//
// Revision 1.6  2004/07/27 16:03:38  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.5  2003/10/09 13:32:16  cocoa
// A few glitches which slipped through the first major merge.
//
// Revision 1.4  2003/10/08 14:08:51  cocoa
// rubbish comment
//
// Revision 1.3  2003/09/30 09:47:09  cocoa
// - removed trailing spaces
//
// Revision 1.2  2003/09/30 09:35:26  cocoa
// - new coding convention "my"
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.16  2003/06/23 17:11:20  abbott
// Minor cleaning prior to public release.
//
// Revision 1.15  2003/05/15 12:03:45  bigatti
// - improved ComputeFScaleAndGScale
// - used IsOne
// - Poly(b) --> poly(b)
//
// Revision 1.14  2003/05/14 16:23:47  bigatti
// - removed IamNormalized
// - new ring syntax
//
// Revision 1.13  2002/11/18 18:06:05  bigatti
// - code for reduction on GCD rings
// - removed default length in some function calls
// - removed flag  IamNormalized
//
// Revision 1.12  2002/09/19 17:12:21  bigatti
// - Optimizations with: myApproxLen, myMaxLen, IamNormalized
//
// Revision 1.11  2002/05/14 15:28:08  bigatti
// - CascadeFrom now cascades only if len exceeds allowed lengths*2
//
// Revision 1.10  2002/04/29 16:10:48  bigatti
// - the "content divisor" (myCoeff) in bucket now works
//   (useful for using geobuckets in GCD rings)
// - some tidying (normalize(), ComputeFScaleAndGScale(..), ...)
//
// Revision 1.9  2002/04/15 15:16:48  bigatti
// - new function: cascade
// - class bucket made public
// - SetLPP adds LC*LPP into myBuckets[0]
// - myBuckets[0] of length gbk_minlen (used to be just 1)
// - works for GCD rings (mul not optimized yet)
//
// Revision 1.8  2002/03/21 15:18:17  bigatti
// - new: bucket class
//
// Revision 1.7  2002/01/10 13:21:46  bigatti
// - new structure of reduction
//
// Revision 1.6  2001/12/07 13:11:18  bigatti
// - applied coding conventions
//
// Revision 1.5  2001/12/04 17:42:28  bigatti
// - updated calls to DMP (coding conventions)
//
// Revision 1.4  2001/11/30 21:13:45  bigatti
// - changed names: myBuckets, myR, myPPM, ...
// - positive experiments with new merge (to be saved in DMP) and
//   precomputed len of reducers (to be stored in GPoly)
// - added init_clear(DMP& f)
//
// Revision 1.3  2001/11/16 19:23:18  bigatti
// added:  using namespace std;
// for compatibility with gcc-3
// and added  buckets.push_back(zero);  in init
//
// Revision 1.2  2001/11/09 20:51:33  abbott
// Fixed a wrong local variable name -- code compiled and works now! :-)
//
// Revision 1.1  2001/11/07 12:54:41  abbott
// Initial revision
//

