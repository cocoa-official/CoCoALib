#ifndef CoCoA_geobucket_H
#define CoCoA_geobucket_H

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


#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H

#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{
 
  class geobucket
  {
  public:

    class bucket
    {
      friend class geobucket;
    public:
      bucket(const SparsePolyRing&, long MaxLen);
      bucket(const bucket&);///< copy constructor to use C++ vectors of buckets

      ~bucket();

    private:
      void myNormalize();    ///< myPoly *= myCoeff; myCoeff = 1
      /**< \li \em weak exception guarantee */
      void myAddClear(RingElem& f, long FLen);    ///< *this += f; f = 0
      /**< \li *this will be normalized
         \li \em weak exception guarantee */
      void myAddClear(bucket& b);                 ///< *this += b; b = 0
      /**< \li *this will be normalized
         \li \em weak exception guarantee */
      void myMulByCoeff(ConstRefRingElem coeff);  ///< *this *= coeff
      /**< \li \em weak exception guarantee */
      void myDivByCoeff(ConstRefRingElem coeff);  ///< *this /= coeff
      /**< \li assumes *this is divisible by coeff
         \li \em weak exception guarantee */

      friend bool IsZero(const bucket& b);
      friend RingElem content(const bucket& b);
      friend ConstRefRingElem poly(bucket& b);    ///< it normalizes the bucket and returns a reference to the polynomial

      ///@name Dirty functions for efficiency
      //@{
      static bool myIsZeroAddLCs(const SparsePolyRing&, bucket& b1, bucket& b2);  ///<  b1 += LM(b2);  b2 -= LM(b2);  return LC(b1)+LC(b2)==0
      /**< \li it assumes LPP(b1) == LPP(b2)
         \li b1 and b2 will be normalized */
      friend void MoveLMToFront(bucket& b1, bucket& b2); ///< b1 += LM(b2); b2 -= LM(b2)
      friend void MoveLMToFront(const SparsePolyRing&, bucket& b1, bucket& b2); ///< b1 += LM(b2); b2 -= LM(b2)
      /**< \li it assumes LPP(b1)<LPP(b2)
         \li b1 and b2 will be normalized */
      //@}

      ///@name Friends functions on geobuckets
      //@{
      friend ConstRefPPMonoidElem LPP(const geobucket& gbk);
      friend void AddClear(RingElem& f, geobucket& gbk);
      friend void ReductionStep(geobucket& gbk, ConstRefRingElem f, long RedLen);
      friend void MoveLMToFront(RingElem& f, geobucket& gbk);
      friend void MoveLMToBack(RingElem& f, geobucket& gbk);

      friend std::ostream& operator<<(std::ostream& out, const geobucket& g);
      friend void PrintLengths(std::ostream& out, const geobucket& g);
      //@}


      ///@name data members of geobucket::bucket
      //@{
    private: // data members
      RingElem myCoeff;  ///< the coefficient factor
      RingElem myPoly;   ///< the polynomial (a "clean" RingElem)
      long myMaxLen;     ///< the maximal length allowed for the polynomial of this bucket
      long myApproxLen;  ///< an upper bound for the current length of the polynomial of this bucket
      //@}

    }; // class geobucket::bucket


  public:
    geobucket(const SparsePolyRing&);
    ~geobucket();

    friend std::ostream& operator<<(std::ostream& out, const geobucket& gbk);
    friend void PrintLengths(std::ostream& out, const geobucket& gbk); ///< just for debugging

    long myLen() const;
    void myPushBackZeroBucket(long MaxLen);
    void myAddClear(RingElem& f, long len);
    void myDeleteLM();
    long myBucketIndex(long len); ///< \return the bucket index for a polynomial of length \a len
    void myAddMulLM(ConstRefRingElem monom, ConstRefRingElem g, long gLen); ///< *this += monom*g
    void myAddMulLM(ConstRefRingElem monom, ConstRefRingElem g, long gLen, SparsePolyRingBase::SkipLMFlag); ///< *this += monom*g

    friend RingElem content(const geobucket& gbk);
    friend void RemoveBigContent(geobucket& gbk);
    void myDivByCoeff(ConstRefRingElem coeff); ///< content MUST be divisible by \a coeff
    void myMulByCoeff(ConstRefRingElem coeff);
    friend const ring& CoeffRing(const geobucket& gbk);
    friend const PPMonoid& PPM(const geobucket& gbk);
    friend void AddClear(RingElem& f, geobucket& gbk);
    friend bool IsZero(const geobucket& gbk);
    friend RingElemAlias LC(const geobucket& gbk);
    friend ConstRefPPMonoidElem LPP(const geobucket& gbk);
    void myCascadeFrom(long i);
    friend void MoveLMToFront(RingElem& g, geobucket& gbk);
    friend void MoveLMToBack(RingElem& g, geobucket& gbk);

    friend void ReductionStep(geobucket& gbk, ConstRefRingElem g, long RedLen);
    friend void ReductionStepGCD(geobucket& gbk, ConstRefRingElem g, RingElem& FScale, long RedLen);

    ///@name data members of geobucket
    //@{
  private: // data members
    SparsePolyRing myPolyRing; ///< the SparsePolyRing gbk lives in
    mutable  bool IhaveLM; ///< true if certified  that LM(gbk) = LM(gbk[0])
    mutable  std::vector<bucket> myBuckets;  ///< the bucket vector
    //@}

    void mySetLM() const;  /**< Sets the LM of *this in the 0-th bucket and set IhaveLM to true
                            *this will be normalized */

  }; // end of class geobucket

  std::ostream& operator<<(std::ostream& out, const geobucket& gbk);
  long len(const geobucket& gbk);
  void AddClear(RingElem& f, geobucket& gbk);
  void ReductionStep(geobucket& gbk, ConstRefRingElem g, long RedLen);
  void ReductionStepGCD(geobucket& gbk, ConstRefRingElem g, RingElem& FScale, long RedLen);

  //----------------------------------------------------------------------//
  // inline functions
  //----------------------------------------------------------------------//

  //----------  bucket functions  ----------//

  inline bool IsZero(const geobucket::bucket& b)
  { return IsZero(b.myPoly); }

  //----------  geobucket functions  ----------//


  inline ConstRefPPMonoidElem LPP(const geobucket& gbk)
  {
    if (!gbk.IhaveLM) gbk.mySetLM();
    return LPP(gbk.myBuckets[0].myPoly);
    //    return gbk.myPolyRing->myLPP(raw(gbk.myBuckets[0].myPoly));
  }


  inline const ring& CoeffRing(const geobucket& gbk)
  { return CoeffRing(gbk.myPolyRing); }


  inline const PPMonoid& PPM(const geobucket& gbk)
  { return PPM(gbk.myPolyRing); }

  // must be inline for efficiency (I checked!)
  inline bool IsZero(const geobucket& gbk)
  {
    if (!gbk.IhaveLM) gbk.mySetLM();
    return IsZero(gbk.myBuckets[0]);
  }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/geobucket.H,v 1.12 2022/02/18 14:11:52 abbott Exp $
// $Log: geobucket.H,v $
// Revision 1.12  2022/02/18 14:11:52  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2018/05/18 16:39:39  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.10  2018/04/09 16:32:54  bigatti
// -- back to original name
//
// Revision 1.1  2018/04/06 15:12:02  bigatti
// -- renamed geobucket.H
//
// Revision 1.8  2016/05/18 12:19:34  abbott
// Summary: Gave consistent names to args in fn prototypes
//
// Revision 1.7  2015/05/04 12:53:59  bigatti
// -- new friend: MoveLMToBack
// -- removed doxygen headers (now in doc)
//
// Revision 1.6  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.5  2014/04/30 16:10:38  bigatti
// -- changed size/size_t into len/long
//
// Revision 1.4  2012/10/24 12:08:41  abbott
// Changed return type of LC.
//
// Revision 1.3  2012/10/16 09:55:27  abbott
// Replaced  RefRingElem  by  RingElem&  (several times)
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/02/12 19:04:49  cocoa
// Very minor mods to geobucket.
//
// Revision 1.5  2006/12/20 18:58:48  cocoa
// -- only white space differences
//
// Revision 1.4  2006/12/06 17:17:44  cocoa
// -- removed #include "config.H"
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.2  2006/06/20 17:25:27  cocoa
// -- added function geobucket::myAddMul(monom, g, gLen);   [without SkipLMFlag]
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.12  2006/04/27 15:57:43  cocoa
// -- minor tidying
//
// Revision 1.11  2006/04/27 14:15:07  cocoa
// -- added constant [ gbk_numbuckets = 20 ] to avoid realloc of geobuckets
//
// Revision 1.10  2006/04/27 13:42:26  cocoa
// Parameter name change for better readability.
//
// Revision 1.9  2006/04/10 10:24:12  cocoa
// -- changed all unsigned int/short into size_t (ex-bug for 4x4)
//
// Revision 1.8  2006/03/24 14:21:49  cocoa
// -- added declarations of ReductionStep, and ReductionStepGCD
//
// Revision 1.7  2006/03/16 12:43:01  cocoa
// -- changed: myMul, myDiv --> myMulByCoeff, myDivByCoeff
//
// Revision 1.6  2006/03/13 16:57:49  cocoa
// -- changed: member data  myPolyRing  is no longer a reference
//
// Revision 1.5  2006/03/07 16:50:44  cocoa
// -- changed: LPP returns a ConstRefPPMonoidElem (as all LPP functions)
//
// Revision 1.4  2006/03/01 14:24:17  cocoa
// -- removed DivMask from geobuckets (not generally useful)
//
// Revision 1.3  2006/01/17 10:23:08  cocoa
// Updated DivMask; many consequential changes.
// A few other minor fixes.
//
// Revision 1.2  2005/12/31 12:22:18  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.3  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//
// Revision 1.2  2005/06/22 14:47:56  cocoa
// PPMonoids and PPMonoidElems updated to mirror the structure
// used for rings and RingElems.  Many consequential changes.
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.8  2004/11/11 13:24:01  cocoa
// -- change: PrintLengths now takes an ostream as first argument
//
// Revision 1.7  2004/11/02 14:49:50  cocoa
// -- just a little fix for doxygen (\return -> return)
//
// Revision 1.6  2004/10/29 16:17:11  cocoa
// -- new function LPPwMask  (a bit dirty for allowing geobuckets without DivMask
// -- new constructor with DivMask::base
//
// Revision 1.5  2004/09/20 15:31:51  cocoa
// -- only an "\a" for doxygen
//
// Revision 1.4  2004/07/27 16:03:38  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.3  2003/10/09 13:32:15  cocoa
// A few glitches which slipped through the first major merge.
//
// Revision 1.2  2003/09/30 09:35:26  cocoa
// - new coding convention "my"
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
// Revision 1.15  2003/06/27 10:59:09  bigatti
// - added short documentation for geobucket fields
//
// Revision 1.14  2003/06/23 17:11:08  abbott
// Minor cleaning prior to public release.
// Reindented.
//
// Revision 1.13  2003/05/15 12:04:39  bigatti
// - Poly(b) --> poly(b)
//
// Revision 1.12  2003/05/14 16:34:49  bigatti
// - new ring syntax
// - doxygen documentation syntax
// - ..LPP --> ..LM
// - cleaning
//
// Revision 1.11  2002/11/18 18:03:26  bigatti
// - code for reduction on GCD rings
// - removed default length in some funcion calls
// - removed flag  IamNormalized
//
// Revision 1.10  2002/09/19 17:14:25  bigatti
// - Optimizations with: myApproxLen, myMaxLen, IamNormalized
//
// Revision 1.9  2002/04/29 16:13:04  bigatti
// - added function: "bucket::normalize()"
// - modified function: "CascadeFrom" [was "cascade"]
//
// Revision 1.8  2002/04/15 15:18:37  bigatti
// - new function: cascade
// - class bucket made public
// - SetLPP adds LC*LPP into myBuckets[0]
// - myBuckets[0] of length gbk_minlen (used to be just 1)
// - works for GCD rings (mul not optimized yet)
//
// Revision 1.7  2002/03/21 15:21:07  bigatti
// - new: bucket class
//
// Revision 1.6  2002/01/10 13:20:23  bigatti
// - new structure of reduction
//
// Revision 1.5  2001/12/07 13:12:09  bigatti
// - applied coding conventions
//
// Revision 1.4  2001/12/04 17:39:23  bigatti
// - changed "IsZero"
//
// Revision 1.3  2001/11/30 21:16:53  bigatti
// - changed names: myBuckets, myR, myPPM
// - added init_clear(DMP& f)
//
// Revision 1.2  2001/11/16 19:21:43  bigatti
// added:  std::
// for compatibility with gcc-3
//
// Revision 1.1  2001/11/07 12:54:33  abbott
// Initial revision
//

#endif
