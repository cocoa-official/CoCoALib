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

#endif
