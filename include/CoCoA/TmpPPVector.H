#ifndef CoCoA_TmpPPVector_H
#define CoCoA_TmpPPVector_H

//   Copyright (c)  2008-2019  John Abbott, Anna M. Bigatti

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

#include "CoCoA/PPWithMask.H"
#include "CoCoA/assert.H"
#include "CoCoA/utils.H"

#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{
  class ConstRefPPMonoidElem;   // fwd decl -- defined in PPMonoid.H
  class PPWithMask;   // fwd decl -- defined in PPWithMask.H
  class PPMonoid;   // fwd decl -- defined in PPMonoid.H
  class RingElem;   // fwd decl -- defined in ring.H
  class SparsePolyRing;   // fwd decl -- defined in SparsePolyRing.H
  class ring;   // fwd decl -- defined in ring.H


  class PPVector
  {
  public:
    PPVector(PPMonoid PPM, DivMaskRule DMR);
    PPVector(PPMonoid PPM, DivMaskRule DMR, const std::vector<RingElem>& v);
    ~PPVector() {};
    PPVector& operator=(const PPVector& rhs);
    PPVector(const PPVector& copy);
  private:
  public:
    //      void myOutputSelf(std::ostream& out) const;
    //      void myOutputSelf_OM(OpenMathOutput& OMOut) const;

    long myLen() const { return len(myVec); } ///< size of the vector as a long
    const PPWithMask& operator[](long i) const { CoCoA_ASSERT(0 <= i && i < myLen()); return myVec[i]; }
    //    PPVector* myZeroClone() const; ///< creates a PPVector of the same kind
    void myReserve(long n); ///< calls "reserve"
    void myClear(); ///< calls "clear"
    bool IamEmpty() const; ///< return true if vector is empty
    void myPushBackPopBack(PPVector& FromPPs); ///< push_back(pop_back(FromPPs))
    void myPushBack(ConstRefPPMonoidElem pp); ///< push_back(pp), if owner(pp)!=PPM(PPs) it maps it
    void myPushBack(const PPWithMask& pp); ///< push_back(pp)
    void mySort(); ///< sort
    bool myDivides(const PPWithMask& pp) const; ///< true if pp is divisible by an element of *this
    bool myDivides(ConstRefPPMonoidElem pp) const; ///< true if pp is divisible by an element of *this
    void myInterreduce(); ///< interreduce PPs (NOT exception clean)
    void myInterreduceSort(); ///< interreduce and sort PPs (NOT exception clean)
    void myLcms(const PPVector& PPs1, const PPVector& PPs2); ///< all the lcm between elements of PPs1 and PPs2, effectively the generators of the intersection ideal
    void myAlexanderDual(const PPVector& PPs); ///< only for squarefree ideals (AMB:2011-06-27) 
    void mySupport(ConstRefPPMonoidElem pp); ///< vector of indets in pp
    void myConvert(std::vector<RingElem>& v, SparsePolyRing P) const; ///< converts PP's into RingElem's

    friend void swap(PPVector& PPs1, PPVector& PPs2);
    friend PPMonoid PPM(const PPVector& PPs);
    friend DivMaskRule DMR(const PPVector& PPs);
    friend std::ostream& operator<<(std::ostream&, PPVector);
//     friend void PushBack(PPVector& PPs, ConstRefPPMonoidElem pp);
//     friend void PushBackPopBack(PPVector& ToPPs, PPVector& FromPPs);
//     friend bool IsDivisible(const PPWithMask& pp, const PPVector& ByL);
//     friend bool IsDivisible(const PPWithMask& pp, const PPVector& ByL);
//     friend void interreduce(PPVector& PPs);
//     friend void InterreduceSort(PPVector& PPs);
//     friend void lcms(PPVector& PPs, const PPVector& PPs1, const PPVector& PPs2);
  private: // member fields
    PPMonoid myPPM;
    DivMaskRule myDMR;
    std::vector<PPWithMask> myVec;
  };

  // inline ----------------------------------------------------
  inline PPMonoid PPM(const PPVector& PPs) { return PPs.myPPM; }
  inline DivMaskRule DMR(const PPVector& PPs) { return PPs.myDMR; }
  inline bool IsEmpty(const PPVector& PPs) { return PPs.IamEmpty(); }
  inline long len(const PPVector& PPs) { return PPs.myLen(); }

  void PushBack(PPVector& PPs, ConstRefPPMonoidElem pp);
  void PushBack(PPVector& PPs, const PPWithMask& pp);
  void PushBackPopBack(PPVector& ToPPs, PPVector& FromPPs); ///< append(ToPPs, last(FromPPs)), and remove the element from FromPPs
  bool IsDivisible(const PPWithMask& pp, const PPVector& ByL); ///< true if pp is divisible by an element of L
  bool IsDivisible(ConstRefPPMonoidElem pp, const PPVector& ByL); ///< true if pp is divisible by an element of L
  void interreduce(PPVector& PPs); ///< interreduce PPs (NOT exception clean)
  void InterreduceSort(PPVector& PPs); ///< interreduce and sort PPs (NOT exception clean)
  void lcms(PPVector& PPs, const PPVector& PPs1, const PPVector& PPs2); ///< all the lcm between elements of PPs1 and PPs2, effectively the generators of the intersection ideal
  void convert(std::vector<RingElem>& v, ring P, const PPVector& PPs); ///< converts PP's into RingElem's
  void convert(PPVector& PPs, const std::vector<RingElem>& v); ///< converts RingElem's into PP's
  void swap(PPVector& PPs1, PPVector& PPs2); ///< swap PPs1 and PPs1


  // void pivot(std::vector<PPWithMask>& PPs); // anna: for Hilbert not in PPs (pp or sp), for mvt in PPs and needs the position for deletion
  // void sort(std::vector<PPWithMask>& PPs);
  // void colon(std::vector<PPWithMask>& PPs, const std::vector<PPWithMask>& PPs1, const std::vector<PPWithMask>& PPs2);
  // void IndetsFrequencies(std::vector<int>& freq, const std::vector<PPWithMask>& PPs);

  //----

// PPVector SplitIndets(PPVector PPs);
// void ReduceAndDivideBySimplePower( PPVector PPs, PPVector *DivTList,
//                                    long PIndex, unsigned PExp);
// void ReduceAndDivideByPivot(PPVector PPs, PPVector *DivTList, PPWithMask Pivot);
// void myBigPivot(int *PivotIndex, int *PivotDeg);
// void myGCD3Pivot(PPWithMask pp);

// void MoveNotCoprimeSP(PPWithMask FromSPList, PPWithMask ToSPList, PPWithMask theTerm);

//void InsInTList(PPVector PPs, PPWithMask pp, int *NewLen);

// inline void PPVector::myEraseNth(int N)
// {
//   swap(myPPs[N], myPPs.last());
//   pop_back();
// }


// /********************    list of simple powers    ********************/
// /*
//     SPList(PPs) is represented with a PPWithMask
//     which is the product of the simple powers.
// */
// void InsInSPList(long Index, unsigned Exp, PPWithMask SPL);


// /**********************************************************************/
}

#endif /* PPVECTOR_H */
