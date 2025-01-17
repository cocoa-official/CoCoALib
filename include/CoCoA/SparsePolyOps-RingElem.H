#ifndef CoCoA_SparsePolyOps_RingElem_H
#define CoCoA_SparsePolyOps_RingElem_H

//   Copyright (c)  2005-2007,2009,2012-2018  John Abbott and Anna M. Bigatti

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


//#include "CoCoA/PPMonoid.H" // included by SparsePolyRing
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"
#include "CoCoA/SparsePolyRing.H" // for IsSparsePolyRing

#include <iosfwd>
// using std::ostream;
#include <vector>
// using std::vector;

namespace CoCoA
{

  class SparsePolyIter;  // forward declaration -- defined in SparsePolyIter.H

  //-- Pseudo constructors for RingElems owned by SparsePolyRing -----------
  RingElem monomial(const SparsePolyRing& P, ConstRefPPMonoidElem pp);
  RingElem monomial(const SparsePolyRing& P, ConstRefRingElem c, ConstRefPPMonoidElem pp);
  RingElem monomial(const SparsePolyRing& P, const MachineInt& n, ConstRefPPMonoidElem pp);
  RingElem monomial(const SparsePolyRing& P, const BigInt& N, ConstRefPPMonoidElem pp);
  RingElem monomial(const SparsePolyRing& P, const BigRat& N, ConstRefPPMonoidElem pp);

  RingElem monomial(const SparsePolyRing& P, const std::vector<long>& expv);
  RingElem monomial(const SparsePolyRing& P, ConstRefRingElem c, const std::vector<long>& expv);
  RingElem monomial(const SparsePolyRing& P, const MachineInt& n, const std::vector<long>& expv);
  RingElem monomial(const SparsePolyRing& P, const BigInt& N, const std::vector<long>& expv);
  RingElem monomial(const SparsePolyRing& P, const BigRat& N, const std::vector<long>& expv);

  RingElem RandomLinearForm(const ring& P);
  RingElem RandomLinearForm(const ring& P, long N);
  //RingElem RandomLinearForm(const ring& P, long lo, long hi); // anon nmsp
  

  RingElem& PushFront(RingElem& f, ConstRefRingElem c, const std::vector<long>& expv); /// SHOULD BE vector<BigInt> ????
  RingElem& PushBack(RingElem& f, ConstRefRingElem c, const std::vector<long>& expv); /// SHOULD BE vector<BigInt> ????
  RingElem& PushFront(RingElem& f, ConstRefRingElem c, ConstRefPPMonoidElem pp); ///< with checks
  RingElem& PushBack(RingElem& f, ConstRefRingElem c, ConstRefPPMonoidElem pp); ///< with checks
  //@}

  RingElem ClearDenom(const SparsePolyRing& ZZx, const RingElem& f); // see also PolyRing.H


  //-- SparsePolyIter -----------------------------
  SparsePolyIter BeginIter(ConstRefRingElem f);
  SparsePolyIter EndIter(ConstRefRingElem f);

  //-- Special functions on RingElems owned by SparsePolyRing --------------
  long UnivariateIndetIndex(ConstRefRingElem f);  // returns -1  if not univariate
  ConstRefPPMonoidElem PPContent(ConstRefRingElem f);     // gcd of all PPs in support
  RingElem IndetsProd(ConstRefRingElem f);             // (monomial) prod of indets appearing in f
  std::vector<long> IndetsIn(ConstRefRingElem f);              // indices of indets appearing in f
  ConstRefPPMonoidElem LPP(ConstRefRingElem f);           // inline
  degree wdeg(ConstRefRingElem f);
  int CmpWDeg(ConstRefRingElem f1, ConstRefRingElem f2);
  int CmpWDegPartial(ConstRefRingElem f1, ConstRefRingElem f2, long GrDim);  // assumes 0 <= GrDim
  bool IsHomog(ConstRefRingElem f);
  bool IsHomogPartial(ConstRefRingElem f, long GrDim);    // assumes 0 <= GrDim
  RingElem homog(ConstRefRingElem f, ConstRefRingElem h); ///< homogenize f with indet h (GrDim=1, wdeg(h)=1)
  RingElem LF(ConstRefRingElem f);
  RingElem CutLF(RingElem& f);  ///< cut off leading form; **MODIFIES** f
  RingElem HomogCompt(ConstRefRingElem f, long d);


  struct CoeffPP
  {
    CoeffPP(const RingElem& c, const PPMonoidElem& t): myCoeff(c), myPP(t) {}
    RingElem myCoeff;
    PPMonoidElem myPP;
  };
  std::ostream& operator<<(std::ostream& out, const CoeffPP& term);

  RingElem ConstantCoeff(ConstRefRingElem f);

  RingElem ContentWRT(ConstRefRingElem f, const std::vector<long>& indets);
  RingElem ContentWRT(ConstRefRingElem f, ConstRefRingElem x);
  void CRTPoly(RingElem& f, BigInt& M,
               ConstRefRingElem f1, const BigInt& M1,
               ConstRefRingElem f2, const BigInt& M2);
  RingElem RatReconstructPoly(ConstRefRingElem fCRT, const BigInt& modulus);


  //-----------------------------------------------
  RingElem CoeffHeight(ConstRefRingElem f);

  bool IsPalindromic(ConstRefRingElem f); // univariate only
  RingElem reverse(ConstRefRingElem f); // univariate only
  RingElem reverse(ConstRefRingElem f, ConstRefPPMonoidElem t);

  RingElem graeffe(ConstRefRingElem f);
  RingElem graeffe3(ConstRefRingElem f);



  //----------------------------------------------------------------------
  // Inline functions

  inline ConstRefPPMonoidElem LPP(ConstRefRingElem f)
  {
    if (!IsSparsePolyRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (IsZero(f))  CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    return SparsePolyRingPtr(owner(f))->myLPP(raw(f));
  }


  //----------------------------------------------------------------------
  // vector<RingElem> functions (implemented in SparsePolyOps-vector.C)

  RingElem IndetsProd(const std::vector<RingElem>& L); // (monomial) prod of indets appearing in elems of (non-empty) L
  std::vector<long> IndetsIn(const std::vector<RingElem>& L);  // indices of indets appearing in elems of (non-empty) L

  std::vector<RingElem> interreduced(std::vector<RingElem> v);
  // NOT DEFINED!!  void interreduce(std::vector<RingElem>& v);
  RingElem NR(ConstRefRingElem f, const std::vector<RingElem>& v); ///< Normal Remainder (remainder of the Groebner Division Algorithm).  If v is a GBasis it gives the Normal Form.
  std::vector<RingElem> TmpDivAlg(ConstRefRingElem f, const std::vector<RingElem>& v);
  // the remainder is in result[len(v)]   (this interface might change)

  std::vector<CoeffPP> CoefficientsWRT(ConstRefRingElem f, const std::vector<long>& indets);
  std::vector<CoeffPP> CoefficientsWRT(ConstRefRingElem f, ConstRefRingElem x);
  std::vector<RingElem> CoeffVecWRT(ConstRefRingElem f, ConstRefRingElem x);
  std::vector<RingElem> CoeffVecWRTSupport(ConstRefRingElem f, ConstRefRingElem basis);
  


} // end of namespace CoCoA

#endif
