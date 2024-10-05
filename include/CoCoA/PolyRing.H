#ifndef CoCoA_PolyRing_H
#define CoCoA_PolyRing_H

//   Copyright (c)  2005,2008,2009,2014,2019  John Abbott, Anna M. Bigatti

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


// Header file for the abstract class PolyRing


#include "CoCoA/error.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingElemInput.H"

#include <iosfwd>
//using std::ostream;
#include <string>
//using std::string;
#include <vector>
//using std::vector;


namespace CoCoA
{

  class PolyRingBase: public RingBase
  {
  private: // This class has no data members.
  public:
    ///@name functions inherited from ring with a single implementation for all PolyRing implementations
    //@{
    const ring& myBaseRing() const override  { return myCoeffRing(); }
    BigInt myCharacteristic() const override; ////  { return characteristic(myCoeffRing()); }
    bool IamCommutative() const override  { return true; } ///< RingWeyl redefines this
    bool3 IamIntegralDomain3(bool QuickMode) const override  { return myCoeffRing()->IamIntegralDomain3(QuickMode); }
    bool IamTrueGCDDomain() const override  { return myCoeffRing()->IamTrueGCDDomain() || myCoeffRing()->IamField(); };
    bool IamField() const override  { return false; } //??? (myNumIndets()==0 && myCoeffRing()->IamField())
    bool IamFiniteField() const override  { return false; } //??? (myNumIndets()==0 && myCoeffRing()->IamFiniteField())
    bool IamExact() const override  { return IsExact(myCoeffRing()); }
    virtual std::string myImplDetails() const = 0;
    void myOutputSelf(std::ostream& out) const override;     // out << R
    void myOutputSelfLong(std::ostream& out) const override; // out << R
//???    virtual void myOutputSelf_OM(OpenMathOutput& OMOut) const = 0; ///< OMOut << R    //@}

    /*----------------------------------------------------------------------
     Member functions every concrete PolyRing must have
     in addition to those of a plain ring.
    ----------------------------------------------------------------------*/

    virtual long myNumIndets() const = 0;
    virtual const ring& myCoeffRing() const = 0;

    ///@name Simple functions on polynomials
    //@{
    bool myIsIrred(ConstRawPtr rawf) const override; ///< common impl
    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< common impl
    virtual const std::vector<RingElem>& myIndets() const = 0;
    virtual const symbol& myIndetSymbol(long idx) const = 0;
    virtual std::vector<RingElem> myIndets(const std::string& s) const; ///< common impl
    virtual void myIndetPower(RawPtr rawf, long var, long exp) const = 0;  ///< assumes exp >= 0
    virtual long myNumTerms(ConstRawPtr rawf) const = 0;
    virtual bool myIsConstant(ConstRawPtr rawf) const = 0;
    virtual bool myIsIndet(long& index, ConstRawPtr rawf) const = 0;
    virtual bool myIsIndetPosPower(ConstRawPtr rawf) const = 0;
    virtual bool myIsMonomial(ConstRawPtr rawf) const = 0;
    virtual bool myIsEvenPoly(ConstRawPtr rawf) const = 0;
    virtual bool myIsOddPoly(ConstRawPtr rawf) const = 0;

    virtual long myStdDeg(ConstRawPtr rawf) const = 0;
    virtual long myDeg(ConstRawPtr rawf, long var) const = 0;
    virtual RingElemAlias myLC(ConstRawPtr rawf) const = 0;
    virtual void myContent(RawPtr rawcontent, ConstRawPtr rawf) const = 0;
    virtual void myContentFrF(RawPtr rawcontent, ConstRawPtr rawf) const = 0;
    virtual void myCommonDenom(RawPtr rawcontent, ConstRawPtr rawf) const = 0;
    virtual void myClearDenom(RawPtr rawres, ConstRawPtr rawf) const = 0;
    virtual void myRemoveBigContent(RawPtr rawf) const = 0;
    virtual bool myIsValid(ConstRawPtr rawf) const = 0; ///< true iff f is a valid representation
    virtual void myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const = 0;
    virtual bool myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const = 0; ///< true iff divisible
    virtual void myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const = 0; ///< lhs = deriv(f, x)
    void myMonic(RawPtr rawmonic, ConstRawPtr rawf) const; // common impl
    //@}


    virtual RingHom myCoeffEmbeddingHomCtor() const = 0;
    virtual RingHom myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const = 0;

    bool myImageLiesInSubfield(const RingHom& phi) const override;

///    void myCheckIndetIndex(long i, const char* where) const; ///< throws ERR::BadIndetIndex if i is negative or too large
    void myCheckIndetIndex(long i, const ErrorContext& ErrCtx) const; ///< throws ERR::BadIndetIndex if i is negative or too large
  };


  class PolyRing: public ring
  {
    // Data member inherited from class ring.
  public:
    /*implicit*/ PolyRing(const ring& R);
///    PolyRing(const ring& R, const char* const FnName);
    PolyRing(const ring& R, const ErrorContext& ErrCtx);
    explicit PolyRing(const PolyRingBase* RingPtr);
    // Default copy ctor & assignment work fine.
    // Default dtor works fine.
  public:
    const PolyRingBase* operator->() const; // allow member fns to be called
  };


  const PolyRing& RingQQt(const MachineInt& NumIndets);

  bool IsPolyRing(const ring& R);                         // inline
  const ring& CoeffRing(const PolyRing& Rx);              // inline
  long NumIndets(const PolyRing& Rx);                     // inline
  RingHom CoeffEmbeddingHom(const PolyRing& Rx);          // inline
  const std::vector<RingElem>& indets(const PolyRing& Rx);// inline

  const RingElem& indet(const PolyRing& P, long i); ///< returns x[i] as a polynomial
  RingElem IndetPower(const PolyRing& P, long i, long n); ///< returns x[i]^n as a polynomial; error if n < 0
  RingElem IndetPower(const PolyRing& P, long i, const BigInt& N); ///< returns x[i]^N as a polynomial; error if N < 0
  long NumTerms(ConstRefRingElem f);
  bool IsMonomial(ConstRefRingElem f);            ///< f == coeff*pp
  bool AreMonomials(const std::vector<RingElem>& v); ///< all elements are monomials
  bool IsConstant(ConstRefRingElem f);            ///< f == coeff
  bool IsIndet(ConstRefRingElem f);               ///< f == x[i] for some index i
  bool IsIndet(long& index, ConstRefRingElem f);  ///< f == x[i]; index=i
  bool IsIndetPosPower(ConstRefRingElem f);       ///< f == x[i]^d for some i and d>0
  bool IsEvenPoly(ConstRefRingElem f);            ///< all exponents are even
  bool IsOddPoly(ConstRefRingElem f);             ///< all exponents are odd

  long StdDeg(ConstRefRingElem f);                // inline
  long deg(ConstRefRingElem f);                   // inline (identical to StdDeg)
  long deg(ConstRefRingElem f, long var);
  RingElemAlias LC(ConstRefRingElem f);           ///< NB result belongs to CoeffRing(owner(f))
  RingElem content(ConstRefRingElem f);           ///< NB result belongs to CoeffRing(owner(f))
  RingElem FixedDivisor(ConstRefRingElem f);      ///< NB result belongs to CoeffRing(owner(f))
  RingElem CommonDenom(ConstRefRingElem f);       ///< NB result belongs to BaseRing(CoeffRing(owner(f)))
  RingElem CommonDenom(const std::vector<RingElem>& v);
  RingElem ClearDenom(ConstRefRingElem f);        ///< NB result belongs to same ring as f
  RingElem monic(ConstRefRingElem f);             ///< f/LC(f)
  RingElem prim(ConstRefRingElem f);              ///< alpha*f for alpha in QQ s.t. coeffs are integer and coprime
  bool IsSqFree(ConstRefRingElem f);              ///< true iff f is square-free (coeffs must be in field)
  RingElem deriv(ConstRefRingElem f, ConstRefRingElem x); ///< derivative of f w.r.t. x, x must be an indeterminate
  RingElem deriv(ConstRefRingElem f, long var);   ///< here var is the index of the indeterminate

  //----- RingHom --------------------------------------------------
  RingHom PolyRingHom(const PolyRing& domain, const ring& codomain, RingHom CoeffHom, const std::vector<RingElem>& IndetImages); // R[x] -> S, defined by CoeffHom: R->S and x->IndetImages
  RingHom PolyRingHom(const PolyRing& Rx, const ring& S, RingHom CoeffHom, const std::string& IndetImages);
  RingHom EvalHom(const PolyRing& domain, const std::vector<RingElem>& IndetImages); // Maps f in R[x_1..x_n] into f(a_1,...,a_n) in R
  RingHom EvalHom(const PolyRing& domain, const MachineInt& n); // Maps f in R[x] into f(n) in R
  RingHom EvalHom(const PolyRing& domain, const BigInt& N);     // Maps f in R[x] into f(N) in R
  RingHom EvalHom(const PolyRing& domain, const BigRat& q);     // Maps f in R[x] into f(q) in R
  RingHom EvalHom(const PolyRing& domain, ConstRefRingElem r);  // Maps f in R[x] into f(r) in R
  RingHom PolyAlgebraHom(const PolyRing& Rx, const ring& Ry, const std::vector<RingElem>& IndetImages); // R[x] -> R or R[x]->R[y], defined by R->R and x->IndetImages
  RingHom PolyAlgebraHom(const PolyRing& Rx, const ring& Ry, const std::string& IndetImages); // R[x] -> R or R[x]->R[y], defined by R->R and x->IndetImages


  //----------------------------------------------------------------------
  // Inline functions

  inline const PolyRingBase* PolyRingPtr(const ring& R)
  { return dynamic_cast<const PolyRingBase*>(R.myRawPtr()); }

  // inline const PolyRingBase* PolyRingPtr(const ring& R, const char* const FnName)
  // {
  //   const PolyRingBase* ptr = PolyRingPtr(R);
  //   if (ptr == nullptr)  CoCoA_THROW_ERROR(ERR::ReqPolyRing, FnName);
  //   return ptr;
  // }

  inline const PolyRingBase* PolyRingPtr(const ring& R, const ErrorContext& ErrCtx)
  {
    const PolyRingBase* ptr = PolyRingPtr(R);
    if (ptr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqPolyRing, ErrCtx);
    return ptr;
  }

  inline bool IsPolyRing(const ring& R)
  {
    return PolyRingPtr(R) != nullptr;
  }


  inline PolyRing::PolyRing(const ring& R):
      ring(PolyRingPtr(R,CoCoA_ERROR_CONTEXT))
  {}

  // inline PolyRing::PolyRing(const ring& R, const char* const FnName):
  //     ring(PolyRingPtr(R, FnName))
  // {}

  inline PolyRing::PolyRing(const ring& R, const ErrorContext& ErrCtx):
      ring(PolyRingPtr(R, ErrCtx))
  {}

  inline PolyRing::PolyRing(const PolyRingBase* RingPtr):
      ring(RingPtr)
  {}


  inline const PolyRingBase* PolyRing::operator->() const
  {
    return static_cast<const PolyRingBase*>(ring::operator->());
  }


  inline const ring& CoeffRing(const PolyRing& Rx)
  {
    return Rx->myCoeffRing();
  }


  inline long NumIndets(const PolyRing& Rx)
  {
    return Rx->myNumIndets();
  }


  inline const std::vector<RingElem>& indets(const PolyRing& Rx)
  {
    return Rx->myIndets();
  }


  inline const symbol& IndetSymbol(const PolyRing& Rx, long idx)
  {
    return Rx->myIndetSymbol(idx);
  }


  inline std::vector<RingElem> indets(const PolyRing& Rx, const std::string& s)
  {
    return Rx->myIndets(s);
  }


  inline long StdDeg(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsZero(f))
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    return PolyRingPtr(owner(f))->myStdDeg(raw(f));
  }

  inline long deg(ConstRefRingElem f)
  {
    return StdDeg(f);
  }


  inline RingElemAlias LC(ConstRefRingElem f)
  {
    if (!IsPolyRing(owner(f)))
      CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsZero(f))
      CoCoA_THROW_ERROR1(ERR::ReqNonZeroRingElem);
    return PolyRingPtr(owner(f))->myLC(raw(f));
  }


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/PolyRing.H,v 1.64 2024/07/31 07:52:46 abbott Exp $
// $Log: PolyRing.H,v $
// Revision 1.64  2024/07/31 07:52:46  abbott
// Summary: Switched to new CoCoA_THROW_ERROR macros (with suffix 1 or 2)
//
// Revision 1.63  2024/07/22 15:43:57  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.62  2024/07/02 20:42:04  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.61  2024/07/02 09:57:56  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.60  2024/03/16 10:15:11  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.59  2024/03/15 19:04:45  abbott
// Summary: Arg checking fns now expect a ErrorContext object (instead of string; related to redmine 1726)
//
// Revision 1.58  2023/03/08 19:43:56  abbott
// Summary: Made myCharacteristic a function (instead of proc with OUT param) redmine 536
//
// Revision 1.57  2022/02/18 14:11:47  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.56  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.55  2022/02/02 09:40:09  abbott
// Summary: Removed discriminant to SparsePolyOps-resultant
//
// Revision 1.54  2021/10/30 19:03:08  abbott
// Summary: Added some more "override" (clang said they were mising; redmine 1625)
//
// Revision 1.53  2021/10/30 17:39:00  abbott
// Summary: Minor correction
//
// Revision 1.52  2021/10/30 17:14:16  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.51  2021/03/01 19:39:07  abbott
// Summary: Added prim (redmine 1580)
//
// Revision 1.50  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.49  2019/10/10 06:46:25  bigatti
// -- added PolyRingHom PolyAlgebraHom with string (using ReadExprVector)
//
// Revision 1.48  2019/09/16 17:21:48  abbott
// Summary: Added IsEvenPoly, IsOddPoly and FixedDivisor
//
// Revision 1.47  2019/03/18 17:13:24  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.46  2018/03/20 11:40:43  bigatti
// -- added CommonDenom for vector of polynomials over QQ
//
// Revision 1.45  2016/10/25 20:54:09  abbott
// Summary: Added new fn IsSqFree (for BigInt and ringelem of PolyRing over field)
//
// Revision 1.44  2016/04/07 14:35:45  abbott
// Summary: Removed code which disabled assignment of rings
//
// Revision 1.43  2014/07/28 15:16:57  abbott
// Summary: Renamed mem fn myCoeffEmbeddingHomCtor
// Author: JAA
//
// Revision 1.42  2014/07/28 15:10:23  abbott
// Summary: Redesign: ringhoms no longer cached in rings (caused ref count trouble)
// Author: JAA
//
// Revision 1.41  2014/07/14 11:39:04  abbott
// Summary: Added convenient "univariate" EvalHom signatures
// Author: JAA
//
// Revision 1.40  2014/07/11 15:41:56  bigatti
// -- added  myOutputSelf (default impl),  myImplDetails()
//
// Revision 1.39  2014/07/08 13:14:40  abbott
// Summary: Removed AsQuotientRing; added new defn of BaseRing
// Author: JAA
//
// Revision 1.38  2014/07/07 17:11:20  abbott
// Summary: [MAJOR CHANGE] Removed AsPolyRing; added PolyRingPtr
// Author: JAA
//
// Revision 1.37  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.36  2014/04/15 13:28:18  abbott
// Summary: Removed useless IsIrred decl (already decl in ring.H)
// Author: JAA
//
// Revision 1.35  2014/04/02 10:59:19  abbott
// Summary: Revised design of IamIntegralDomain3
// Author: JAA
//
// Revision 1.34  2014/03/27 17:17:31  abbott
// Summary: Added new fn IsIntegralDomain3 (and mem fn IamIntegralDomain3)
// Author: JAA
//
// Revision 1.33  2013/06/18 12:27:37  abbott
// Renamed HibertSeriesPolyRing to RingQQt.
//
// Revision 1.32  2013/03/26 14:58:08  abbott
// Added IndetPower allowing exp to be a BigInt (handy elsewhere).
//
// Revision 1.31  2013/01/17 13:45:22  abbott
// Changed name: IndetsCalled --> indets.
// Added virt mem fn  myImageLiesInSubfield
//
// Revision 1.30  2012/10/24 13:30:07  abbott
// Changed return type of LC/myLC.
//
// Revision 1.29  2012/05/22 10:02:38  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.28  2012/04/27 15:02:23  abbott
// Added myIsFiniteField
//
// Revision 1.27  2012/03/30 16:19:30  bigatti
// -- removed doxygen header
//
// Revision 1.26  2011/11/07 10:29:09  bigatti
// -- added AreMonomials
//
// Revision 1.25  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.24  2011/06/27 13:43:39  bigatti
// -- allowing CoeffRing(Rx) -> CoeffRing(S) as CoeffHom in PolyRingHom
//
// Revision 1.23  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.22  2011/05/26 11:57:48  bigatti
// -- swapped args in IndetsCalled
//
// Revision 1.21  2011/03/16 15:39:34  bigatti
// -- added myIsIndetPosPower(f), IsIndetPosPower(f)
//
// Revision 1.20  2011/03/16 13:17:21  abbott
// Removed indet(PolyRing,ZZ).
// Added commented out IndetPower with ZZ exponent.
//
// Revision 1.19  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.18  2011/03/01 14:10:47  bigatti
// -- added ClearDenom/myClearDenom
//
// Revision 1.17  2011/01/24 17:16:25  bigatti
// -- added monic/myMonic
//
// Revision 1.16  2011/01/18 14:35:38  bigatti
// -- moved **_forC5 functions into CoCoA-5/CoCoALibSupplement:
//    myMonomials_forC5, mySupport_forC5, monomials_forC5, support_forC5,
//    LPP_forC5, LT_forC5, LM_forC5
//
// Revision 1.15  2010/11/30 11:34:01  bigatti
// -- added unique implementation for IndetsCalled
// -- added virtual myIndetSymbol
//
// Revision 1.14  2010/11/25 12:30:13  bigatti
// -- added IndetsCalled/myIndetsCalled
//
// Revision 1.13  2010/11/05 15:54:26  bigatti
// -- added myMonomials_forC5, mySupport_forC5
//
// Revision 1.12  2010/11/02 16:03:12  bigatti
// -- added indet(*, ZZ) [for CoCoA-5]
// -- indet(PPM, i) now returns "const PPMonoidElem&"
//
// Revision 1.11  2009/10/02 13:41:02  bigatti
// -- unique implementation of myDiv in PolyRing.C
// -- myDivByCoeff now returns a bool (if exact division)
//
// Revision 1.10  2009/07/24 12:26:43  abbott
// Added CommonDenom function for polynomials.
//
// Revision 1.9  2008/04/21 12:50:09  abbott
// Added EvalHom, an easy way to create evaluation homomorphisms for poly rings.
//
// Revision 1.8  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.7  2008/03/12 16:02:48  bigatti
// -- added: IsIrred
//
// Revision 1.6  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.5  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.4  2007/05/31 16:32:58  bigatti
// -- added function IsValid for sanity checks
//
// Revision 1.3  2007/05/31 15:26:23  bigatti
// -- default implementation for IamCommutative, IamIntegralDomain,
//    IamGCDDomain, IamField, myCharacteristic
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.13  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.12  2007/03/07 13:59:45  bigatti
// -- added PolyAlgebraHom(Rx, Ry, IndetImages)
//
// Revision 1.11  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.10  2007/02/28 13:58:35  bigatti
// -- added documentation for IsMonomial, IsConstant, IsIndet
//
// Revision 1.9  2007/02/28 13:51:59  bigatti
// -- added function IsMonomial
//
// Revision 1.8  2007/01/11 14:06:15  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.7  2006/12/06 17:15:26  cocoa
// -- removed #include "config.H"
//
// Revision 1.6  2006/11/21 18:09:24  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.5  2006/11/20 15:55:03  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.4  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.3  2006/10/06 14:04:15  cocoa
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
// Revision 1.3  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.2  2006/02/20 22:41:20  cocoa
// All forms of the log function for power products now return SmallExponent_t
// (instead of int).  exponents now resizes the vector rather than requiring
// the user to pass in the correct size.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.5  2005/08/08 16:36:33  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.4  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
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
// Revision 1.8  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from ZZ).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.7  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.6  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.5  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.4  2005/03/11 18:36:36  cocoa
// -- new: NewPolyRingHom Rx-->Ry ; code is still commented out
//
// Revision 1.3  2005/02/11 16:45:24  cocoa
// Removed the useless and misleading functions myInit and myKill
// from the SmallFp*Impl classes; various consequential changes.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.15  2004/12/16 13:35:01  cocoa
// Added new monomial functions (with more coherent signatures).
//
// Revision 1.14  2004/11/25 16:14:21  cocoa
// (1) Fixed definition of specialization of std::swap template function
//     so that it compiles with gcc 3.4.3
// (2) Implemented monomial function for polynomial rings.
// (3) Added one(PPM) and PPM->myOne() functions.
//
// Revision 1.13  2004/11/19 16:15:51  cocoa
// (a) Removed unused error message about degree of zero;
//     replaced it by a more general message about needing a
//     non-zero polynomial (for various ops such as LC, LPP).
// (b) Added some missing arg checking in LC, LPP and deg
//     (for elems of a PolyRing).
// (c) Updated some commented out code in GPair and GPoly.
//
// Revision 1.12  2004/11/19 15:44:27  cocoa
// Changed names of "casting" functions which convert a ring into
// one with a more special structure (e.g. FractionField).  These
// functions now have names starting with "As".  There were several
// consequential changes.
//
// Revision 1.11  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.10  2004/11/08 14:04:47  cocoa
// -- small changes for doxygen support
//
// Revision 1.9  2004/10/29 16:05:03  cocoa
// -- changed PPOrdering::ExpvElem --> SmallExponent_t
//
// Revision 1.8  2004/07/27 16:03:39  cocoa
// Added IsCommutative test and IamCommutative member function
// to all rings.  Tidied geobuckets a little.
//
// Revision 1.7  2004/07/20 15:04:06  cocoa
// The next step in the new "ring element" conversion process:
// handling the case of creating a "const RefRingElem" object
// (since C++ refuses to do this properly itself).
//
// Revision 1.6  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.5  2004/05/24 15:52:14  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.4  2004/01/28 15:27:55  cocoa
// Added function CmpDeg and related member functions myCmpDeg.
//
// Revision 1.3  2003/10/17 10:51:06  cocoa
// Major cleaning, and new naming convention.
//
// Revision 1.2  2003/10/09 12:18:42  cocoa
// New coding convention for rings.
//
// Revision 1.6  2003/06/23 17:04:39  abbott
// Minor cleaning prior to public release.
//
// Revision 1.5  2003/05/30 12:24:03  abbott
// Added CoeffEmbeddingHom.  Added declaration of plain multiplication
// which had become hidden by C++'s "extraordinary" visibility rules.
// ["extraordinary" is the abridged adjective]
//
// Revision 1.4  2003/05/14 16:50:48  abbott
// WIDESPREAD CHANGES, some consequential, and some due to restructuring.
// Updated AbstractRing to RingBase (hmm, should have checked in
// ages ago), and numerous other related and consequential changes.
// Eliminated PolyRing::owner since it cannot work.
// Renamed the old class PolyRing to PolyRingBase, and created a
// new class called PolyRing.
//
// Revision 1.3  2002/11/18 17:59:03  bigatti
// - added  DontSkipLMg  in  SkipLMFlag  enum
// - added  ReductionStep  and  ReductionStepGCD
//
// Revision 1.2  2002/11/15 15:14:47  abbott
// Revised according to the renaming in ring.H.
// Extensively modified and reorganized.  Still in a state of flux.
//
// Revision 1.1  2002/06/22 16:55:27  abbott
// Initial revision
//

#endif
