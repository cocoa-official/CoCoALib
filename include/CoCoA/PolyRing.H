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
    const ring& myBaseRing() const override final  { return myCoeffRing(); }
    BigInt myCharacteristic() const override final; ////  { return characteristic(myCoeffRing()); }
    bool IamCommutative() const override  { return true; } ///< RingWeyl redefines this
    bool3 IamIntegralDomain3(bool QuickMode) const override final  { return myCoeffRing()->IamIntegralDomain3(QuickMode); }
    bool IamTrueGCDDomain() const override /*final*/  { return myCoeffRing()->IamTrueGCDDomain() || myCoeffRing()->IamField(); };
    bool IamField() const override final  { return false; } //??? (myNumIndets()==0 && myCoeffRing()->IamField())
    bool IamFiniteField() const override final  { return false; } //??? (myNumIndets()==0 && myCoeffRing()->IamFiniteField())
    bool IamExact() const override final  { return IsExact(myCoeffRing()); }
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

#endif
