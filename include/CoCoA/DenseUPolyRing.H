#ifndef CoCoA_DenseUPolyRing_H
#define CoCoA_DenseUPolyRing_H

//   Copyright (c)  2007  Anna Bigatti

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


// Header file for the abstract class DenseUPolyRingBase etc.


#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"

#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;
#include <memory>
// using std::auto_ptr;


namespace CoCoA
{

  class symbol;  // fwd decl -- defined in symbol.H

  class DenseUPolyRingBase;


  class DenseUPolyRing: public PolyRing
  {
  public:
    /*implicit*/ DenseUPolyRing(const ring& R);
    DenseUPolyRing(const ring& R, const ErrorContext& ErrCtx);
    explicit DenseUPolyRing(const DenseUPolyRingBase* RingPtr);
    // Default copy ctor & assignment work fine.
    // Default dtor works fine.
  public:
    const DenseUPolyRingBase* operator->() const; ///< allow member fns to be called
  };


  class DenseUPolyRingBase: public PolyRingBase
  {
  public:
    ///@name  Member functions every concrete DenseUPolyRing implementation must have in addition to those of PolyRingBase.
    //@{
    using PolyRingBase::myIndetSymbol;    // disable warnings of overloading
    virtual const symbol& myIndetSymbol() const = 0;
    virtual std::string myImplDetails() const override = 0;
    virtual long mySize(ConstRawPtr rawf) const = 0;            ///< returns (memory) size of f
    virtual long myDegPlus1(ConstRawPtr rawf) const = 0;        ///< returns deg of f +1 (0 for zero poly)
    virtual RingElemAlias myCoeff(ConstRawPtr rawf, long d) const = 0;
    virtual void myAddMulLM(RawPtr rawf, ConstRawPtr rawc, long d, ConstRawPtr rawg) const = 0; ///<  f += c*indet^d*g
    virtual void myMulByXExp(RawPtr rawf, unsigned long d) const = 0;
    virtual void myMulBy1MinusXExp(RawPtr rawf, unsigned long d) const = 0;
    //@}

    ///@name  Functions for creating/building polynomials
    //@{
    virtual RingElem myMonomial(ConstRawPtr rawc, unsigned long d) const; ///< c is element of CoeffRing!!
    virtual void myAssignZeroCoeff(RawPtr rawf, long d) const = 0; ///< f_d = 0, no check on size
    virtual void myAssignNonZeroCoeff(RawPtr rawf, ConstRawPtr rawc, long d) const = 0; ///< f_d = c, no check on size
    virtual void myResize(RawPtr rawf, long NewSize) const = 0;
    virtual void myResetDeg(RawPtr rawf) const = 0;
    //@}

    ///@name  Functions inherited from PolyRing with a unique implementation for all DenseUPolyRing implementations
    //@{
    long myNumIndets() const override  {return 1;}
    void myIndetPower(RawPtr rawf, long var, long exp) const override;
    long myNumTerms(ConstRawPtr rawf) const override;
    bool myIsConstant(ConstRawPtr rawf) const override;               ///< true iff rawf == c (possibly 0)
    bool myIsIndet(long& IndetIndex, ConstRawPtr rawf) const override;///< true iff rawf is x[IndetIndex]
    bool myIsIndetPosPower(ConstRawPtr rawf) const override;          ///< true iff rawf is x[i]^d for some i and d>0
    bool myIsMonomial(ConstRawPtr rawf) const override;               ///< is it useful for univariate?
    bool myIsEvenPoly(ConstRawPtr rawf) const override;
    bool myIsOddPoly(ConstRawPtr rawf) const override;

    long myStdDeg(ConstRawPtr rawf) const override;                   ///< standard degree of f, error for zero poly
    long myDeg(ConstRawPtr rawf, long index) const override;          ///< max exponent of index (same as myStdDeg)
    RingElemAlias myLC(ConstRawPtr rawf) const override;
    void myContent(RawPtr rawcontent, ConstRawPtr rawf) const override;
    void myContentFrF(RawPtr rawcontent, ConstRawPtr rawf) const override;
    void myCommonDenom(RawPtr rawcontent, ConstRawPtr rawf) const override;
    void myClearDenom(RawPtr rawcontent, ConstRawPtr rawf) const override;
    void myRemoveBigContent(RawPtr rawf) const override;
    bool myIsValid(ConstRawPtr rawf) const override;                    ///< true iff f is a valid representation, for debugging
    //    bool myIsHomog(ConstRawPtr rawf) const override;                  ///< true iff f is homogeneous
    const symbol& myIndetSymbol(long idx) const override;
    //@}

    ///@name  Functions inherited from ring with a single implementation for all DenseUPolyRing implementations
    //@{
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = x*y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                ///< true iff x is invertible constant
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = gcd(x,y) if TrueGCDDomain;
    void myExgcd(RawPtr rawlhs, RawPtr rawxcofac, RawPtr rawycofac, ConstRawPtr rawx, ConstRawPtr rawy) const; ///< lhs = gcd(x,y) = xcofac*x + ycofac*y  if TrueGCDDomain;
    void mySymbols(std::vector<symbol>& SymList) const override;         ///< append symbols in ring to SymList
    RingElem mySymbolValue(const symbol& sym) const override;            ///< returns the RingElem corresponding to sym
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;   ///< out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;                 ///< <indet>^n may be printed without parentheses
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;///< OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                       ///< x == 0
    bool myIsOne(ConstRawPtr rawf) const override;                        ///< true iff rawf == 1
    bool myIsMinusOne(ConstRawPtr rawf) const override;                   ///< true iff rawf == -1
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;         ///< true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;        ///< true iff x is rational
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;    ///< true iff rawx == rawy
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;///< lhs = x^n, n>1, x not -1,0,1
    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;
    //@}

    ///@name  Functions for homomorphisms
    //@{
    RingHom myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const override;
    //    RingHom myHomCtor(const ring& codomain, const RingHom& CoeffHom, ConstRefRingElem IndetImage) const;
    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override;
    //@}

  protected: // Ideal class for this type of ring.

    class IdealImpl: public IdealBase
    {
    public:
      IdealImpl(const DenseUPolyRing& P, const std::vector<RingElem>& gens);
      // default copy ctor is OK
      IdealBase* myClone() const override;
//???    virtual void swap(ideal& other);

      const DenseUPolyRing& myRing() const override { return myP; }
      bool IamZero() const override;
      // virtual bool IamOne() const;  // default impl
      void myReduceMod(RingElemRawPtr rawx) const override; ///< x elem of P, where I is ideal in P
      bool IhaveElem(RingElemConstRawPtr rawx) const override;
      void myAdd(const ideal&) override;
      void myMul(const ideal&) override;
      void myIntersect(const ideal&) override;
      void myColon(const ideal&) override;
      bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const override; ///< ??? lhs = num/den modulo the ideal  (lhs = 0 if quotient does not exist)

      const std::vector<RingElem>& myGens() const override { return myGensValue; } ///< gens as specified by user
      const std::vector<RingElem>& myTidyGens(const CpuTimeLimit& CheckForTimeout) const override; ///< tidier set of gens
    protected:
      void myTestIsMaximal() const override;
      void myTestIsPrimary() const override;
      void myTestIsPrime() const override;
      void myTestIsRadical() const override;

    protected: // Data members of DenseUPolyRingBase::IdealImpl
      DenseUPolyRing myP;
      std::vector<RingElem> myGensValue;
      mutable bool myTidyGensIsValid;
      mutable std::vector<RingElem> myTidyGensValue;
    protected: // Auxiliary functions.
      static const DenseUPolyRingBase::IdealImpl* GetPtr(const ideal& J);
    };  // class IdealImpl

  protected:    //  General homomorphism class for this type of ring.
    class HomImpl: public RingHomBase
    {
    public:
      HomImpl(const DenseUPolyRing& domain, const ring& codomain, const RingHom& CoeffHom, ConstRefRingElem IndetImage);
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return IsPartial(myCoeffHom); }
    private:
      void myOutputSelfDetails(std::ostream& out) const override;

    private: // Data members of DenseUPolyRingBase::HomImpl
      const RingHom myCoeffHom;
      RingElem myIndetImage;
    };  // class HomImpl

  protected:    //  Special homomorphism class for this type of ring.
    class CoeffEmbeddingHomImpl: public RingHomEmbeddingBase
    {
    public:
      CoeffEmbeddingHomImpl(const DenseUPolyRing& P);
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
    };  // class CoeffEmbeddingHomImpl
  };


  bool IsDenseUPolyRing(const ring& R);                         // inline
  bool IsGoodIndetName(const ring& CoeffRing, const symbol& IndetName);

  ///@name  Pseudo constructors for RingElems owned by DenseUPolyRing
  //@{
  RingElem monomial(const DenseUPolyRing& P, ConstRefRingElem c, const MachineInt& exp);      ///< c*indet^exp
  RingElem monomial(const DenseUPolyRing& P, const BigRat& q, const MachineInt& exp);             ///< q*indet^exp
  RingElem monomial(const DenseUPolyRing& P, const BigInt& N, const MachineInt& exp);         ///< N*indet^exp
  RingElem monomial(const DenseUPolyRing& P, const MachineInt& n, const MachineInt& exp); ///< n*indet^exp
  //@}

  ///@name  Special functions on RingElems owned by DenseUPolyRing
  //@{
  RingElemAlias coeff(ConstRefRingElem f, long d);
  //@}

  //----------------------------------------------------------------------
  // Inline functions

  inline const DenseUPolyRingBase* DenseUPolyRingPtr(const ring& R)
  {
    return dynamic_cast<const DenseUPolyRingBase*>(R.myRawPtr());
  }

  // inline const DenseUPolyRingBase* DenseUPolyRingPtr(const ring& R, const char* const FnName)
  // {
  //   const DenseUPolyRingBase* ptr = DenseUPolyRingPtr(R);
  //   if (ptr == nullptr) CoCoA_THROW_ERROR(ERR::NotDenseUPolyRing, FnName);
  //   return ptr;
  // }

  inline const DenseUPolyRingBase* DenseUPolyRingPtr(const ring& R, const ErrorContext& ErrCtx)
  {
    const DenseUPolyRingBase* ptr = DenseUPolyRingPtr(R);
    if (ptr == nullptr)
      CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::NotDenseUPolyRing, ErrCtx);
    return ptr;
  }

  inline bool IsDenseUPolyRing(const ring& R)
  {
    return DenseUPolyRingPtr(R) != nullptr;
  }


  inline DenseUPolyRing::DenseUPolyRing(const ring& R):
      PolyRing(DenseUPolyRingPtr(R, CoCoA_ERROR_CONTEXT))
  {}

  inline DenseUPolyRing::DenseUPolyRing(const ring& R, const ErrorContext& ErrCtx):
      PolyRing(DenseUPolyRingPtr(R, ErrCtx))
  {}

  inline DenseUPolyRing::DenseUPolyRing(const DenseUPolyRingBase* RingPtr):
    PolyRing(RingPtr)
  {}


  inline const DenseUPolyRingBase* DenseUPolyRing::operator->() const
  {
    return static_cast<const DenseUPolyRingBase*>(ring::operator->());
  }
  

} // end of namespace CoCoA

#endif
