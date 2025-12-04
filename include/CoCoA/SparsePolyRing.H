#ifndef CoCoA_SparsePolyRing_H
#define CoCoA_SparsePolyRing_H

//   Copyright (c)  2005-2007,2009,2012  John Abbott, and Anna M. Bigatti

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


// Header file for the abstract class SparsePolyRingBase etc.

#include "CoCoA/PPMonoid.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"
#include "CoCoA/ideal.H"
#include "CoCoA/FinGenModule.H"
//#include "CoCoA/MatrixForOrdering.H"
//#include "CoCoA/TmpUniversalInvolutiveBasisContainer.H"

#include <iosfwd>
// using std::ostream;
#include <map>
#include <memory>
// using std::auto_ptr;
#include <vector>
// using std::vector;

namespace CoCoA
{

  class CpuTimeLimit;   // fwd decl -- defined in CpuTimeLimit.H
  class symbol;         // fwd decl -- defined in symbol.H
  class matrix;         // fwd decl -- defined in matrix.H
  class SparsePolyIter; // fwd decl -- defined in SparsePolyIter.H

  namespace Involutive
  {
    class UniversalInvolutiveBasisContainer; // fwd decl -- defined in TmpUniversalInvolutiveBasisContainer.H
  }


  class SparsePolyRingBase;

  class SparsePolyRing: public PolyRing
  {
  public:
    /*implicit*/ SparsePolyRing(const ring& R);
    SparsePolyRing(const ring& R, const ErrorContext& ErrCtx);
    explicit SparsePolyRing(const SparsePolyRingBase* RingPtr);
    // Default copy ctor & assignment work fine.
    // Default dtor works fine.
  public:
    const SparsePolyRingBase* operator->() const; ///< allow member fns to be called
  };


  class SparsePolyRingBase: public PolyRingBase
  {
  public:
    ///@name Member functions every concrete SparsePolyRing implementation must have in addition to those of PolyRingBase.
    //@{
    virtual const PPMonoid& myPPM() const = 0;
    virtual long myGradingDim() const { return GradingDim(myPPM()); }
    virtual bool IamStdGraded() const { return IsStdGraded(myPPM()); }
//    virtual void myWDeg(degree& d, ConstRawPtr rawf) const = 0;
//    virtual int myCmpWDeg(ConstRawPtr rawf, ConstRawPtr rawg) const = 0; //pure???
    virtual ConstRefPPMonoidElem myLPP(ConstRawPtr rawf) const = 0;
    virtual void myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const = 0;
    virtual bool myIsHomog(ConstRawPtr rawf) const; ///< true iff f is homogeneous
    virtual bool myIsHomogPartial(ConstRawPtr rawf, long n) const; ///< true iff f is homogeneous wrt the first n rows of weights (assumes 0 <= n <= GradingDim)
    virtual void myHomog(RawPtr rawfHom, ConstRawPtr rawf, ConstRawPtr rawh) const;  ///< fHom = f homogenized with h (indet of wdeg 1, GrDim=1)
    //@}

    ///@name   Functions for creating/building polynomials
    //@{
    virtual RingElem myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const = 0; // c is element of CoeffRing!!
    virtual SparsePolyIter myBeginIter(ConstRawPtr rawf) const = 0;
    virtual SparsePolyIter myEndIter(ConstRawPtr rawf) const = 0;
    virtual void myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const = 0; /// SHOULD BE vector<BigInt> ????
    virtual void myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const = 0; /// SHOULD BE vector<BigInt> ????
    virtual void myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const = 0;
    virtual void myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const = 0;
    //@}

    ///@name functions inherited from PolyRing with a unique implementation for all SparsePolyRing implementations
    //@{
    long myStdDeg(ConstRawPtr rawf) const override;                     ///< returns standard degree of f
    long myDeg(ConstRawPtr rawf, long index) const override;            ///< returns the degree of f in x[index]
    void myContent(RawPtr rawcontent, ConstRawPtr rawf) const override; ///< exception clean
    void myContentFrF(RawPtr rawcontent, ConstRawPtr rawf) const override; ///< exception clean
    void myCommonDenom(RawPtr rawcontent, ConstRawPtr rawf) const override; ///< exception clean
    void myClearDenom(RawPtr rawres, ConstRawPtr rawf) const override;  ///< exception clean
    void myRemoveBigContent(RawPtr rawf) const override;                ///< exception guarantee: same as myDivByCoeff
    bool myIsValid(ConstRawPtr rawf) const override;                    ///< true iff f is a valid representation, FOR DEBUGGING
    const symbol& myIndetSymbol(long idx) const override  {return myPPM()->myIndetSymbol(idx);}
    //@}

    ///@name functions inherited from ring with a unique implementation for all SparsePolyRing implementations
    //@{
//????    using PolyRingBase::myOutputSelf; // disable warnings of overloading
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;
    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = x*y
    //    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = x/y, if divisible
    void myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const override; ///< lhs = deriv(f, x)
    void mySymbols(std::vector<symbol>& SymList) const override;            ///< append symbols in ring to SymList
    RingElem mySymbolValue(const symbol& sym) const override;               ///< returns the RingElem corresponding to sym
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;      ///< out << x
    bool myIsPrintAtom(ConstRawPtr rawx) const override;                    ///< x^n may be printed without parentheses
    bool myIsPrintedWithMinus(ConstRawPtr rawx) const override;             ///< x^n may be printed without parentheses
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;  ///< OMOut << x
    bool myIsOne(ConstRawPtr rawf) const override;                          ///< true iff rawf == 1
    bool myIsMinusOne(ConstRawPtr rawf) const override;                     ///< true iff rawf == -1
    bool myIsConstant(ConstRawPtr rawf) const override;                     ///< true iff rawf == c (possibly 0)
    bool myIsIndet(long& IndetIndex, ConstRawPtr rawf) const override;      ///< true iff rawf is x[IndetIndex]
    bool myIsIndetPosPower(long& IndetIndex, BigInt& EXP, ConstRawPtr rawf) const /*override*/;                ///< true iff rawf is x[i]^EXP for i=IndexIndex and EXP>0
    bool myIsIndetPosPower(ConstRawPtr rawf) const override;                ///< true iff rawf is x[i]^d for some i and d>0
    bool myIsEvenPoly(ConstRawPtr rawf) const override;
    bool myIsOddPoly(ConstRawPtr rawf) const override;
    bool myIsInvertible(ConstRawPtr rawx) const override;                   ///< true iff x is invertible
    bool myIsZeroDivisor(ConstRawPtr rawx) const override;
    void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override; ///< lhs = gcd(x,y) if TrueGCDDomain;
    void myNormalizeFracNoGcd(RawPtr rawnum, RawPtr rawden) const override; ///< normalize fraction num/den
    bool myIsInteger(BigInt& N, ConstRawPtr rawx) const override;           ///< true iff x is integer
    bool myIsRational(BigRat& Q, ConstRawPtr rawx) const override;          ///< true iff x is rational
    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;
    void mySquare(RawPtr rawlhs, ConstRawPtr rawx) const override;          ///< lhs = x^2
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override; ///< lhs = x^n, n>1, x not -1,0,1
    //@}

    ///@name  Functions for homomorphisms
    //@{
    RingHom myCoeffEmbeddingHomCtor() const override;
    RingHom myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const override;
    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override;

    bool myImageLiesInSubfield(const RingHom& phi) const override;
    //@}

    ///@name Special functions on polynomials needed for implementing Buchberger's Algorithm
    //@{
    virtual bool myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const = 0; ///< f+=LM(g); g-=LM(g); assumes LPP(f)==LPP(g); returns LC(f)+LC(g)==0
    virtual void myMoveLMToFront(RawPtr rawf, RawPtr rawg) const = 0; ///< f+=LM(g); g-=LM(g); assumes LM(f)<LM(g)
    virtual void myMoveLMToBack(RawPtr rawf, RawPtr rawg) const = 0; ///< f+=LM(g); g-=LM(g); assumes last PP of f>LM(g)
    virtual void myDeleteLM(RawPtr rawf) const = 0;            // ????? right interface
    virtual void myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const = 0; ///< lhs=div(LM(f),LM(g)); assumes f!=0,g!=0
    virtual int  myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const = 0; ///< cmp(LPP(f),LPP(g)); assumes f!=0,g!=0
    virtual void myAddClear(RawPtr rawf, RawPtr rawg) const = 0;     ///< f+=g; g=0;
    virtual void myAppendClear(RawPtr rawf, RawPtr rawg) const = 0;  ///< f+=g; g=0; appends g to f with no checks

    enum SkipLMFlag { SkipLMg, DontSkipLMg }; // used only for AddMul
    virtual void myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const = 0; ///<  f += LM(h)*g
    virtual void myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag) const = 0; ///<  f += LM(h)*g
    virtual void myReductionStep(RawPtr rawf, ConstRawPtr rawg) const = 0;
    // ??? aggiungere coefficiente
    virtual void myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& FScale) const = 0;
    // should it all be in ReductionStep ??? ANNA
    //@}

// need to be friend for SparsePolyRingBase::IdealImpl::ourGetPtr(I);
// not nice, not nice... (AMB 2011-04-19)
//     friend bool HasGBasis(const ideal& I);
//     friend bool AreGensMonomial(const ideal& I);
//     friend bool AreGensSqFreeMonomial(const ideal& I);
//     friend void SetGBasisAsGens(const ideal& I);
//     friend const std::vector<RingElem>& MinGens(const ideal& I);
//     friend std::vector<RingElem> GBasisSelfSatCore(const ideal& I);
//     friend std::vector<RingElem> GBasisRealSolve(const ideal& I);


    //  protected: // Ideal class for this type of ring.
    //  now public!! is this OK, or just a workaround (ourGetPtr)? //AMB 2018-03

    //------------------------------------------------------------------
    //--- SparsePolyRingBase::IdealImpl --------------------------------

    class IdealImpl: public IdealBase
    {
      friend std::vector<RingElem> GBasisTrunc(const ideal& I, const long TruncDeg, const CpuTimeLimit&);
      friend const std::vector<RingElem>& GBasisByHomog(const ideal& I, const CpuTimeLimit&);
      friend std::vector<RingElem> GBasisSelfSatCore(const ideal& I);
      friend std::vector<RingElem> GBasisRealSolve(const ideal& I);
      friend const std::vector<RingElem>& MinGens(const ideal& I);
    public:
      IdealImpl(const SparsePolyRing& P, const std::vector<RingElem>& gens);
      // default copy ctor is OK
      IdealBase* myClone() const override;
//???    virtual void swap(ideal& other);

    public: // functions every ideal must implement
      const SparsePolyRing& myRing() const override;
      bool IamZero() const override;
      // bool IamOne() const override;  // default impl
      virtual bool IhaveMonomialGens() const;
      virtual bool IhaveSqFreeMonomialGens() const;
      virtual bool IhaveGBasis() const;
      virtual void mySetGBasisAsGens() const;
      virtual void myReset() const; // reset all member fields but myGens
      void myReduceMod(RingElemRawPtr rawx) const override; // r elem of R, where I is ideal in R
      bool IhaveElem(RingElemConstRawPtr rawx) const override;
      virtual bool IamZeroDim() const;
      void myAdd(const ideal&) override;
      void myMul(const ideal&) override;
      void myIntersect(const ideal&) override;
      void myColon(const ideal&) override;
      void mySaturate(const ideal&) override;
      void myMinimalize() override;
      //      ideal myElim(const std::vector<RingElem>& ELimIndets) override;
      void myAssignElim(const ideal& I, const std::vector<RingElem>& ELimIndets) override;
      bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const override; // lhs = num/den modulo the ideal  (lhs = 0 if quotient does not exist)
      virtual const std::vector<RingElem>& myJanetBasis() const;

      const std::vector<RingElem>& myGens() const override; // gens as specified by user
      const std::vector<RingElem>& myTidyGens(const CpuTimeLimit&) const override; // tidier set of gens


      //--- involutive --------------------------------
      // all the involutive bases stuff always with prefix Inv to avoid some naming collisions
      virtual bool InvIamDeltaRegular() const;
      virtual bool InvIamMonomial() const;
      virtual bool InvIamHomogeneous() const;
      virtual bool InvIamCohenMacaulay() const;
      virtual void InvMultVars(std::map<PPMonoidElem, std::vector<bool> >& MultVars) const;
      virtual void InvNonMultVars(std::map<PPMonoidElem, std::vector<bool> >& NonMultVars) const;
      virtual void InvHilbertPol(RingElem& HilbertPol, ConstRefRingElem var) const;
      virtual void InvHilbertSeries(RingElem& HilbertSeries, ConstRefRingElem var) const;
      virtual void InvFirstSyzygy(FinGenModule& syz) const;
      virtual long InvDimension() const;
      virtual void InvComplementaryDecomposition(std::vector<std::pair<PPMonoidElem, std::vector<bool> > >& CompDecomp) const;
      virtual long InvDepth() const;
      virtual long InvProjDim() const;
      virtual void InvSocle(std::vector<RingElem>& socle) const;
      virtual void InvExtremalBettiNumbers(std::map<std::pair<long, long>, long>& ExtremalBettiNumbers) const;
      virtual void InvRegularSequence(std::vector<RingElem>& RegSeq) const;
      virtual void InvMaximalStronglyIndependentSet(std::vector<RingElem>& MaxSet) const;
      virtual long InvRegularity() const;
      virtual long InvSatiety() const;
      virtual void InvSaturation(std::vector<RingElem>& saturation) const;

      //---- data members -----------------------------------------------
    protected: // Data members of SparsePolyRingBase::IdealImpl
      SparsePolyRing myP;
      std::vector<RingElem> myGensValue;
      SmartPtrIRC<Involutive::UniversalInvolutiveBasisContainer> myInvBasisContainerPtr;
      mutable bool3 IhaveMonomialGens3Value; ///< the given gens are monomials
      mutable bool3 IhaveSqFreeMonomialGens3Value; ///< the given gens are square free monomials
      mutable bool IhaveGBasisValue;
      mutable std::vector<RingElem> myGBasisValue;
      mutable std::vector<RingElem> myMinGensValue;
      //      mutable HPSeries myHPSeriesValue;

    protected: // Auxiliary functions.
      virtual void myGBasis_EasyCases() const;
      virtual const std::vector<RingElem>& myGBasis(const CpuTimeLimit&) const;
      virtual std::vector<RingElem> myGBasisTrunc(const long TruncDeg, const CpuTimeLimit&) const;
      virtual std::vector<RingElem> myGBasisTrunc_compute(const long TruncDeg, const CpuTimeLimit&) const;
      virtual const std::vector<RingElem>& myGBasisByHomog(const CpuTimeLimit&) const;
      virtual std::vector<RingElem> myGBasisSelfSatCore() const;
      virtual std::vector<RingElem> myGBasisRealSolve() const;
      virtual const std::vector<RingElem>& myMinGens() const;
      void myTestIsMaximal() const override;
      void myTestIsPrimary() const override;
      void myTestIsPrime() const override;
      void myTestIsRadical() const override;
      std::vector<ideal> myPrimaryDecomposition() const override;

    private:
      //-- 0 dimensional ideals ------------------
      void myTestIsMaximal_0dim() const;
      void myTestIsPrimary_0dim() const;
      void myTestIsRadical_0dim() const;
      ideal myRadical_0dimDRL() const;
      std::vector<ideal> myPrimaryDecomposition_0dim() const;
      void myPrimaryDecompositionCore_0dim(bool& IsCertified, std::vector<ideal>& Q) const;

      //--- monomial ideals ----------------------
      const std::vector<RingElem>& myGBasis_MonId() const;
      const std::vector<RingElem>& myMinGens_MonId() const;
      void myTestIsRadical_MonId() const;
      ideal myRadical_MonId() const;
      std::vector<ideal> myPrimaryDecomposition_MonId() const;
      void myMul_MonId(const ideal&);
      void myIntersect_MonId(const ideal&);
      void myColon_MonId(const ideal&);
      //      void myElim_MonId(const std::vector<RingElem>& ELimIndets);
      void myAssignElim_MonId(const ideal& I, const std::vector<RingElem>& ELimIndets);

    public:
      static const SparsePolyRingBase::IdealImpl* ourGetPtr(const ideal& J);

    // need to be friend for SparsePolyRingBase::IdealImpl::myTestIsPrimary_0dim
    friend bool IsPrimary_0dim(const ideal& I);
    // need to be friend for SparsePolyRingBase::IdealImpl::myRadical_0dimDRL
    friend ideal radical_0dimDRL(const ideal& I);
    // need to be friend for SparsePolyRingBase::IdealImpl::myRadical_MonId
    friend ideal radical_MonId(const ideal& I);
    };  // class IdealImpl


    //----- HomImpl ---------------------------------------------------
  protected:  // General homomorphism class for this type of ring.

    class HomImpl: public RingHomBase
    {
    public:
      HomImpl(const SparsePolyRing& domain, const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages);
      void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const override;
      bool IamPartial() const override  { return IsPartial(myCoeffHom); }
    private:
      void myOutputSelfDetails(std::ostream& out) const override;

    private: // Data members of SparsePolyRingBase::HomImpl
      const RingHom myCoeffHom;
      const std::vector<RingElem> myIndetImages;
    };  // class HomImpl

  protected:  // Special homomorphism class for this type of ring.

    class CoeffEmbeddingHomImpl: public RingHomEmbeddingBase
    {
    public:
      CoeffEmbeddingHomImpl(const SparsePolyRing& P);
      virtual void myApply(RingElemRawPtr rawimage, RingElemConstRawPtr rawarg) const;
    };  // class CoeffEmbeddingHomImpl

  };


  SparsePolyRing NewPolyRing(const ring& CoeffRing, const PPMonoid& PPM);
  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);
  //  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetNames, ConstMatrixView M, long GrDim);
  SparsePolyRing NewPolyRing(const ring& CoeffRing, const std::vector<symbol>& IndetNames); ///< Default ordering is StdDegRevLex
  SparsePolyRing NewPolyRingWeights(const ring& CoeffRing, const std::vector<symbol>& IndetNames, ConstMatrixView Ws);

  bool IsSparsePolyRing(const ring& R);                  // inline
  const PPMonoid& PPM(const SparsePolyRing& Rx);         // inline
  long GradingDim(const SparsePolyRing& Rx);             // inline
  bool IsStdGraded(const SparsePolyRing& Rx);            // inline
  ConstMatrixView OrdMat(const SparsePolyRing& Rx);
  ConstMatrixView GradingMat(const SparsePolyRing& Rx);
  bool HasPositiveGrading(const SparsePolyRing& Rx);

  bool AreGoodIndetNames(const ring& CoeffRing, const std::vector<symbol>& IndetNames);


  //----------------------------------------------------------------------
  // Inline functions

  inline const SparsePolyRingBase* SparsePolyRingPtr(const ring& R)
  { return dynamic_cast<const SparsePolyRingBase*>(R.myRawPtr()); }


  inline const SparsePolyRingBase* SparsePolyRingPtr(const ring& R, const ErrorContext& ErrCtx)
  {
    const SparsePolyRingBase* ptr = SparsePolyRingPtr(R);
    if (ptr == nullptr)  CoCoA_THROW_ERROR_WITH_CONTEXT2(ERR::ReqSparsePolyRing, ErrCtx);
    return ptr;
  }


  inline bool IsSparsePolyRing(const ring& R)
  { return SparsePolyRingPtr(R) != nullptr; }


  inline SparsePolyRing::SparsePolyRing(const ring& R):
      PolyRing(SparsePolyRingPtr(R, CoCoA_ERROR_CONTEXT))
  {}

  inline SparsePolyRing::SparsePolyRing(const ring& R, const ErrorContext& ErrCtx):
      PolyRing(SparsePolyRingPtr(R, ErrCtx))
  {}

  inline SparsePolyRing::SparsePolyRing(const SparsePolyRingBase* RingPtr):
      PolyRing(RingPtr)
  {}


  inline const SparsePolyRingBase* SparsePolyRing::operator->() const
  { return static_cast<const SparsePolyRingBase*>(ring::operator->()); }

  inline const PPMonoid& PPM(const SparsePolyRing& Rx)
  { return Rx->myPPM(); }

  inline long GradingDim(const SparsePolyRing& Rx)
  { return Rx->myGradingDim(); }

  inline bool IsStdGraded(const SparsePolyRing& Rx)
  { return Rx->IamStdGraded(); }

  inline bool HasLex(const SparsePolyRing& P)
  { return IsLex(ordering(PPM(P))); }
  
  inline bool HasStdDegLex(const SparsePolyRing& P)
  { return IsStdDegLex(ordering(PPM(P))); }

  inline bool HasStdDegRevLex(const SparsePolyRing& P)
  { return IsStdDegRevLex(ordering(PPM(P))); }

  inline bool HasMatrixOrdering(const SparsePolyRing& P)
  { return IsMatrixOrdering(ordering(PPM(P))); }

} // end of namespace CoCoA

#endif
