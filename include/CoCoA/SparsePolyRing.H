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
#include "CoCoA/FGModule.H"
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
      virtual void InvFirstSyzygy(FGModule& syz) const;
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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyRing.H,v 1.145 2024/08/02 13:28:15 bigatti Exp $
// $Log: SparsePolyRing.H,v $
// Revision 1.145  2024/08/02 13:28:15  bigatti
// Summary: just spacing
//
// Revision 1.144  2024/08/01 21:22:41  abbott
// Summary: Fixed so that it compiles: changed name os some error IDs
//
// Revision 1.143  2024/07/23 13:22:35  bigatti
// Summary: added myGBasisTrunc_compute and myMinGens_MonId
//
// Revision 1.142  2024/07/22 15:43:57  abbott
// Summary: Renamed macros CoCoA_THROW_ERROR_WITH_CONTEXT so that last char is num of args
//
// Revision 1.141  2024/05/30 09:41:21  bigatti
// Summary: removed MinGens_almost (redmine #1740)
//
// Revision 1.140  2024/05/28 12:53:40  bigatti
// Summary: added GBasisTrunc, myGBasisTrunc, ComputeGBasisTrunc
//
// Revision 1.139  2024/05/24 14:38:41  bigatti
// Summary: added MinGens_almost
//
// Revision 1.138  2024/03/16 10:15:11  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.137  2024/03/15 19:04:45  abbott
// Summary: Arg checking fns now expect a ErrorContext object (instead of string; related to redmine 1726)
//
// Revision 1.136  2024/03/15 15:35:06  bigatti
// Summary: changed myElim and myElim_MonId into myAssignElim and myAssignElim_MonId
//
// Revision 1.135  2023/07/06 16:32:34  bigatti
// Summary: added NewPolyRing with weight matrix
//
// Revision 1.134  2022/02/18 14:11:49  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.133  2022/02/08 20:18:53  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.132  2021/10/30 19:03:08  abbott
// Summary: Added some more "override" (clang said they were mising; redmine 1625)
//
// Revision 1.131  2021/10/30 17:26:35  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.130  2021/01/07 14:44:46  abbott
// Summary: Corrected copyright
//
// Revision 1.129  2020/06/17 15:49:20  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.128  2020/02/12 09:01:47  bigatti
// -- changed myTestIsMaximal etc to return void (and consequences)
//
// Revision 1.127  2019/12/28 17:55:49  abbott
// Summary: Added myIsIndetPosPower (with 3 args)
//
// Revision 1.126  2019/10/03 13:33:56  bigatti
// -- implemented radical for monomial ideals (and used where useful)
//
// Revision 1.125  2019/09/25 14:35:05  bigatti
// -- removed old commented out code about RingElem
//
// Revision 1.124  2019/09/25 14:28:30  bigatti
// -- added HasPositiveGrading
//
// Revision 1.123  2019/09/16 17:41:41  abbott
// Summary: Added myIsEvenPoly, myIsOddPoly
//
// Revision 1.122  2019/03/18 17:13:25  abbott
// Summary: Used nullptr instead of 0
//
// Revision 1.121  2018/10/02 09:45:30  abbott
// Summary: Moved pseudo-ctors NewPolyRing(CoeffRing, NumIndets, ...) to obsolescent
//
// Revision 1.120  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.119  2018/08/06 09:38:29  bigatti
// -- renamed GBasisViaHomog --> GBasisByHomog
//
// Revision 1.118  2018/08/06 08:57:48  bigatti
// -- added timeout for GBasisViaHomog
// -- now using GBasisViaHomog in IsPrimary and radical
//
// Revision 1.117  2018/06/27 12:15:18  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.116  2018/06/15 08:48:46  abbott
// Summary: Added NewPolyRing with args CoeffRing, NumIndets, PPord
//
// Revision 1.115  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.114  2018/05/18 16:34:52  bigatti
// -- split SparsePolyOps-RingElem from SparsePolyRing
//
// Revision 1.113  2018/05/17 15:22:21  bigatti
// -- SparsePolyIter separated from SparsePolyRing
//
// Revision 1.112  2018/04/16 21:46:16  bigatti
// -- added IamZeroDim
// -- added myPrimaryDecomposition_0dim
//
// Revision 1.111  2018/04/10 14:51:44  bigatti
// -- added virtual myPrimaryDecomposition (with default implementation)
//
// Revision 1.110  2018/04/10 14:20:10  bigatti
// -- started work on primary decomposition
//
// Revision 1.109  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.108  2018/03/20 11:39:56  bigatti
// -- changed ***MonId --> ***_MonId
//
// Revision 1.107  2018/03/20 11:20:44  bigatti
// -- changes names IamRadicalOdimDRL --> myTestIsRadical_0dimDRL (etc)
//
// Revision 1.106  2018/03/15 14:21:55  bigatti
// -- new files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.105  2017/11/29 17:40:59  bigatti
// -- added GBasisRealSolveCore
// -- changed GBasisSelfSatCore: now DOES NOT store GBasis
//
// Revision 1.104  2017/11/24 17:46:40  bigatti
// -- renamed GBasisSelfSat --> GBasisSelfSatCore
// -- added GBasisSelfSat in cpkg5
//
// Revision 1.103  2017/11/23 12:32:28  bigatti
// -- added GBasisSelfSat
//
// Revision 1.102  2017/09/25 12:36:58  abbott
// Summary: Added QuotientBasisSorted
//
// Revision 1.101  2017/07/23 15:31:33  abbott
// Summary: Added GBasisTimeout (just for ideals)
//
// Revision 1.100  2017/07/03 19:52:37  abbott
// Summary: Added CutLF
//
// Revision 1.99  2017/04/26 15:56:00  bigatti
// -- added IdealOfGBasis, IdealOfMinGens
//
// Revision 1.98  2017/03/08 15:33:15  bigatti
// -- just some comments on "sections" (removed doxygen)
//
// Revision 1.97  2016/11/18 18:08:29  abbott
// Summary: Reordered prototypes for some monomial fns
//
// Revision 1.96  2016/11/07 13:54:21  bigatti
// ++ added AreGensSqFreeMonomial
// ++ replaced all SquareFree into SqFree
//
// Revision 1.95  2016/11/07 12:11:21  bigatti
// -- added HasGBasis / IhaveGBasis
//
// Revision 1.94  2016/10/27 13:55:47  bigatti
// -- added myRadicalMonId
//
// Revision 1.93  2016/09/08 14:01:28  bigatti
// -- added test functions for maximal and radical
//
// Revision 1.92  2016/06/29 13:13:49  abbott
// Summary: Changed/corrected rtn type of SparsePolyRingPtr
//
// Revision 1.91  2016/06/24 14:27:41  bigatti
// -- renamed CRT_poly --> CRTPoly
//
// Revision 1.90  2016/04/27 06:59:53  bigatti
// -- new: myReset()
// -- removed: myClearGBasis()
//
// Revision 1.89  2016/04/07 14:35:45  abbott
// Summary: Removed code which disabled assignment of rings
//
// Revision 1.88  2015/12/08 13:56:09  abbott
// Summary: Updated Mario's code!  Very many changes!
//
// Revision 1.87  2015/11/30 21:53:55  abbott
// Summary: Major update to matrices for orderings (not yet complete, some tests fail)
//
// Revision 1.86  2015/10/01 10:13:10  bigatti
// -- added CRT_poly, RatReconstructPoly
//
// Revision 1.85  2015/06/30 12:51:10  abbott
// Summary: Added new fn IndetsIn
// Author: JAA
//
// Revision 1.84  2015/06/11 16:55:16  bigatti
// -- added monomial(ring, pp) and monomial(ring, expv)
//
// Revision 1.83  2015/04/30 08:44:18  bigatti
// -- added myIsZeroDivisor
// -- removed some newlines in function definition
//
// Revision 1.82  2015/04/24 15:40:59  bigatti
// -- renamed: myAddMul --> myAddMulLM
// -- renamed: myMoveLM --> myMoveLMToFront
// -- new myMoveLMToBack (used in ReductionCog --> bug in test-TmpMorseGraph??)
//
// Revision 1.81  2014/07/31 13:10:46  bigatti
// -- GetMatrix(PPO) --> OrdMat(PPO)
// -- added OrdMat and GradingMat to PPOrdering, PPMonoid, SparsePolyRing
//
// Revision 1.80  2014/07/30 14:01:14  abbott
// Summary: Changed myAmbientRing into myRing
// Author: JAA
//
// Revision 1.79  2014/07/28 15:41:47  abbott
// Summary: Added myCoeffEmbeddingHomCtor
// Author: JAA
//
// Revision 1.78  2014/07/11 15:42:54  bigatti
// -- added  myOutputSelf (default impl),  myImplDetails()
//
// Revision 1.77  2014/07/07 17:11:57  abbott
// Summary: [MAJOR CHANGE] Removed AsSparsePolyRing; added SparsePolyRingPtr
// Author: JAA
//
// Revision 1.76  2014/05/06 13:20:41  abbott
// Summary: Changed names (my)MaxExponents into (my)Deg
// Author: JAA
//
// Revision 1.75  2014/04/15 13:28:55  abbott
// Summary: Added new ClearDenom fn (with 2 args)
// Author: JAA
//
// Revision 1.74  2014/03/27 14:57:23  bigatti
// -- added myMinimalize
//
// Revision 1.73  2014/01/28 13:08:32  bigatti
// -- added LF for ideal
//
// Revision 1.72  2013/07/30 15:01:51  bigatti
// -- added IsStdGraded
// -- removed empty lines in inline functions
//
// Revision 1.71  2013/06/28 17:03:51  abbott
// Modified semantics of IdealBase::myDivMod;
// it now returns a boolean.
// Several consequential changes.
//
// Revision 1.70  2013/06/12 08:49:17  bigatti
// -- added MinGens
//
// Revision 1.69  2013/05/30 13:13:58  bigatti
// -- added SetGBasisAsGens
//
// Revision 1.68  2013/05/28 13:29:54  abbott
// Added printing for struct CoeffPP.
//
// Revision 1.67  2013/02/21 12:50:44  abbott
// Added new fn UnivariateIndetIndex.
//
// Revision 1.66  2013/01/21 14:07:08  bigatti
// -- commented out "IamPartial"  (for the time being)
//
// Revision 1.65  2013/01/17 15:11:19  abbott
// Added new fn CoeffVecWRT.
// Added new virt mem fns myImageLiesInSubfield & IamPartial.
//
// Revision 1.64  2012/10/24 13:31:08  abbott
// Changed return type of coeff/myCoeff.
// Minor change to signature of CoefficientsWRT & ContentWRT.
//
// Revision 1.63  2012/10/17 09:40:16  abbott
// Replaced  RefRingElem  by  RingElem&
// (plus a few consequential changes)
//
// Revision 1.62  2012/10/05 10:21:39  bigatti
// -- added LF (leading form)
//
// Revision 1.61  2012/10/02 16:44:04  bigatti
// -- added homog for ideal
//
// Revision 1.60  2012/05/30 13:44:11  bigatti
// -- renamed IhaveMonomialGensB3Value --> IhaveMonomialGens3Value
//
// Revision 1.59  2012/05/22 10:02:37  abbott
// Removed IsGCDDomain; substituted by IsTrueGCDDomain.
// Added IsFractionFieldOfGCDDomain.
//
// Revision 1.58  2012/04/03 16:12:46  abbott
// Changed slightly the fn signatures of CoefficientsWRT.
// Added some things to *.C file -- will complete later.
//
// Revision 1.57  2012/03/02 14:21:46  bigatti
// -- added ContentWRT, CoefficientsWRT(f, x)
//
// Revision 1.56  2012/02/10 17:09:38  abbott
// Added new fns  indets, CoefficientsWRT, ContentWRT.
//
// Revision 1.55  2011/11/09 13:50:01  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.54  2011/08/24 10:21:09  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.53  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.52  2011/07/05 15:02:17  bigatti
// -- added AlexanderDual
// -- added ad-hoc functions for colon, elim on monomial ideals
//
// Revision 1.51  2011/06/27 13:30:38  bigatti
// -- new file for monomial ideals:
// -- some functions moved there, but some have to be declared in SparsePolyRing.H
//
// Revision 1.50  2011/06/23 16:04:47  abbott
// Added IamExact mem fn for rings.
// Added myRecvTwinFloat mem fn for rings.
// Added first imple of RingHom from RingTwinFloat to other rings.
//
// Revision 1.49  2011/05/20 16:04:31  bigatti
// -- added QuotientBasis
//
// Revision 1.48  2011/05/19 13:56:12  abbott
// Removed comment out decl of old form of myIsRational.
//
// Revision 1.47  2011/04/19 13:59:49  bigatti
// -- added AreGensMonomial
//
// Revision 1.46  2011/04/12 09:52:15  bigatti
// -- added IsHomog(ideal), LT(ideal)
//
// Revision 1.45  2011/03/16 15:39:44  bigatti
// -- added myIsIndetPosPower(f), IsIndetPosPower(f)
//
// Revision 1.44  2011/03/16 13:21:23  abbott
// Added comments for myIsHomogPartial & myCmpWDegPartial.
// Cleaned up impls of myIsHomog & myIsHomogPartial.
// Corrected typo in impl of myElim.
//
// Revision 1.43  2011/03/11 10:53:38  bigatti
// -- added mySaturate
//
// Revision 1.42  2011/03/10 16:39:35  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.41  2011/03/01 14:10:47  bigatti
// -- added ClearDenom/myClearDenom
//
// Revision 1.40  2011/02/28 14:16:08  bigatti
// -- added GBasis(ideal)  -- only for SparsePolyRing
// -- error for myGcd when CoeffRing is not a field
//
// Revision 1.39  2011/01/28 17:58:07  bigatti
// -- added myElim
//
// Revision 1.38  2011/01/28 11:41:42  bigatti
// -- added IsPrintedWithMinus
// -- improved myOutput
// -- fixed bug in IsMinusOne
//
// Revision 1.37  2011/01/18 14:35:38  bigatti
// -- moved **_forC5 functions into CoCoA-5/CoCoALibSupplement:
//    myMonomials_forC5, mySupport_forC5, monomials_forC5, support_forC5,
//    LPP_forC5, LT_forC5, LM_forC5
//
// Revision 1.36  2010/11/30 11:34:37  bigatti
// -- added virtual myIndetSymbol
//
// Revision 1.35  2010/11/25 12:30:49  bigatti
// -- added myIndetsCalled
//
// Revision 1.34  2010/11/05 16:10:58  bigatti
// -- added myMonomials_forC5, mySupport_forC5
//
// Revision 1.33  2010/10/01 15:20:33  bigatti
// -- added mySymbolValue
// -- added RingElem(R, sym)
//
// Revision 1.32  2010/09/02 13:06:52  bigatti
// -- changed comment "fwd decl"
//
// Revision 1.31  2010/07/14 09:00:44  bigatti
// -- added comment
//
// Revision 1.30  2010/07/13 12:50:55  bigatti
// -- added comment
//
// Revision 1.29  2010/06/10 08:00:02  bigatti
// -- fixed naming conventions
//
// Revision 1.28  2010/03/18 13:55:56  abbott
// Added pseudo-ctors for monomials with BigRat coeffs.
//
// Revision 1.27  2010/03/05 18:43:48  abbott
// Added pseudo-ctors allowing polynomial rings to be created specifying
// the ordering using a PPOrderingCtor object.
//
// Revision 1.26  2010/02/04 09:57:11  bigatti
// -- added "mul" for ideals.  Implemented only for SparsePolyRing
//
// Revision 1.25  2010/01/21 14:12:28  bigatti
// -- swapped position of IhaveMonomialGens and myGBasisIsValid
//
// Revision 1.24  2009/11/26 17:21:34  bigatti
// -- added PushFront/PushBack(f, c, pp)
// -- added in .C inline functions:
// --   CheckCompatible, CheckElemSparsePolyRing, CheckCoeffExpv, CheckCoeffPP
//
// Revision 1.23  2009/10/02 13:27:26  bigatti
// -- unique implementation of myDiv in PolyRing.C
//
// Revision 1.22  2009/09/28 16:19:43  bigatti
// -- unique implementation for myDeriv
//
// Revision 1.21  2009/09/25 13:02:09  bigatti
// -- myDiv with one implementation in SparsePolyRing
//
// Revision 1.20  2009/09/22 14:07:33  bigatti
// -- added CmpWDegPartial and IsHomogPartial
//
// Revision 1.19  2009/07/24 12:26:43  abbott
// Added CommonDenom function for polynomials.
//
// Revision 1.18  2009/07/02 16:32:11  abbott
// Consequential changes stemming from new class BigRat, and modified interface to the member
// function RingBase::myIsRational.  Also some new conversion functions.
//
// Revision 1.17  2008/12/17 12:11:52  abbott
// Changed type from long to MachineInt in operations which use a machine integer
// in place of a RingElem.  The change is "superficial" but affects many files.
//
// Revision 1.16  2008/11/18 15:20:10  bigatti
// -- added const to myGBasis return value
// -- added myIdealCtor to RingWeyl for proper inheritance
//
// Revision 1.15  2008/10/09 15:15:05  bigatti
// -- added "bool3 IhaveMonomialGensValue" member field for ideals
//
// Revision 1.14  2008/04/21 12:32:54  abbott
// Corrected size_t into std::size_t in several header files; in some cases,
// replaced size_t with MachineInt (with consequent changes to impl files).
//
// Revision 1.13  2008/04/15 14:54:28  bigatti
// -- added mySquare
//
// Revision 1.12  2008/04/10 15:00:24  bigatti
// -- added  myPushBack/Front(RawPtr, ConstRawPtr, PPMonoidElemConstRawPtr)
//
// Revision 1.11  2008/03/12 16:35:18  bigatti
// -- changed: IsHomogeneous --> IsHomog
// -- changed: ERR:ZeroPoly --> ERR::ZeroRingElem
//
// Revision 1.10  2007/12/07 15:27:01  bigatti
// -- default implementation of "IamOne" in ideal.C
//
// Revision 1.9  2007/12/05 11:06:24  bigatti
// -- changed "size_t StdDeg/myStdDeg(f)" into "long"  (and related functions)
// -- changed "log/myLog(f, i)" into "MaxExponent/myMaxExponent(f, i)"
// -- fixed bug in "IsOne(ideal)" in SparsePolyRing.C
//
// Revision 1.8  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.7  2007/10/05 15:17:13  bigatti
// -- just a proper sorting on functions
//
// Revision 1.6  2007/05/31 16:31:53  bigatti
// -- adjusted comments
//
// Revision 1.5  2007/05/31 15:43:56  bigatti
// -- added mySymbols and AreGoodIndetNames
//
// Revision 1.3  2007/05/22 22:45:14  abbott
// Changed fn name IsUnit to IsInvertible.
//
// Revision 1.2  2007/03/12 16:00:29  bigatti
// -- moved myLog(F, index) into unique implementation in SparsePolyRing
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.23  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.22  2007/03/08 11:07:12  cocoa
// Made pseudo ctors for polynomial rings more uniform.  This allowed me to
// remove an include of CoCoA/symbol.H  from the RingDistrM*.H files, but then
// I had to put the include in several .C files.
//
// Revision 1.21  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.20  2007/01/20 14:07:25  bigatti
// -- moved code for homomorphism into common implementation in SparsePolyRing
//
// Revision 1.19  2007/01/15 13:33:34  cocoa
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.18  2006/12/07 17:36:19  cocoa
// -- migrated  myRemoveBigContent myContent myPowerSmallExp  into
//    single implementation in SparsePolyRing
// -- removed  content  from DistrMPoly(..)
//
// Revision 1.17  2006/12/06 17:25:43  cocoa
// -- removed #include "config.H"
//
// Revision 1.16  2006/11/23 17:43:53  cocoa
// -- minor change
//
// Revision 1.15  2006/11/22 17:51:31  cocoa
// -- moved printing functions into unified implementation in SparsePolyRing
//
// Revision 1.14  2006/11/21 18:09:24  cocoa
// -- added myIsMonomial
// -- implemented myIsOne, myIsMinusOne, myIsConstant, myIsIndet in SparsePolyRing
// -- removed the 4 functions from DistrMPoly(..) and RingDistrMPoly(..)
// -- changed all names of RawPtr arguments into "raw(..)"
//
// Revision 1.13  2006/11/20 15:55:02  cocoa
// ring is now a class again.  Improved definitions of operator-> in derived classes.
//
// Revision 1.12  2006/11/14 17:36:49  cocoa
// -- fixed implementation for ideal in RingWeyl
//
// Revision 1.11  2006/11/09 17:46:58  cocoa
// -- version 0.9712:
// --   IdealImpl moved to SparsePolyRing from concrete rings
// -- PolyList in GTypes is now vector<RingElem> (was list)
// -- "my" coding convention applied to DistrMPoly
//
// Revision 1.10  2006/11/03 14:01:46  cocoa
// -- changed: reference counting in ring, PPMonoids and OrdvArith now
//    uses SmartPtrIRC
//
// Revision 1.9  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.8  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.7  2006/08/17 09:45:07  cocoa
// -- added: homogenization
//
// Revision 1.6  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
// Revision 1.5  2006/07/20 17:06:08  cocoa
// -- moved myStdDeg into SparsePolyRing
//
// Revision 1.4  2006/07/20 16:51:38  cocoa
// -- added common implementation of myStdDeg
//
// Revision 1.3  2006/07/17 11:05:53  cocoa
// -- added: myIsValid, myIsHomogeneous, IsHomogeneous
//
// Revision 1.2  2006/06/08 16:45:28  cocoa
// -- RingDistrMPoly*.H  have been "moved" into RingDistrMPoly*.C
// -- some coding conventions fixed in DistrMPoly*
// -- functions wdeg and CmpWDeg have a common implementation in SparsePolyRing
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.10  2006/05/29 16:22:37  cocoa
// Third time lucky???
// Added myIsInteger member function to all rings (NYI for RingFloat).
//
// Revision 1.9  2006/05/12 17:01:28  cocoa
// -- added myIsUnit, myGcd whose implementation in RingDistr*** were identical
//
// Revision 1.8  2006/04/27 13:06:19  cocoa
// -- just spaces
//
// Revision 1.7  2006/04/26 16:44:53  cocoa
// -- myMul has now a single implementation in SparsePolyRing
// -- myMul and mul in RingDistrMPoly* and DistrMPoly* have been disabled
//
// Revision 1.6  2006/04/21 14:58:04  cocoa
// Removed myWDeg member function: it is no longer needed.
//
// Revision 1.5  2006/03/17 18:13:25  cocoa
// -- changed: myMul --> myMulByPP
//
// Revision 1.4  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.3  2006/03/07 10:11:06  cocoa
// -- fixed:  PPMonoidElem myLPP/LPP now return ConstRefPPMonoidElem
//
// Revision 1.2  2006/01/19 16:34:42  cocoa
// -- added NF, myReduceMod functions (not yet tested)
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.5  2005/09/30 15:03:39  cocoa
// Minor cleaning and tidying.
// DistrMPolyInlPP: use of summands now rather cleaner.
//
// Revision 1.4  2005/08/08 16:36:32  cocoa
// Just checking in before going on holiday.
// Don't really recall what changes have been made.
// Added IsIndet function for RingElem, PPMonoidElem,
// and a member function of OrdvArith.
// Improved the way failed assertions are handled.
//
// Revision 1.3  2005/07/19 15:30:20  cocoa
// A first attempt at iterators over sparse polynomials.
// Main additions are to SparsePolyRing, DistrMPoly*.
// Some consequential changes to PPMonoid*.
//
// Revision 1.2  2005/07/15 16:34:33  cocoa
// Added iterators for sparse polynomials.
// The code compiles (and the old tests still run).
// It'd Friday evening -- I'm going home before
// getting any ideas about making the iterator code run.
//
// Revision 1.1  2005/07/01 16:08:16  cocoa
// Friday check-in.  Major change to structure under PolyRing:
// now SparsePolyRing and DUPolyRing are separated (in preparation
// for implementing iterators).
//
// A number of other relatively minor changes had to be chased through
// (e.g. IndetPower).
//

#endif
