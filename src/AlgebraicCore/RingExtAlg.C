//   Copyright (c)  2018,2021  John Abbott, Anna Bigatti

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


// Source code for RingExtAlg (modified from RingWeylImpl)

#include "CoCoA/DistrMPolyClean.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/MatrixForOrdering.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PPOrdering.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingWeyl.H"
#include "CoCoA/TmpGReductor.H"
#include "CoCoA/TmpUniversalInvolutiveBasisContainer.H"
#include "CoCoA/assert.H"
#include "CoCoA/convert.H"
#include "CoCoA/ideal.H"
#include "CoCoA/matrix.H"
#include "CoCoA/symbol.H"
#include "CoCoA/utils.H"

#include <memory>
using std::unique_ptr;
#include <iostream>
using std::ostream;
using std::endl;   // just for debugging
//#include <vector>
using std::vector;



namespace CoCoA
{

  class RingExtAlgImpl: public SparsePolyRingBase  //???? should be NCRingBase  non-commutative!!!
  {

  public:
    //    RingWeylImpl(const ring& CoeffRing, long NumIndets, const std::vector<long>& ElimIndets);
    RingExtAlgImpl(const ring& CoeffRing, const std::vector<symbol>& names);
    ~RingExtAlgImpl();

  private: // data members of RingWeyl
    const SparsePolyRing myReprRing;
//    long myNumTrueIndetsValue;
    std::unique_ptr<RingElem> myZeroPtr;  ///< Every ring stores its own zero.
    std::unique_ptr<RingElem> myOnePtr;   ///< Every ring stores its own one.
    std::vector<RingElem> myIndetVector;
//    std::vector<RingElem> myDerivationVector;

  private:
    int myCoeffProdPP(PPMonoidElemConstRawPtr rawpp1, PPMonoidElemConstRawPtr rawpp2) const;

  public:
    //----------------------------------------------------------------------
    // Functions which every ring must implement:
    //----------------------------------------------------------------------
    bool IamCommutative() const override;
    bool3 IamIntegralDomain3(bool) const override;
    bool IamTrueGCDDomain() const override;
    ConstRefRingElem myZero() const override;
    ConstRefRingElem myOne() const override;
    RingElemRawPtr myNew() const override;
    RingElemRawPtr myNew(const MachineInt& n) const override;
    RingElemRawPtr myNew(const BigInt& N) const override;
    RingElemRawPtr myNew(ConstRawPtr rawcopy) const override;
    void myDelete(RawPtr rawx) const override;                             // destroys x (incl all resources)
    void mySwap(RawPtr rawx, RawPtr rawy) const override;                  // swap(x, y)
    void myAssign(RawPtr rawlhs, ConstRawPtr rawx) const override;         // lhs = x
    void myAssign(RawPtr rawlhs, const MachineInt& n) const override;      // lhs = n
    void myAssign(RawPtr rawlhs, const BigInt& N) const override;          // lhs = N
    void myAssignZero(RawPtr rawlhs) const override;                       // lhs = 0
    void myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const override;
    void myNegate(RawPtr rawlhs, ConstRawPtr rawx) const override;         // lhs = -x
    void myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x+y
    void mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x-y
    //    void myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x*y
    //    void myDiv(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;        // lhs = x/y
    bool myIsDivisible(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;// lhs = x/y, if divisible
    bool myIsInvertible(ConstRawPtr rawx) const override;                                // true iff x is invertible
    //    bool myIsPrintAtom(ConstRawPtr rawx) const override;
     void myGcd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const override;  // not a TrueGCDDomain
    void myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const override;// lhs = x^n, n>1, x not -1,0,1
    void mySymbols(std::vector<symbol>& SymList) const override;   // appends ring's symbols to SymList
    void myOutput(std::ostream& out, ConstRawPtr rawx) const override;      // out << x
    //    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;                // OMOut << R
    std::string myImplDetails() const override {return "RingWeyl";}
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const override;  // OMOut << x
    bool myIsZero(ConstRawPtr rawx) const override;                         // x == 0
    bool myIsOne(ConstRawPtr rawx) const override;                          // x == 1
    bool myIsMinusOne(ConstRawPtr rawx) const override;                     // x == -1
    //    bool myIsZeroAddMul: use default definition
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override;      // x == y

    ideal myIdealCtor(const std::vector<RingElem>& gens) const override;

    RingHom myCompose(const RingHom& phi, const RingHom& theta) const override; // phi(theta(...))


    // functions every PolyRing must implement
    long myNumIndets() const override;
    const ring& myCoeffRing() const override;
    const std::vector<RingElem>& myIndets() const override;
    void myIndetPower(RawPtr rawf, long var, long exp) const override;

    ///@name Simple functions on polynomials
    //@{
    long myNumTerms(ConstRawPtr rawf) const override;
    bool myIsConstant(ConstRawPtr rawf) const override;
    bool myIsIndet(long& index, ConstRawPtr rawf) const override;
    bool myIsMonomial(ConstRawPtr rawf) const override;
    long myStdDeg(ConstRawPtr rawf) const override;
    long myDeg(ConstRawPtr rawf, long var) const override;
    RingElemAlias myLC(ConstRawPtr rawf) const override;
    void myContent(RawPtr rawcontent, ConstRawPtr rawf) const override;
    void myRemoveBigContent(RawPtr rawf) const override;
    void myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< WEAK EXCEPTION GUARANTEE
    bool myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const override; ///< WEAK EXCEPTION GUARANTEE
    void myDeriv(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawx) const override; ///< lhs = deriv(f, x)
    //@}

//     friend const RingElem& indet(const RingWeyl& RW, long var);
//     friend const RingElem& derivation(const RingWeyl& RW, long var);

    RingHom myHomCtor(const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages) const override;

    //----------------------------------------------------------------------
    // Functions which every SparsePolyRing must implement:
    //----------------------------------------------------------------------

    const PPMonoid& myPPM() const override;

    ///@name   Functions for creating/building polynomials
    //@{
    RingElem myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const override;
    SparsePolyIter myBeginIter(ConstRawPtr rawf) const override;
    SparsePolyIter myEndIter(ConstRawPtr rawf) const override;

    void myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const override;
    void myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const override;
    void myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const override;
    void myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const override;
    //@}

    ///@name Special functions on polynomials needed for implementing Buchberger's Algorithm
    //@{
//     void myWDeg(degree& d, ConstRawPtr rawf) const override;
//     int myCmpWDeg(ConstRawPtr rawf, ConstRawPtr rawg) const override;
    ConstRefPPMonoidElem myLPP(ConstRawPtr rawf) const override;
    void myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const override;
    bool myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const override; ///< f+=LM(g); g-=LM(g); assumes LPP(f)==LPP(g); returns LC(f)+LC(g)==0
    void myMoveLMToFront(RawPtr rawf, RawPtr rawg) const override;
    void myMoveLMToBack(RawPtr rawf, RawPtr rawg) const override;
    void myDeleteLM(RawPtr rawf) const override; // ????? right interface
    void myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const override; ///< lhs=div(LM(f),LM(g)); assumes f!=0,g!=0
    int  myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const override; ///< cmp(LPP(f),LPP(g)); assumes f!=0,g!=0
    void myAddClear(RawPtr rawf, RawPtr rawg) const override; ///< f+=g; g=0;
    void myAppendClear(RawPtr rawf, RawPtr rawg) const override; ///< f+=g; g=0; appends g to f with no checks

    void myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const override; ///<  f += LM(h)*g
    void myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag) const override; ///<  f += LM(h)*g
    void myReductionStep(RawPtr rawf, ConstRawPtr rawg) const override;
    // ??? aggiungere coefficiente
    void myReductionStepGCD(RawPtr rawf, ConstRawPtr rawg, RingElem& FScale) const override;
    // should it all be in ReductionStep ??? ANNA
    //@}


    /////////////////////
    // >>> WARNING <<< //
    /////////////////////
    // DO NOT USE THE HOM CLASS -- I have not checked it
  private: // HomImpl class: overwrite some implementations of SparsePolyRing functions
    class HomImpl: public SparsePolyRingBase::HomImpl
    {
    public:
      HomImpl(const SparsePolyRing& domain, const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages);
    };

  private: // Ideal class: overwrite some implementations of SparsePolyRing functions
    class IdealImpl: public SparsePolyRingBase::IdealImpl
    {
    public:
      IdealImpl(const SparsePolyRing& P, const std::vector<RingElem>& gens);
      // default copy ctor is OK

      // functions every ideal must implement
      void myReduceMod(RingElemRawPtr rawr) const override; // r elem of R, where I is ideal in R
      bool IhaveElem(RingElemConstRawPtr rawr) const override;
      void myIntersect(const ideal&) override;
      void myColon(const ideal&) override;
      bool myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const override; // result is true iff quot exists & is unique; if true set lhs = num/den modulo I
      const std::vector<RingElem>& myGBasis(const CpuTimeLimit& CheckForTimeout) const override;

    protected:
      void myTestIsMaximal() const override;
      void myTestIsPrimary() const override;
      void myTestIsPrime() const override;
      void myTestIsRadical() const override;
    };
  };

  //----------------------------------------------------------------------

  // std::vector<symbol> WANAMES(long NumIndets)
  // {
  //   vector<symbol> ans = SymbolRange("x", 0, NumIndets-1);
  //   vector<symbol> d   = SymbolRange("d", 0, NumIndets-1);
  //   ans.insert(ans.end(), d.begin(), d.end());
  //   return ans;
  // }

  // std::vector<symbol> WANAMES(const std::vector<symbol>& names)
  // {
  //   vector<symbol> ans=names;
  //   const long NumNames = len(names);
  //   ans.reserve(2*NumNames);
  //   for ( long i=0 ; i < NumNames ; ++i )
  //   {
  //     if ( NumSubscripts(names[i])!=0 )
  //       CoCoA_THROW_ERROR("names must have no subscripts",
  //                   "WANAMES(const vector<symbol>& names)");
  //     ans.push_back(symbol("d"+head(names[i])));
  //   }
  //   return ans;
  // }


  RingExtAlgImpl::RingExtAlgImpl(const ring& CoeffRing, const std::vector<symbol>& VarNames):
    myReprRing(NewPolyRing(CoeffRing, VarNames))
//    myNumTrueIndetsValue(len(WANames)/2)
  {
//    const long NumTrueIndets=len(WANames)/2;
    CoCoA_ASSERT(IsCommutative(CoeffRing));
    myRefCountInc();  // this is needed for exception cleanliness, in case one of the lines below throws
    myZeroPtr.reset(new RingElem(ring(this)));
    myOnePtr.reset(new RingElem(ring(this), 1));
    myIndetVector.resize(NumIndets(myReprRing), *myZeroPtr);
    //    myDerivationVector.resize(myNumIndetsValue, *myZeroPtr);
    //    for (long i=0; i < myNumIndetsValue; ++i)
//    myIndetVector.resize(NumIndets(myReprRing), *myZeroPtr);
    vector<long> expv(NumIndets(myReprRing));

    for (long i=0; i < NumIndets(myReprRing); ++i)
    {
      expv[i] = 1;
      myPushFront(raw(myIndetVector[i]), raw(one(CoeffRing)), expv);
      expv[i] = 0;
    }
    myRefCountZero(); // otherwise it is 2 + NumIndets and won't be destroyed
  }


  RingExtAlgImpl::~RingExtAlgImpl()
  {}


//   inline const RingExtAlgImpl* RingWeyl::operator->() const
//   {
//     return static_cast<const RingExtAlgImpl*>(myRingPtr());
//   }


//   inline WeylPoly& ring_weyl::AsWeylPoly(RawPtr rawx) const
//   {
//     return *x.WeylPolyPtr;
//     //    return *static_cast<ring_weyl::WeylPoly*>(x.WeylPolyPtr);
//   }


//   inline const WeylPoly& ring_weyl::AsWeylPoly(ConstRawPtr rawx) const
//   {
//     return *x.WeylPolyPtr;
//     //    return *static_cast<ring_weyl::WeylPoly*>(x.WeylPolyPtr);
//   }


//   RingWeyl::RingWeyl(const RingExtAlgImpl* RingPtr):
//       SparsePolyRing(RingPtr)
//   {}


  //----------------------------------------------------------------------
  // Functions which every ring must implement:
  //----------------------------------------------------------------------

  bool RingExtAlgImpl::IamCommutative() const
  {
    return myNumIndets() == 0; // we are assuming the coeff ring is commutative
  }


  bool3 RingExtAlgImpl::IamIntegralDomain3(bool QuickMode) const
  {
    if (myNumIndets() > 0) return false3;
    return myReprRing->IamIntegralDomain3(QuickMode); //??? I think this is right
  }


  bool RingExtAlgImpl::IamTrueGCDDomain() const
  {
    return false; // I have no clue how to compute GCDs even if they do exist
  }


  ConstRefRingElem RingExtAlgImpl::myZero() const
  {
    return *myZeroPtr;
  }


  ConstRefRingElem RingExtAlgImpl::myOne() const
  {
    return *myOnePtr;
  }


  RingElemRawPtr RingExtAlgImpl::myNew() const
  {
    return myReprRing->myNew();
  }


  RingElemRawPtr RingExtAlgImpl::myNew(const MachineInt& n) const
  {
    return myReprRing->myNew(n);
  }


  RingElemRawPtr RingExtAlgImpl::myNew(const BigInt& N) const
  {
    return myReprRing->myNew(N);
  }


  RingElemRawPtr RingExtAlgImpl::myNew(ConstRawPtr rawcopy) const
  {
    return myReprRing->myNew(rawcopy);
  }


  void RingExtAlgImpl::myDelete(RawPtr rawx) const
  {
    myReprRing->myDelete(rawx);
  }


  void RingExtAlgImpl::mySwap(RawPtr rawx, RawPtr rawy) const
  {
    myReprRing->mySwap(rawx, rawy);
  }


  void RingExtAlgImpl::myAssign(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myReprRing->myAssign(rawlhs, rawx);
  }


  void RingExtAlgImpl::myAssign(RawPtr rawlhs, const MachineInt& n) const
  {
    myReprRing->myAssign(rawlhs, n);
  }


  void RingExtAlgImpl::myAssign(RawPtr rawlhs, const BigInt& N) const
  {
    myReprRing->myAssign(rawlhs, N);
  }


  void RingExtAlgImpl::myAssignZero(RawPtr rawlhs) const
  {
    myReprRing->myAssignZero(rawlhs);
  }


  void RingExtAlgImpl::myRecvTwinFloat(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    CoCoA_ASSERT(!IsExact(myReprRing));
    myReprRing->myRecvTwinFloat(rawlhs, rawx);
  }

  void RingExtAlgImpl::myNegate(RawPtr rawlhs, ConstRawPtr rawx) const
  {
    myReprRing->myNegate(rawlhs, rawx);
  }


  void RingExtAlgImpl::myAdd(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->myAdd(rawlhs, rawx, rawy);
  }


  void RingExtAlgImpl::mySub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    myReprRing->mySub(rawlhs, rawx, rawy);
  }


  //??? ANNA: I believe the general implementation in SparsePolyRing works for RingWeyl
//   void RingExtAlgImpl::myMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
//   {
//     if (myReprRing->myIsZero(x) || myReprRing->myIsZero(y)) { myReprRing->myAssignZero(rawlhs); return; }
//     RingElem f(myReprRing);
//     myReprRing->myAssign(raw(f), rawx);
//     ConstRefRingElem g(myReprRing, rawy);
// //    std::clog<<"f="<<f<<std::endl;
// //    std::clog<<"g="<<g<<std::endl;
//     RingElem ans(myReprRing); // compute answer in a temporary for exception safety
//     while (!CoCoA::IsZero(f))
//     {
//       RingElem LMf(myReprRing);
//       myMoveLMToFront(raw(LMf), raw(f));
//       PPMonoidElem LPPf = LPP(LMf);
//       RingElem LMfg(g);
//       myMulByPP(raw(LMfg), raw(LPPf));
// //      std::clog<<"LMf="<<LMf<<std::endl;
// //      std::clog<<"LMfg="<<LMfg<<std::endl;
//       myReprRing->myMulByCoeff(raw(LMfg), raw(LC(LMf))); //????UGLY!!!!
//       ans += LMfg;
//     }
//     mySwap(rawlhs, raw(ans));// really an assignment -- is this safe????
//   }


//   void RingExtAlgImpl::myDiv(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
//   {
//     CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::myDiv");
//     return;//???
// //     bool OK;                                            // Two lines o/w compiler complains that
// //     OK = CoCoA::WeylDiv(AsDMPI(rawlhs), AsDMPI(x), AsDMPI(y)); // OK is unused when debugging is off.
// //     CoCoA_ASSERT(OK);
//   }


  bool RingExtAlgImpl::myIsDivisible(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::myIsDivisible");
    return false;
//    return CoCoA::WeylDiv(AsDMPI(rawlhs), AsDMPI(x), AsDMPI(y));
  }


  bool RingExtAlgImpl::myIsInvertible(ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::myIsInvertible");  //??? d[1]*x[1] = 1
    return false;//???
  }


//   bool RingExtAlgImpl::myIsPrintAtom(ConstRawPtr rawx) const
//   {
//     return myReprRing->myIsPrintAtom(rawx);  //??? default always false
//   }


  void RingExtAlgImpl::myGcd(RawPtr /*rawlhs*/, ConstRawPtr /*rawx*/, ConstRawPtr /*rawy*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::myGcd");
  }


  void RingExtAlgImpl::myPowerSmallExp(RawPtr rawlhs, ConstRawPtr rawx, long n) const
  {
    // Assert that we have a genuinely non-trivial case.
    CoCoA_ASSERT(n > 1);
    CoCoA_ASSERT(!myIsZero(rawx) && !myIsOne(rawx) && !myIsMinusOne(rawx));

    mySequentialPower(rawlhs, rawx, n); // probably better than myBinaryPower, I guess.
  }


  void RingExtAlgImpl::mySymbols(std::vector<symbol>& SymList) const
{ //??? ANNA: which symbols are in a RingExtAlgImpl??
    myReprRing->mySymbols(SymList);
  }


  void RingExtAlgImpl::myOutput(std::ostream& out, ConstRawPtr rawx) const
  {
    if (!out) return;  // short-cut for bad ostreams
    myReprRing->myOutput(out, rawx);
  }


//   void RingExtAlgImpl::myOutputSelf(std::ostream& out) const
//   {
//     out << "RingExtAlgImpl(" << CoeffRing(myReprRing) << ", " << NumIndets(myReprRing) << ")";
//     //?????????  WHAT ABOUT THE ORDERING AND GRADING????
//   }


//   void RingExtAlgImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
//   {
//     OMOut->mySendApplyStart();
//     OMOut << OpenMathSymbol("polyd", "poly_ring_d"); //??? weyl?
//     OMOut << CoeffRing(myReprRing);
//     OMOut << NumIndets(myReprRing); //???? losing the ordering and grading info here!!!
//     OMOut->mySendApplyEnd();
//   }


  void RingExtAlgImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawx) const
  {
    myReprRing->myOutput_OM(OMOut, rawx);
  }


  bool RingExtAlgImpl::myIsZero(ConstRawPtr rawx) const
  {
    return myReprRing->myIsZero(rawx);
  }


  bool RingExtAlgImpl::myIsOne(ConstRawPtr rawx) const
  {
    return myReprRing->myIsOne(rawx);
  }


  bool RingExtAlgImpl::myIsMinusOne(ConstRawPtr rawx) const
  {
    return myReprRing->myIsMinusOne(rawx);
  }


  // bool RingExtAlgImpl::myIsZeroAddMul(RawPtr rawlhs, ConstRawPtr rawy, ConstRawPtr rawz) const; // use default definition


  bool RingExtAlgImpl::myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const
  {
    return myReprRing->myIsEqual(rawx, rawy);
  }


  RingHom RingExtAlgImpl::myCompose(const RingHom& /*phi*/, const RingHom& theta) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::compose");
    return theta; // just to keep compiler quiet
  }

  //----------------------------------------------------------------------

//   long NumIndets(const RingWeyl& RW)
//   {
//     return RW->myNumIndets();
//   }

  //----------------------------------------------------------------------
  // Functions which every PolyRing must implement
  //----------------------------------------------------------------------

  long RingExtAlgImpl::myNumIndets() const
  {
    return myReprRing->myNumIndets();
  }


  const ring& RingExtAlgImpl::myCoeffRing() const
  {
    return myReprRing->myCoeffRing();
  }


  const std::vector<RingElem>& RingExtAlgImpl::myIndets() const
  {
    return myIndetVector;
  }


  void RingExtAlgImpl::myIndetPower(RawPtr rawf, long var, long exp) const
  {
    myReprRing->myIndetPower(rawf, var, exp);
  }


  long RingExtAlgImpl::myNumTerms(ConstRawPtr rawx) const
  {
    return myReprRing->myNumTerms(rawx);
  }


  bool RingExtAlgImpl::myIsConstant(ConstRawPtr rawf) const
  {
    return myReprRing->myIsConstant(rawf);
  }


  bool RingExtAlgImpl::myIsIndet(long& index, ConstRawPtr rawf) const
  {
    return myReprRing->myIsIndet(index, rawf);
  }


  bool RingExtAlgImpl::myIsMonomial(ConstRawPtr rawf) const
  {
    return myReprRing->myIsMonomial(rawf);
  }


  long RingExtAlgImpl::myStdDeg(ConstRawPtr rawf) const
  {
    return myReprRing->myStdDeg(rawf);
  }


  long RingExtAlgImpl::myDeg(ConstRawPtr rawf, long var) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return myReprRing->myDeg(rawf, var);//????????
  }


  RingElemAlias RingExtAlgImpl::myLC(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return myReprRing->myLC(rawf);//???????????
  }


  void RingExtAlgImpl::myContent(RawPtr rawdest, ConstRawPtr rawf) const
  {
    myReprRing->myContent(rawdest, rawf);
  }


  void RingExtAlgImpl::myRemoveBigContent(RawPtr rawf) const
  {
    myReprRing->myRemoveBigContent(rawf);
  }


  void RingExtAlgImpl::myMulByCoeff(RawPtr rawf, ConstRawPtr rawc) const // WEAK EXCEPTION GUARANTEE
  {
    myReprRing->myMulByCoeff(rawf, rawc);
  }


  bool RingExtAlgImpl::myDivByCoeff(RawPtr rawf, ConstRawPtr rawc) const // WEAK EXCEPTION GUARANTEE
  {
    return myReprRing->myDivByCoeff(rawf, rawc);
  }


  void RingExtAlgImpl::myDeriv(RawPtr /*rawlhs*/, ConstRawPtr /*rawf*/, ConstRawPtr /*rawx*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::myDeriv");
  }


  RingHom RingExtAlgImpl::myHomCtor(const ring& /*codomain*/, const RingHom& /*CoeffHom*/, const std::vector<RingElem>& /*IndetImages*/) const
  {
    CoCoA_THROW_ERROR("DOES NOT EXIST", "RingExtAlgImpl::myHomCtor");
    return IdentityHom(myReprRing); // just to keep compiler quiet
  }



  //----------------------------------------------------------------------
  // Functions which every SparsePolyRing must implement:
  //----------------------------------------------------------------------

  const PPMonoid& RingExtAlgImpl::myPPM() const
  {
    return myReprRing->myPPM();
  }


  RingElem RingExtAlgImpl::myMonomial(ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    RingElem ans(ring(this));
    RingElem m = myReprRing->myMonomial(rawc, rawpp);
    mySwap(raw(ans), raw(m));
    return ans;
  }


  SparsePolyIter RingExtAlgImpl::myBeginIter(ConstRawPtr rawf) const
  {
    return myReprRing->myBeginIter(rawf);
  }


  SparsePolyIter RingExtAlgImpl::myEndIter(ConstRawPtr rawf) const
  {
    return myReprRing->myEndIter(rawf);
  }


  void RingExtAlgImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    myReprRing->myPushFront(rawf, rawc, expv);//???  how many elements does expv have???
  }


  void RingExtAlgImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
  {
    myReprRing->myPushBack(rawf, rawc, expv);//???  how many elements does expv have???
  }


  void RingExtAlgImpl::myPushFront(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    myReprRing->myPushFront(rawf, rawc, rawpp);
  }


  void RingExtAlgImpl::myPushBack(RawPtr rawf, ConstRawPtr rawc, PPMonoidElemConstRawPtr rawpp) const
  {
    myReprRing->myPushBack(rawf, rawc, rawpp);
  }


  ConstRefPPMonoidElem RingExtAlgImpl::myLPP(ConstRawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    return myReprRing->myLPP(rawf);
  }


    int RingExtAlgImpl::myCoeffProdPP(PPMonoidElemConstRawPtr rawpp1, PPMonoidElemConstRawPtr rawpp2) const
      {
        const PPMonoid& PParith = PPM(myReprRing);
        if (!PParith->myIsCoprime(rawpp1, rawpp2)) return 0;
        const int n = NumIndets(PParith);
        vector<long> expv(n);
        PParith->myExponents(expv, rawpp1);
        vector<long> Rcount(n);
        for (int i=0; i < n; ++i)
        {
          if (expv[i] == 0) continue;
          for (int j=0; j <= i; ++j)
            ++Rcount[j];
        }
        PParith->myExponents(expv, rawpp2);
        int count = 0;
        for (int i=0; i < n; ++i)
          if (expv[i] > 0)
            count += Rcount[i];
        if ((count&1) == 1) return -1;
        return 1;
      }

  // f = pp*f  NOTE: pp is on the LEFT and f is on the RIGHT
  void RingExtAlgImpl::myMulByPP(RawPtr rawf, PPMonoidElemConstRawPtr rawpp) const
  {
    if (myReprRing->myIsZero(rawf)) return;

    RingElem g(myReprRing);
//    myReprRing->myAssign(raw(g), rawf);
//???    std::clog<<"Alias(f)="<<g<<std::endl;
//    const long nvars = myNumIndets();
    const PPMonoid& PParith = PPM(myReprRing);
    PPMonoidElem ProdPP(PParith);
    for (SparsePolyIter it = myReprRing->myBeginIter(rawf); !IsEnded(it); ++it)
    {
      const int PPcoeff = myCoeffProdPP(rawpp, raw(PP(it)));
      if (PPcoeff == 0) continue;
      const RingElem c = PPcoeff*coeff(it);
      PParith->myMul(raw(ProdPP), rawpp, raw(PP(it)));
      PushBack(g, c, ProdPP);
    }

    mySwap(rawf, raw(g));// really an assignment -- is this safe????
  }


  bool RingExtAlgImpl::myIsZeroAddLCs(RawPtr rawf, RawPtr rawg) const
  {
    return myReprRing->myIsZeroAddLCs(rawf, rawg);
  }


  void RingExtAlgImpl::myMoveLMToFront(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    myReprRing->myMoveLMToFront(rawf, rawg);
  }


  void RingExtAlgImpl::myMoveLMToBack(RawPtr rawf, RawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawg));
    myReprRing->myMoveLMToBack(rawf, rawg);
  }


  void RingExtAlgImpl::myDeleteLM(RawPtr rawf) const
  {
    CoCoA_ASSERT(!myIsZero(rawf));
    myReprRing->myDeleteLM(rawf);
  }


  void RingExtAlgImpl::myDivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    myReprRing->myDivLM(rawlhs, rawf, rawg);
  }


  int  RingExtAlgImpl::myCmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const
  {
    CoCoA_ASSERT(!myIsZero(rawf) && !myIsZero(rawg));
    return myReprRing->myCmpLPP(rawf, rawg);
  }


  void RingExtAlgImpl::myAddClear(RawPtr rawf, RawPtr rawg) const
  {
    myReprRing->myAddClear(rawf, rawg);
  }


  void RingExtAlgImpl::myAppendClear(RawPtr rawf, RawPtr rawg) const
  {
    myReprRing->myAppendClear(rawf, rawg);
  }


  void RingExtAlgImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg) const //???? delete me???
  {
    myAddMulLM(rawf, rawh, rawg, DontSkipLMg);
  }


  void RingExtAlgImpl::myAddMulLM(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, SkipLMFlag skip) const
  {
    CoCoA_ASSERT(myNumTerms(rawh)==1);
    RingElem prod(ring(this), myNew(rawg));
    myMulByCoeff(raw(prod), raw(myLC(rawh)));
    myMulByPP(raw(prod), raw(myLPP(rawh)));
    myReprRing->myAddMulLM(rawf, raw(myOne()), raw(prod), skip);
  }


  void RingExtAlgImpl::myReductionStep(RawPtr rawf, ConstRawPtr rawg) const
  {
    PPMonoidElem q = myLPP(rawf)/myLPP(rawg);
    RingElem c = -myLC(rawf)/myLC(rawg); //    RingElem c = -LC(repf)/LC(repg);
    RingElem qg(myReprRing); myReprRing->myAssign(raw(qg), rawg); // qg is copy of rawg
    myMulByPP(raw(qg), raw(q));  // qg = q * g  //??? should we use myAddMul
    myMulByCoeff(raw(qg), raw(c));
    myAdd(rawf, rawf, raw(qg)); //  repf += qg;


/*
//???    ConstRefRingElem aliasg(ring(this), rawg);
//???    std::clog<<"g="<<aliasg<<std::endl;
    PPMonoidElem LPPf = myReprRing->myLPP(f);
    PPMonoidElem LPPg = myReprRing->myLPP(g);
    PPMonoidElem q = LPPf/LPPg; // quotient exists
    RingElem LCf(CoeffRing(myReprRing));
    CoeffRing(myReprRing)->myAssign(raw(LCf), myReprRing->myRawLC(f));
    RingElem LCg(CoeffRing(myReprRing));
    CoeffRing(myReprRing)->myAssign(raw(LCg), myReprRing->myRawLC(g));
    RingElem c = -LCf/LCg;
///    RingElem c(myCoeffRing());
///    myCoeffRing()->myDiv(raw(c), myReprRing->myRawLC(f), myReprRing->myRawLC(g));
    RingElem g1(myReprRing);
//???    std::clog<<"ReductionStep: q="<<q<<std::endl;
//???    std::clog<<"ReductionStep: c="<<c<<std::endl;
    myReprRing->myAssign(raw(g1), rawg);
//???    std::clog<<"ReductionStep: g="<<g1<<std::endl;
    myMul(raw(g1), raw(q));  // g1 = q*g;
    myReprRing->myMulByCoeff(raw(g1), raw(c));
//???    std::clog<<"ReductionStep: prod="<<g1<<std::endl;
    {
//???      ConstRefRingElem aliasf(ring(this),f);
//???      std::clog<<"REDUCING f="<<aliasf<<std::endl;
//???      std::clog<<"REDN by  g="<<aliasg<<std::endl;
      myAdd(f, f, raw(g1));
//???    std::clog<<"RESULT is  "<<aliasf<<std::endl<<std::endl;
    }
*/
  }


  void RingExtAlgImpl::myReductionStepGCD(RawPtr /*rawf*/, ConstRawPtr /*rawg*/, RingElem& /*fscale*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::ReductionStepGCD");
  }


  //---------------------------------------------------------------------------
  // Functions to do with RingExtAlgImpl::HomImpl


  RingExtAlgImpl::HomImpl::HomImpl(const SparsePolyRing& domain, const ring& codomain, const RingHom& CoeffHom, const std::vector<RingElem>& IndetImages):
    SparsePolyRingBase::HomImpl(domain, codomain, CoeffHom, IndetImages)
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::HomImpl::HomImpl");
  }



  //---------------------------------------------------------------------------
  // Functions for the class RingExtAlgImpl::IdealImpl

  // inheritance is delicate here: this function is **necessary**
  ideal RingExtAlgImpl::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new IdealImpl(SparsePolyRing(this), gens)); //??? ugly ???
  }


  RingExtAlgImpl::IdealImpl::IdealImpl(const SparsePolyRing& P, const std::vector<RingElem>& gens):
    SparsePolyRingBase::IdealImpl(P, gens)
  {}


  void RingExtAlgImpl::IdealImpl::myTestIsMaximal() const
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::IdealImpl::myTestIsMaximal()"); }


  void RingExtAlgImpl::IdealImpl::myTestIsPrimary() const
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::IdealImpl::myTestIsPrimary()"); }


  void RingExtAlgImpl::IdealImpl::myTestIsPrime() const
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::IdealImpl::myTestIsPrime()"); }


  void RingExtAlgImpl::IdealImpl::myTestIsRadical() const
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::IdealImpl::myTestIsRadical()"); }


  void RingExtAlgImpl::IdealImpl::myReduceMod(RingElemRawPtr /*rawr*/) const
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::WeylIdeal::IdealImpl::myReduceMod"); }


  bool RingExtAlgImpl::IdealImpl::IhaveElem(RingElemConstRawPtr /*rawr*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::WeylIdeal::IhaveElem");
    return false; /* just to keep compiler quiet */
  }


  void RingExtAlgImpl::IdealImpl::myIntersect(const ideal& /*J*/)
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::WeylIdeal::myIntersect"); }//???


  void RingExtAlgImpl::IdealImpl::myColon(const ideal& /*J*/)
  { CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::WeylIdeal::myColon"); }//???


  bool RingExtAlgImpl::IdealImpl::myDivMod(RingElemRawPtr /*rawlhs*/, RingElemConstRawPtr /*rawnum*/, RingElemConstRawPtr /*rawden*/) const
  {
    CoCoA_THROW_ERROR(ERR::NYI, "RingExtAlgImpl::WeylIdeal::myDivMod");
    return false; /* just to keep compiler quiet */
  }


  const std::vector<RingElem>& RingExtAlgImpl::IdealImpl::myGBasis(const CpuTimeLimit& CheckForTimeout) const
  {
    if (IhaveGBasis()) return myGBasisValue;
    CoCoA_ASSERT(myGBasisValue.empty());
    vector<RingElem> GensList, GBList;
    GensList.insert(GensList.end(), myGensValue.begin(), myGensValue.end());
    bool IsHomog=false;
    bool IsSatAlg=false;
    GRingInfo GRI(myP, IsHomog, IsSatAlg, NewDivMaskEvenPowers(), CheckForTimeout);
    GBCriteria criteria(GBCriteria::DontUseCoprime, GBCriteria::UseGM,
                        GBCriteria::UseBack, GBCriteria::DontUseDiv);
    GReductor GR(GRI, GensList,
                 GReductor::AffineAlg,
                 //                 Reductors::DontUseBorel, GReductor::DontUseDynamicAlg,
                 criteria);
    GR.myDoGBasis();
    //    GR.myDoAFFGBasis();
    GR.myCopyGBasis(GBList);
    myGBasisValue.insert(myGBasisValue.end(), GBList.begin(), GBList.end());
    IhaveGBasisValue = true;
    return myGBasisValue;
  }


  //----------------------------------------------------------------------
  // Pseudo-ctors for (sparse) polynomial rings.


  SparsePolyRing NewExtAlgebra(const ring& CoeffRing, long NumIndets)
  {
    return SparsePolyRing(new RingExtAlgImpl(CoeffRing, SymbolRange("x",1,NumIndets)));
  }


  SparsePolyRing NewExtAlgebra(const ring& CoeffRing, const std::vector<symbol>& names)
  {
    return SparsePolyRing(new RingExtAlgImpl(CoeffRing, names));
  }


//   const RingElem& indet(const RingWeyl& RW, long var)
//   {
//     if (var > CoCoA::NumIndets(RW)) CoCoA_THROW_ERROR("Indeterminate index too large", "indet(RW,var)");
//     return RW->myIndetVector[var];
//   }


//   const RingElem& derivation(const RingWeyl& RW, long var)
//   {
//     if (var > CoCoA::NumIndets(RW)) CoCoA_THROW_ERROR("Indeterminate index too large", "derivation(RW,var)");
//     return RW->myDerivationVector[var];
//   }


} // end of namespace CoCoA

// #error "JUNK AFTER THIS LINE"
// ?????????????????????????????????????????????????????????????????????????????
// #include <algorithm>
// #include <vector>
// #include <iostream>

// using std::max;
// using std::swap;
// using std::vector;
// using std::ostream;
// using std::endl;   // just for debugging


// #include "CoCoA/deriv.H"
// #include "CoCoA/ring_weyl.H"

// namespace  // unnamed, for file local functions
// {

//   using namespace CoCoA;


// //----------------------------------------------------------------------


// //  // f = pp*g where the product is in the Weyl algebra
// //  void act(RingElem& f, const RingElem& pp, const RingElem& g)
// //  {
// //    ASSERT(IsPolyRing(owner(f)));
// //    const PolyRing& P = PolyRing::owner(f);
// //    ASSERT(&owner(g) == &P && &owner(pp) == &P);
// //    ASSERT(!IsZero(pp) && NumTerms(pp) == 1);
// //    //  std::clog<<"entered ucha's action" << endl;
// //    //  std::clog << "pp=" << pp << endl;
// //    f = g;
// //    if (IsZero(f)) return;
// //    const long nvars = NumIndets(P);

// //    vector<long> expv(nvars);
// //    PPM(P).exponents(expv, raw(LPP(pp)));

// //    for (long var = nvars/2; var < nvars; ++var)
// //    {
// //      //    std::clog << "doing D var=" << var << endl;
// //      const long n = expv[var];
// //      //    std::clog << "order = " << n << endl;
// //      RingElem derf = f;
// //      f *= P.IndetPower(var, n);
// //      for (long i=1; i <= n; ++i)
// //      {
// //        derf = deriv(derf, var-nvars/2);
// //        f += binomial(n, i)*derf*P.IndetPower(var, n-i); // *IndetPower(h, 2*i); // for homog case
// //        std::clog<<"binomial("<<n<<","<<i<<")="<<binomial(n,i)<<std::endl;
// //      }
// //    }
// //    { // f *= var^deg(pp, var); for the early vars
// //      for (long var = nvars/2; var < nvars; ++var)
// //        expv[var] = 0;
// //      PPMonoid::elem qq(PPM(P), expv);
// //      f *= P.monomial(RingElem(CoeffRing(P), 1), qq);
// //    }
// //  }
// //----------------------------------------------------------------------
// } // end of unnamed namespace


// namespace CoCoA
// {




//   RingExtAlgImpl::ring_weyl(const AbstractRing& R, const PPMonoid& PPM):
//     myCoeffRing(R),
//     myPPM(PPM),
//     myWeylPolyPool(sizeof(WeylPoly), "RingExtAlgImpl::myWeylPolyPool"),
//     myNumIndets(NumIndets(PPM)),
//     myOrdvWords(OrdvWords(ordering(PPM))),
//     mySummandSize(sizeof(WeylPoly::summand) + sizeof(int)*(myOrdvWords-1)),
//     mySummandPool(mySummandSize, "RingExtAlgImpl::mySummandPool")
//   {
//     myNumPolys = 0;
//     myZero = new RingElem(*this);
//     myIndetVector.resize(myNumIndets, RingElem(*this));
//     vector<long> expv(myNumIndets);
//     RingElem one(myCoeffRing, 1);
//     for (long i=0; i < myNumIndets; ++i)
//     {
//       expv[i] = 1;
//       PushFront(raw(myIndetVector[i]), raw(one), expv);
//       expv[i] = 0;
//     }
//   }


//   RingExtAlgImpl::~ring_weyl()
//   {
//     myIndetVector.clear();
//     delete myZero;
//     ASSERT(myNumPolys == 0);
//   }


//   void RingExtAlgImpl::MakeWritable(RawPtr rawx) const
//   {
//     if (x.WeylPolyPtr->myRefCount == 1) return;
//     --x.WeylPolyPtr->myRefCount;
//     WeylPoly* copy = static_cast<WeylPoly*>(myWeylPolyPool.alloc(sizeof(WeylPoly)));
//     const WeylPoly* orig = x.WeylPolyPtr;
//     x.WeylPolyPtr = copy;
//     copy->myRefCount = 1;
//     copy->mySummands = nullptr;
//     copy->myEnd = &copy->mySummands;
//     for (const WeylPoly::summand* it = orig->mySummands; it != nullptr; it = it->myNext)
//     {
//       // more or less hand inlined body of PushBack -- NB using PushBack is "circular"
//       *copy->myEnd = CopySummand(it);
//       copy->myEnd = &((*copy->myEnd)->myNext);
//     }

//     ///  myCoeffRing.init(copy->myDenom);
// //    new (&copy->myLC) RingElem(myCoeffRing, AbstractRing::elem::ALIAS);
//   }


//   inline WeylPoly::summand* RingExtAlgImpl::AllocSummand() const
//   {
//     return static_cast<WeylPoly::summand*>(mySummandPool.alloc(mySummandSize));
//   }


//   WeylPoly::summand* RingExtAlgImpl::InitSummand(WeylPoly::summand* ptr) const
//   {
//     myCoeffRing.init(ptr->myCoeff);
//     ptr->myNext = nullptr;
//     ptr->myModulePosn = 0;
//     return ptr;
//   }


//   WeylPoly::summand* RingExtAlgImpl::InitSummand(WeylPoly::summand* ptr, ConstRawPtr rawc, const vector<long>& expv) const
//   {
//     myCoeffRing.init(ptr->myCoeff, c);
//     ptr->myNext = nullptr;
//     ptr->myModulePosn = 0;
//     ordering(myPPM).ComputeOrdv(ptr->myOrdv, &expv[0]); // &expv[0] converts vector<T> to T*
//     return ptr;
//   }


//   WeylPoly::summand* RingExtAlgImpl::CopySummand(const WeylPoly::summand* original) const
//   {
//     WeylPoly::summand* copy = AllocSummand();
//     copy->myNext = nullptr;
//     myCoeffRing.init(copy->myCoeff, original->myCoeff);
//     copy->myModulePosn = original->myModulePosn;
//     for (long i=0; i < myOrdvWords; ++i)
//       copy->myOrdv[i] = original->myOrdv[i];
//     return copy;
//   }


//   void RingExtAlgImpl::SetSummandMOrdv(WeylPoly::summand* dest, const WeylPoly::summand* src) const
//   {
//     dest->myModulePosn = src->myModulePosn;
//     for (long i=0; i < myOrdvWords; ++i)
//       dest->myOrdv[i] = src->myOrdv[i];
//   }


//   void RingExtAlgImpl::DeleteSummands(WeylPoly::summand* ptr) const
//   {
//     WeylPoly::summand* next;
//     while (ptr != 0)
//     {
//       next = ptr->myNext;
//       myCoeffRing.kill(ptr->myCoeff);
//       mySummandPool.free(ptr, mySummandSize);
//       ptr = next;
//     }
//   }


//   bool RingExtAlgImpl::EqualSummands(const WeylPoly::summand& lhs, const WeylPoly::summand& x) const
//   {
//     if (lhs.myModulePosn != x.myModulePosn) return false;
//     const PPOrdering::OrdvElem* const lordv = lhs.myOrdv;
//     const PPOrdering::OrdvElem* const rordv = x.myOrdv;
//     for (long i = 0; i < myOrdvWords; ++i)
//       if (lordv[i] != rordv[i]) return false;
//     return myCoeffRing.IsEqual(lhs.myCoeff, x.myCoeff);
//   }


//   inline void RingExtAlgImpl::MulOrdv(PPOrdering::OrdvElem* ov, const PPOrdering::OrdvElem* ov1, const PPOrdering::OrdvElem* ov2) const
//   {
//     for (long i=0; i < myOrdvWords; ++i)
//       ov[i] = ov1[i]+ov2[i];
//   }


//   inline void RingExtAlgImpl::DivOrdv(PPOrdering::OrdvElem* ov, const PPOrdering::OrdvElem* ov1, const PPOrdering::OrdvElem* ov2) const
//   {
//     for (long i=0; i < myOrdvWords; ++i)
//       ov[i] = ov1[i]-ov2[i];
//   }




//   long RingExtAlgImpl::NumIndets() const
//   {
//     return myNumIndets;
//   }


//   const AbstractRing& RingExtAlgImpl::CoeffRing() const
//   {
//     return myCoeffRing;
//   }


//   const PPMonoid& RingExtAlgImpl::PPM() const
//   {
//     return myPPM;
//   }


//   long RingExtAlgImpl::NumTerms(ConstRawPtr rawx) const
//   {
//     long nsummands = 0;
//     for (const WeylPoly::summand* it = AsWeylPoly(rawx).mySummands; it != nullptr; it = it->myNext) ++nsummands;
//     return nsummands;
//   }


//   PolyIter RingExtAlgImpl::BeginIter(AbstractRing::RawPtr rawx) const
//   {
//     return PolyIter(*this, &AsWeylPoly(rawx).mySummands);
//   }


//   PolyIter RingExtAlgImpl::EndIter(AbstractRing::RawPtr rawx) const
//   {
//     return PolyIter(*this, AsWeylPoly(rawx).myEnd);
//   }


//   const RingElem& RingExtAlgImpl::indet(long var) const
//   {
//     ASSERT(var < myNumIndets);
//     return myIndetVector[var];
//   }


//   RingElem RingExtAlgImpl::IndetPower(long var, long n) const
//   {
//     ASSERT(0 <= var && var < myNumIndets);
//     return monomial(CoCoA::IndetPower(myPPM, var, n));
//   }


//   RingElem RingExtAlgImpl::monomial(const RingElem& c, const PPMonoid::alias& pp) const
//   {
//     vector<long> expv(myNumIndets);
//     myPPM.exponents(expv, raw(pp));
//     RingElem ans(*this);
//     PushFront(raw(ans), raw(c), expv);
//     return ans;
//   }


//   RingElem RingExtAlgImpl::monomial(const PPMonoid::alias& pp) const
//   {
//     RingElem c(myCoeffRing, 1);
//     vector<long> expv(myNumIndets);
//     myPPM.exponents(expv, raw(pp));
//     RingElem ans(*this);
//     PushFront(raw(ans), raw(c), expv);
//     return ans;
//   }


//   RingElem RingExtAlgImpl::monomial(const RingElem& c) const
//   {
//     vector<long> expv(myNumIndets);
//     RingElem ans(*this);
//     PushFront(raw(ans), raw(c), expv);
//     return ans;
//   }


//   long RingExtAlgImpl::deg(ConstRawPtr rawf) const
//   {
//     if (IsZero(f)) CoCoA_THROW_ERROR("RingExtAlgImpl::deg: cannot compute degree of zero polynomial");
//     if (GradingDim(myPPM) > 0)
//       return ordering(myPPM).deg(AsWeylPoly(f).mySummands->myOrdv); //// BUG????  "valid" only if grading dim == 1
//     //    return ???; the vector in Z^0 -- ring not graded!!!
//     return -1;//??????????  BUG BUG INCOMPLETE
//   }


//   int RingExtAlgImpl::deg(ConstRawPtr rawf, long var) const
//   {
//     if (IsZero(f)) CoCoA_THROW_ERROR("RingExtAlgImpl::deg: cannot compute degree of zero polynomial");
//     long d = 0;
//     for (WeylPoly::summand* it = AsWeylPoly(f).mySummands; it != nullptr; it = it->myNext)
//       d = max(d, ordering(myPPM).exponent(it->myOrdv, var));
//     return d;
//   }


// //    int RingExtAlgImpl::deg(ConstRawPtr rawf, const PPGrading& G) const
// //    {
// //      if (IsZero(f)) CoCoA_THROW_ERROR("RingExtAlgImpl::deg: cannot compute degree of zero polynomial");
// //      int d = 0;
// //      for (WeylPoly::summand* it = AsWeylPoly(f).mySummands; it != nullptr; it = it->myNext)
// //        d = max(d, deg(it, G)); // CANNOT JUST USE MAX HERE!!  MUST USE LEX MAX FOR VECTORS
// //      return d;
// //    }



//   const RingElem RingExtAlgImpl::LC(ConstRawPtr rawf) const
//   {
//     ASSERT(!IsZero(f));
// //    if (IsZero(f)) return ::zero(myCoeffRing);
//     return RingElem(myCoeffRing, AsWeylPoly(f).mySummands->myCoeff);
// //    AsWeylPoly(f).myLC.reseat(AsWeylPoly(f).mySummands->myCoeff);
// //    return AsWeylPoly(f).myLC;
//   }


//   AbstractRing::RawPtr RingExtAlgImpl::RawLC(RawPtr rawf) const
//   {
//     MakeWritable(f);//?????????
//     ASSERT(!IsZero(f));
//     return AsWeylPoly(f).mySummands->myCoeff;
//   }


//   const AbstractRing::RawPtr RingExtAlgImpl::RawLC(ConstRawPtr rawf) const
//   {
//     ASSERT(!IsZero(f));
//     return AsWeylPoly(f).mySummands->myCoeff;
//   }


//   RingElem RingExtAlgImpl::content(ConstRawPtr rawf) const
//   {
//     ASSERT(::IsTrueGCDDomain(myCoeffRing));
//     ASSERT(!IsZero(f));
//     RingElem cont(myCoeffRing);
//     for (WeylPoly::summand* it = AsWeylPoly(f).mySummands; it != nullptr; it = it->myNext)
//       myCoeffRing.gcd(raw(cont), raw(cont), it->myCoeff);  // be clever if cont == 1??
//     return cont;
//   }



//   void RingExtAlgImpl::MoveLMToFront(RawPtr rawf, RawPtr rawg) const
//   {
//     ASSERT(!IsZero(g));
//     ASSERT(f.WeylPolyPtr->myRefCount == 1);
//     //    ASSERT(g.WeylPolyPtr->myRefCount == 1);
//     //    MakeWritable(f);//????????
//     MakeWritable(g);//????????
//     WeylPoly& G = AsWeylPoly(g);
//     WeylPoly::summand* ltg = G.mySummands;
//     G.mySummands = G.mySummands->myNext;
//     if (G.mySummands == nullptr) G.myEnd = &(G.mySummands);
//     PushFront(f, ltg);
//   }


//   void RingExtAlgImpl::DeleteLM(RawPtr rawf) const
//   {
//     ASSERT(!IsZero(f));
//     //    ASSERT(f.WeylPolyPtr->myRefCount == 1);
//     MakeWritable(f);//????????
//     WeylPoly& F = AsWeylPoly(f);
//     WeylPoly::summand* old_ltf = F.mySummands;
//     F.mySummands = old_ltf->myNext;
//     if (F.mySummands == nullptr) F.myEnd = &F.mySummands;
//     old_ltf->myNext = nullptr;
//     DeleteSummands(old_ltf);
//   }


//   void RingExtAlgImpl::AddMul(RawPtr rawf, const WeylPoly::summand* s, ConstRawPtr rawg, bool SkipLM) const
//   {
// ////    std::clog << "AddMul: Doing funny product of the following two polys" << endl;
// ////    output(std::clog, rawg);
//     ASSERT(f.WeylPolyPtr->myRefCount == 1);
//     //    MakeWritable(f);

//     RingElem ppg(*this);
//     assign(raw(ppg), rawg);
//     vector<long> expv(myNumIndets);
//     ordering(myPPM).ComputeExpv(&expv[0], s->myOrdv);
// ////    std::clog << "expv: "; for (int i=0; i<myNumIndets;++i) std::clog << expv[i] << "  "; std::clog << endl;
//     for (long var = myNumIndets/2; var < myNumIndets; ++var)
//     {
//       const long n = expv[var];
//       if (n == 0) continue;
// ////      std::clog << "AddMul: doing D variable with index " << n << endl;
//       RingElem der(*this);
//       assign(raw(der), raw(ppg));
// //      ppg *= IndetPower(var, n);  CANNOT DO THIS --> INFINITE RECURSION
//       {
//         PlainMul(raw(ppg),raw(ppg),raw(IndetPower(var, n)));
//       }
// //      mul(raw(ppg), raw(ppg), raw(IndetPower(var, n)));
//       for (long i=1; i <= n; ++i)
//       {
//         deriv(raw(der), raw(der), var-myNumIndets/2);
// ////        std::clog << "der(" << i << ")="; output(std::clog, raw(der)); std::clog << endl;
// //        ppg += binomial(n, i)*der*IndetPower(var, n-i); // *IndetPower(h, 2*i); // for homog case
//         {
//           vector<long> expv2(myNumIndets);
//           expv2[var] = n-i;
//           RingElem shift(*this);
//           PushBack(raw(shift), raw(RingElem(myCoeffRing, binomial(n, i))), expv2);
//           RingElem tmp(*this);
//           PlainMul(raw(tmp), raw(shift), raw(der));
// ////          std::clog << "AddMul: adding "; output(std::clog, raw(tmp)); std::clog << endl;
//           add(raw(ppg),raw(ppg),raw(tmp));
//         }
//       }
//     }
//     { // f *= var^deg(pp, var); for the early vars
//       for (long var = myNumIndets/2; var < myNumIndets; ++var)
//         expv[var] = 0;
//       PPMonoid::elem qq(myPPM, expv);
//       RingElem c(myCoeffRing);
//       myCoeffRing.assign(raw(c), s->myCoeff);
//       PlainMul(raw(ppg), raw(ppg), raw(monomial(c, qq)));
// ///      f *= P.monomial(RingElem(CoeffRing(P), 1), qq);
//     }
// ////    std::clog << "AddMul: OUTPUT "; output(std::clog, raw(ppg)); std::clog << endl;
//     add(f, f, raw(ppg));
//   }


//   void RingExtAlgImpl::AddMul2(RawPtr rawf, ConstRawPtr rawh, ConstRawPtr rawg, bool SkipLM) const //???? delete me???
//   {                                                 //???
//     AddMul(f, (AsWeylPoly(h).mySummands), rawg, SkipLM);     //???
//   }                                                 //???

//   void RingExtAlgImpl::PlainMul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
//   {
//     if (NumTerms(x) > NumTerms(y)) { PlainMul(lhs, rawy, rawx); return; }

//     RingElem ans(*this);

//     for (const WeylPoly::summand* xterm = AsWeylPoly(x).mySummands; xterm; xterm = xterm->myNext)
//       PlainAddMul(raw(ans), xterm, rawy, false);//?????

//     swap(raw(ans), lhs); // really an assignment
//   }
//   void RingExtAlgImpl::PlainAddMul(RawPtr rawf, const WeylPoly::summand* s, ConstRawPtr rawg, bool SkipLM) const
//   {
// //    ASSERT(f.DMPIPtr->myRefCount == 1);
//     //    MakeWritable(f);
// //          std::clog << "\nF := "; output(std::clog,f);
// //          std::clog<<";\nG := "; output(std::clog,g);
// //          std::clog<<endl;
// //          std::clog<<"s = ";
// //          myCoeffRing.output(std::clog, s->myCoeff);
// //          std::clog<<"x^(";
// //          for(long i=0;i<myOrdvWords;++i)std::clog<<s->myOrdv[i]<<" ";
// //          std::clog<<")"<<endl;
//     const AbstractRing& R = myCoeffRing;
//     typedef WeylPoly::summand summand;
//     const summand* g_smnd = AsWeylPoly(g).mySummands;
//     if (SkipLM)    g_smnd = g_smnd->myNext;
//     summand** f_prev = &(AsWeylPoly(f).mySummands);
//     summand*  f_smnd = *f_prev;

//     summand* tmp_smnd = InitSummand(AllocSummand()); // just sets coeff = 0

//     int CMP=0;

//     //    bool qIsOne = (myPPM.cmp(q, tmpPP)==0);
//     bool qIsOne = false;

//     for (; f_smnd != nullptr && g_smnd != nullptr; g_smnd = g_smnd->myNext)
//     {
//       if (qIsOne)
//         while (f_smnd != nullptr && (CMP=ordering(myPPM).CmpOrdvs(f_smnd->myOrdv,g_smnd->myOrdv)) >0)
//           f_smnd = *(f_prev = &f_smnd->myNext);
//       else
//       {
//         MulOrdv(tmp_smnd->myOrdv, s->myOrdv, g_smnd->myOrdv);
//         while (f_smnd != nullptr && (CMP=ordering(myPPM).CmpOrdvs(f_smnd->myOrdv,tmp_smnd->myOrdv)) >0)
//           f_smnd = *(f_prev = &f_smnd->myNext);
//       }
//       if (f_smnd == nullptr)
//       {
//         R.mul(tmp_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff);
//         PushBack(f, tmp_smnd);
//         tmp_smnd = InitSummand(AllocSummand());
//         g_smnd = g_smnd->myNext;
//         break;
//       }
//       if (CMP == 0)
//       {
//         if (R.IsZeroAddMul(f_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff))
//           RemoveSmnd(f, f_prev);  // f_prev = f_prev;
//         else
//           f_prev = &f_smnd->myNext;
//         f_smnd = *f_prev;
//       }
//       else // (CMP < 0)
//       {
//         R.mul(tmp_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff);
//         InsertSmnd(f, tmp_smnd, f_prev);
//         tmp_smnd = InitSummand(AllocSummand());
//         f_prev = &(*f_prev)->myNext;
//         // f_smnd = f_smnd;
//       }
//     }
//     for (;g_smnd != nullptr; g_smnd = g_smnd->myNext)
//     {
//       MulOrdv(tmp_smnd->myOrdv, s->myOrdv, g_smnd->myOrdv);
//       R.mul(tmp_smnd->myCoeff, s->myCoeff, g_smnd->myCoeff);
//       PushBack(f, tmp_smnd);
//       tmp_smnd = InitSummand(AllocSummand());
//     }
//     DeleteSummands(tmp_smnd); // next ptr will be zero (set by InitSummand)
// //      std::clog << "AddMul: produced f=";output(std::clog,f);std::clog<<endl;
// //      std::clog << "------------------------------------------------------"<<endl;

//   }





//   void RingExtAlgImpl::ReductionStep(RawPtr rawf, ConstRawPtr rawg) const
//   {
//  //             std::clog << "\nRingExtAlgImpl::reduce" << endl;
// //              std::clog << "\nF := "; output(std::clog,f);
// //              std::clog<<";\nG := "; output(std::clog,g);
// //              std::clog<<";"<<endl;
//     ASSERT(&g!=&f);
//     const AbstractRing& R = myCoeffRing;
//     const WeylPoly::summand* g_smnd = AsWeylPoly(g).mySummands;
//     const WeylPoly::summand* f_smnd = AsWeylPoly(f).mySummands;
//     WeylPoly::summand* tmp_smnd = InitSummand(AllocSummand()); // just sets coeff = 0

//     DivOrdv(tmp_smnd->myOrdv, f_smnd->myOrdv, g_smnd->myOrdv);
//     R.div(tmp_smnd->myCoeff, f_smnd->myCoeff, g_smnd->myCoeff);
//     R.negate(tmp_smnd->myCoeff, tmp_smnd->myCoeff);

//  //             std::clog<<"--  S := ";
// //              myCoeffRing.output(std::clog, tmp_smnd->myCoeff);
// //              std::clog<<"x^(";
// //              for(long i=0;i<myOrdvWords;++i)std::clog<<tmp_smnd->myOrdv[i]<<" ";
// //              std::clog<<")"<<endl;

//     DeleteLM(f);
//     AddMul(f, tmp_smnd, g, true /*SkipLM*/);

//     DeleteSummands(tmp_smnd); // next ptr will be zero (set by InitSummand)
//  //             std::clog << "H := "; output(std::clog,f);
// //              std::clog << ";\n--------------------------------------------------"<<endl;
//  }


//   void RingExtAlgImpl::AddClear(RawPtr rawf, RawPtr rawg) const
//   {
//     ASSERT(f.WeylPolyPtr->myRefCount == 1);
//     //    MakeWritable(f);
//     MakeWritable(g);
//     // input polynomial are copied
//     //if (g.WeylPolyPtr->myRefCount != 1)
//     //      std::clog << "AddClear: g.myRefCount == " << g.WeylPolyPtr->myRefCount << endl;
//     const AbstractRing& R = myCoeffRing;
//     typedef WeylPoly::summand summand;
//     WeylPoly& F = AsWeylPoly(f);
//     WeylPoly& G = AsWeylPoly(g);
//     summand*  g_smnd = G.mySummands;
//     summand** f_prev = &(F.mySummands);
//     summand*  f_smnd = *f_prev;
//     int CMP=0;
//     ASSERT(*(G.myEnd)==nullptr);//BUG HUNTING  ???

//     //    std::clog << "input f = "; output(std::clog, f) ;std::clog << endl;
//     while ( f_smnd!=nullptr && g_smnd!=nullptr )
//     {
//       while (f_smnd!=nullptr &&
//              (CMP=ordering(myPPM).CmpOrdvs(f_smnd->myOrdv,g_smnd->myOrdv)) >0)
//         f_smnd = *(f_prev = &f_smnd->myNext);
//       if (f_smnd == nullptr)  break;
//       //std::clog <<   "(AddClear error: should never happen for Basic Reduction)" << endl;
//       G.mySummands = G.mySummands->myNext;
//       g_smnd->myNext = nullptr;
//       if (CMP == 0)
//       {
//         R.add(f_smnd->myCoeff, f_smnd->myCoeff, g_smnd->myCoeff);
//         if (R.IsZero(f_smnd->myCoeff))
//           RemoveSmnd(f, f_prev);
//         DeleteSummands(g_smnd);
//       }
//       else // (CMP < 0)
//       {
//         InsertSmnd(f, g_smnd, f_prev);
//         f_prev = &(*f_prev)->myNext;
//       }
//       f_smnd = *f_prev;
//       g_smnd = G.mySummands;
//     }
//     if (G.mySummands!=nullptr)
//     {
//       *(F.myEnd) = G.mySummands;
//       F.myEnd = G.myEnd;
//       G.mySummands = nullptr;
//     }
//     G.myEnd = &G.mySummands;
//     //    if (rare) {std::clog << "f2 = "; output(std::clog, f) ;std::clog << endl;}
//   }


//   void RingExtAlgImpl::AppendClear(RawPtr rawf, RawPtr rawg) const
//   {
//     ASSERT(f.WeylPolyPtr->myRefCount == 1);
//     //    MakeWritable(f);
//     MakeWritable(g);
//     // input polynomial are copied
//     //if (g.WeylPolyPtr->myRefCount != 1)
//     //      std::clog << "AppendClear: g.myRefCount == " << g.WeylPolyPtr->myRefCount << endl;
//     WeylPoly& F = AsWeylPoly(f);
//     WeylPoly& G = AsWeylPoly(g);
//     if (G.mySummands!=nullptr)
//     {
//       *(F.myEnd) = G.mySummands;
//       F.myEnd = G.myEnd;
//       G.mySummands = nullptr;
//     }
//     G.myEnd = &G.mySummands;
//     //    if (rare) {std::clog << "f2 = "; output(std::clog, f) ;std::clog << endl;}
//   }


//   int  RingExtAlgImpl::CmpLPP(ConstRawPtr rawf, ConstRawPtr rawg) const
//   {
//     ASSERT(!IsZero(f));
//     ASSERT(!IsZero(g));
//     return ordering(myPPM).CmpOrdvs(AsWeylPoly(f).mySummands->myOrdv,AsWeylPoly(g).mySummands->myOrdv);
//   }


//   void RingExtAlgImpl::DivLM(RawPtr rawlhs, ConstRawPtr rawf, ConstRawPtr rawg) const
//   {
//     //    std::clog << "DivLM" << endl;
//     const AbstractRing& R = myCoeffRing;
//     typedef WeylPoly::summand summand;
//     const summand* f_smnd = AsWeylPoly(f).mySummands;
//     const summand* g_smnd = AsWeylPoly(g).mySummands;
//     MakeWritable(lhs);
//     assign(lhs,0);
//     summand* SpareSummand = InitSummand(AllocSummand());
//     R.div(SpareSummand->myCoeff, f_smnd->myCoeff, g_smnd->myCoeff);
//     DivOrdv(SpareSummand->myOrdv, f_smnd->myOrdv, g_smnd->myOrdv);
//     PushBack(lhs, SpareSummand);
//   }


//   void RingExtAlgImpl::mul(RawPtr rawf, const PPMonoid::elem& pp) const
//   {
//     //    bool qIsOne = (myPPM.cmp(q, tmpPP)==0);
//     MakeWritable(f);
//     std::vector<long> v(myNumIndets);
//     myPPM.exponents(v, raw(pp));

//     WeylPoly::summand* f_smnd = AsWeylPoly(f).mySummands;
//     WeylPoly::summand* s = InitSummand(AllocSummand());

//     ordering(myPPM).ComputeOrdv(s->myOrdv, &v[0]); // &v[0] converts vector<T> to T*

//     for (; f_smnd != nullptr ; f_smnd = f_smnd->myNext)
//       MulOrdv(f_smnd->myOrdv, f_smnd->myOrdv, s->myOrdv);
//   }


//   void RingExtAlgImpl::PushFront(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
//   {
//     if (CoeffRing().IsZero(c)) return;
//     MakeWritable(f);
//     WeylPoly::summand* t = AllocSummand();
//     InitSummand(t, c, expv);
//     t->myNext = f.WeylPolyPtr->mySummands;
//     if (f.WeylPolyPtr->mySummands == nullptr) f.WeylPolyPtr->myEnd = &t->myNext;
//     f.WeylPolyPtr->mySummands = t;
//   }


//   void RingExtAlgImpl::PushBack(RawPtr rawf, ConstRawPtr rawc, const std::vector<long>& expv) const
//   {
//     if (CoeffRing().IsZero(c)) return;
//     MakeWritable(f);
//     WeylPoly::summand* t = AllocSummand();
//     InitSummand(t, c, expv);
//     *(f.WeylPolyPtr->myEnd) = t;
//     f.WeylPolyPtr->myEnd = &t->myNext;
//   }


//   void RingExtAlgImpl::PushFront(RawPtr rawx, WeylPoly::summand* t) const
//   {
//     MakeWritable(x);
//     WeylPoly& f = AsWeylPoly(x);
//     t->myNext = f.mySummands;
//     f.mySummands = t;
//     if (f.myEnd == &f.mySummands) f.myEnd = &t->myNext;
//   }


//   void RingExtAlgImpl::PushBack(RawPtr rawx, WeylPoly::summand* t) const
//   {
//     MakeWritable(x);
//     WeylPoly& f = AsWeylPoly(x);
//     *f.myEnd = t;
//     f.myEnd = &t->myNext;
//   }


//   void RingExtAlgImpl::RemoveSmnd(RawPtr rawx, WeylPoly::summand** prev_link) const
//   {
//     ASSERT(x.WeylPolyPtr->myRefCount == 1);
//     //    MakeWritable(x);
//     WeylPoly::summand* tmp = *prev_link;
//     ASSERT(tmp != 0);
//     if (tmp->myNext==nullptr) // f.myEnd == &(tmp->myNext)
//     {
//       WeylPoly& f = AsWeylPoly(x);
//       f.myEnd = prev_link;
//     }
//     *prev_link = tmp->myNext;
//     tmp->myNext = nullptr;
//     DeleteSummands(tmp);
//   }


//   void RingExtAlgImpl::InsertSmnd(RawPtr rawx, WeylPoly::summand* s, WeylPoly::summand** prev_link) const
//   {
//     ASSERT(x.WeylPolyPtr->myRefCount == 1);
//     //    MakeWritable(x);
//     WeylPoly& f = AsWeylPoly(x);
//     s->myNext = (*prev_link);
//     (*prev_link) = s;
//     if (f.myEnd == prev_link) f.myEnd = &(s->myNext);
//   }


//   PPMonoid::elem RingExtAlgImpl::LPP(ConstRawPtr rawf) const
//   {
//     ASSERT(!IsZero(f));
//     return PPMonoid::elem(PPM(), PPMonoid::elem::FromOrdv, AsWeylPoly(f).mySummands->myOrdv);
//   }


//   bool RingExtAlgImpl::IsZeroAddLCs(RawPtr rawf, RawPtr rawg) const
//   {
//     ASSERT(!IsZero(f) && !IsZero(g));
//     ASSERT( CmpLPP(f,g) == 0);
//     WeylPoly& F = AsWeylPoly(f);
//     WeylPoly& G = AsWeylPoly(g);
//     ASSERT(F.myRefCount==1 && G.myRefCount==1);
//     myCoeffRing.add(F.mySummands->myCoeff, F.mySummands->myCoeff, G.mySummands->myCoeff);
//     DeleteLM(g);
//     if (!myCoeffRing.IsZero(F.mySummands->myCoeff)) return false;
//     DeleteLM(f);
//     return true;
//   }


//   void RingExtAlgImpl::deriv(RawPtr rawdest, ConstRawPtr rawf, long var) const
//   {
//     const WeylPoly& F = AsWeylPoly(f);
//     RingElem ans(*this);
//     vector<long> expv(myNumIndets);
//     for (WeylPoly::summand* i=F.mySummands; i; i = i->myNext)
//     {
//       ordering(myPPM).ComputeExpv(&expv[0], i->myOrdv);
//       if (expv[var] == 0) continue;
//       RingElem c(myCoeffRing, expv[var]);
//       if (CoCoA::IsZero(c)) continue;
//       myCoeffRing.mul(raw(c), raw(c), i->myCoeff);
//       if (CoCoA::IsZero(c)) continue;
//       --expv[var];
//       PushBack(raw(ans), raw(c), expv);
//     }
//     swap(raw(ans), dest);
//   }


//   void RingExtAlgImpl::negate(RawPtr rawlhs, ConstRawPtr rawx) const
//   {
//     if (lhs.WeylPolyPtr == x.WeylPolyPtr)
//     {
//       MakeWritable(lhs);
//       typedef WeylPoly::summand summand;
//       //      std::clog << "-- negate"; output(std::clog, lhs); std::clog << endl;
//       for (summand* smnd = AsWeylPoly(lhs).mySummands; smnd!=nullptr; smnd=smnd->myNext )
//         myCoeffRing.negate(smnd->myCoeff, smnd->myCoeff);
//       //      std::clog << "-> negate"; output(std::clog, lhs); std::clog << endl;
//     }
//     else
//       std::clog << "RingExtAlgImpl::negate: not yet implemented" << endl;
//   }


//   void RingExtAlgImpl::add(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
//   {
//     const AbstractRing& R = myCoeffRing;
//     RawValue ans;
//     //    ans.WeylPolyPtr = static_cast<WeylPoly*>(myWeylPolyPool.alloc(sizeof(WeylPoly)));
//     init(ans);
//     ///    WeylPoly& sum = AsWeylPoly(ans);
//     typedef WeylPoly::summand summand;
//     const summand* gterm = AsWeylPoly(x).mySummands;
//     const summand* hterm = AsWeylPoly(y).mySummands;
//     summand* SpareSummand = InitSummand(AllocSummand()); // just sets coeff = 0
//     while (gterm != null && hterm != nullptr)
//     {
//       const int cmp = ordering(myPPM).CmpOrdvs(gterm->myOrdv, hterm->myOrdv);

//       if (cmp < 0)
//       {
// 	summand* hcopy = CopySummand(hterm);
// 	PushBack(ans, hcopy);
// 	hterm = hterm->myNext;
// 	continue;
//       }

//       if (cmp > 0)
//       {
// 	summand* gcopy = CopySummand(gterm);
// 	PushBack(ans, gcopy);
// 	gterm = gterm->myNext;
// 	continue;
//       }

//       // Must have cmp == 0 here.
//       // The leading PPs are the same, so we must sum the coeffs.
//       R.add(SpareSummand->myCoeff, gterm->myCoeff, hterm->myCoeff);
//       if (!R.IsZero(SpareSummand->myCoeff))
//       {
// 	SetSummandMOrdv(SpareSummand, gterm); // set module posn and PP
// 	PushBack(ans, SpareSummand);
// 	SpareSummand = InitSummand(AllocSummand());// just sets coeff = 0
//       }
//       gterm = gterm->myNext;
//       hterm = hterm->myNext;
//     }
//     while (gterm != nullptr)
//     {
//       summand* gcopy = CopySummand(gterm);
//       PushBack(ans, gcopy);
//       gterm = gterm->myNext;
//     }
//     while (hterm != nullptr)
//     {
//       summand* hcopy = CopySummand(hterm);
//       PushBack(ans, hcopy);
//       hterm = hterm->myNext;
//     }
//     DeleteSummands(SpareSummand); // next ptr will be zero (set by InitSummand)
//     swap(lhs, ans); // really an assignment
//     kill(ans);
//   }


//   void RingExtAlgImpl::sub(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
//   {
//     // This code copied from RingExtAlgImpl::add...

//     const AbstractRing& R = myCoeffRing;
//     RawValue ans;
//     //    ans.WeylPolyPtr = static_cast<WeylPoly*>(myWeylPolyPool.alloc(sizeof(WeylPoly)));
//     init(ans);
//     ///    WeylPoly& sum = AsWeylPoly(ans);
//     typedef WeylPoly::summand summand;
//     const summand* gterm = AsWeylPoly(x).mySummands;
//     const summand* hterm = AsWeylPoly(y).mySummands;
//     summand* SpareSummand = InitSummand(AllocSummand()); // just sets coeff = 0
//     while (gterm != nullptr && hterm != nullptr)
//     {
//       const int ord = ordering(myPPM).CmpOrdvs(gterm->myOrdv, hterm->myOrdv);

//       if (ord < 0)
//       {
// 	summand* hcopy = CopySummand(hterm);
// 	R.negate(hcopy->myCoeff, hcopy->myCoeff);
// 	//	sum.push_front(hcopy);
// 	PushBack(ans, hcopy);
// 	hterm = hterm->myNext;
// 	continue;
//       }

//       if (ord > 0)
//       {
// 	summand* gcopy = CopySummand(gterm);
// 	//	sum.push_front(gcopy);
// 	PushBack(ans, gcopy);
// 	gterm = gterm->myNext;
// 	continue;
//       }

//       // The leading PPs are the same, so we must sum the coeffs.
//       R.sub(SpareSummand->myCoeff, gterm->myCoeff, hterm->myCoeff);
//       if (!R.IsZero(SpareSummand->myCoeff))
//       {
// 	SetSummandMOrdv(SpareSummand, gterm); // set module posn and PP
// 	//	sum.push_front(SpareSummand);
// 	PushBack(ans, SpareSummand);
// 	SpareSummand = InitSummand(AllocSummand());// just sets coeff = 0
//       }
//       gterm = gterm->myNext;
//       hterm = hterm->myNext;
//     }
//     while (gterm != nullptr)
//     {
//       summand* gcopy = CopySummand(gterm);
//       //      sum.push_front(gcopy);
//       PushBack(ans, gcopy);
//       gterm = gterm->myNext;
//     }
//     while (hterm != nullptr)
//     {
//       summand* hcopy = CopySummand(hterm);
//       //      sum.push_front(hcopy);
//       R.negate(hcopy->myCoeff, hcopy->myCoeff);
//       PushBack(ans, hcopy);
//       hterm = hterm->myNext;
//     }
//     DeleteSummands(SpareSummand); // next ptr will be zero (set by InitSummand)
//     /// NO LONGER NEEDED    sum.reverse();
//     swap(lhs, ans); // really an assignment
//     kill(ans);
//   }


//   void RingExtAlgImpl::mul(RawPtr rawlhs, ConstRawPtr rawx, ConstRawPtr rawy) const
//   {
// // NO!! NOT COMMUTATIVE!!    if (NumTerms(x) > NumTerms(y)) { mul(lhs, y, x); return; }

// ////    std::clog << "MUL on "; output(std::clog, x); std::clog << " and "; output(std::clog, y); std::clog << endl;
//     RingElem ans(*this);

//     for (const WeylPoly::summand* xterm = AsWeylPoly(x).mySummands; xterm; xterm = xterm->myNext)
//       AddMul(raw(ans), xterm, y, false);//?????

//     swap(raw(ans), lhs); // really an assignment
//   }
