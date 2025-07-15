//   Copyright (c)  2005-2013  John Abbott, Anna M. Bigatti

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


// Implementation file for the class SubmoduleImpl

#include "CoCoA/submodule.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseMatrix.H" // for GensAsRows, GensAsCols
#include "CoCoA/FreeModule.H"
#include "CoCoA/MemPool.H"
#include "CoCoA/OpenMath.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/TmpGOperations.H"  // for ComputeGBasis
#include "CoCoA/VectorOps.H"  // for HasUniqueOwner
#include "CoCoA/ideal.H" // for syzygies
#include "CoCoA/matrix.H" // for ConstMatrixView
#include "CoCoA/ring.H"
//#include "CoCoA/MatrixView.H"

#include <iostream>
using std::ostream;
//#include <vector>
using std::vector;


namespace CoCoA
{

  class SubmoduleImpl: public FGModuleBase
  {
    // Two typedefs to save typing.
    typedef ModuleBase::RawPtr RawPtr;
    typedef const ModuleBase::RawPtr& ConstRawPtr;

  public:
    SubmoduleImpl(const module& M, const std::vector<ModuleElem>& gens);
    long myNumCompts() const override  {return NumCompts(myM);}
    const ring& myRing() const override  {return RingOf(myM);}
    const FreeModule& myAmbientFreeModule() const override  {return myM;}
    const std::vector<ModuleElem>& myGens() const override  {return myGensValue;}
    const std::vector<ModuleElem>& myMinGens(const CpuTimeLimit& CheckForTimeout) const override;
    const std::vector<ModuleElem>& myTidyGens(const CpuTimeLimit& CheckForTimeout) const override;
    const std::vector<ModuleElem>& myGBasis(const CpuTimeLimit& CheckForTimeout) const; // for SparsePolyRing
    
    const ModuleElem& myZero() const override  {return zero(myM);}
    void myNew(RawPtr& rawv) const override    {myM->myNew(rawv);}
    void myNew(RawPtr& rawv, ConstRawPtr rawt) const override  {myM->myNew(rawv, rawt);}
    void myDelete(RawPtr& rawv) const override  {myM->myDelete(rawv);}  // destroys v (incl all resources)
    void mySwap(RawPtr& rawv, RawPtr& raww) const override  {myM->mySwap(rawv, raww);}
    void myAssign(RawPtr& rawlhs, ConstRawPtr rawv) const override  {myM->myAssign(rawlhs, rawv);} // lhs = v;
    ConstRefRingElem myCompt(const RawPtr& rawv, long pos) const override;            ///< v[pos] (READ ONLY)
    void myNegate(RawPtr& rawlhs, ConstRawPtr rawv) const override  {myM->myNegate(rawlhs, rawv);} // lhs = -v
    void myAdd(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override  {myM->myAdd(rawlhs, rawv, raww);} // lhs = v+w;
    void mySub(RawPtr& rawlhs, ConstRawPtr rawv, ConstRawPtr raww) const override  {myM->mySub(rawlhs, rawv, raww);} // lhs = v-w;

    void myMul(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override  {myM->myMul(rawlhs, rawx, rawv);} // lhs = r*v;
    void myDiv(RawPtr& rawlhs, RingElemConstRawPtr rawx, ConstRawPtr rawv) const override  {myM->myDiv(rawlhs, rawx, rawv);} // lhs = (1/r)*v;
    void myOutput(std::ostream& out, ConstRawPtr rawv) const override  {myM->myOutput(out, rawv);} // out << v
    void myOutputSelf(std::ostream& out) const override;                   // out << M
    void myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const override; // OMOut << v
    void myOutputSelf_OM(OpenMathOutput& OMOut) const override;               // OMOut << M
    bool myIsZero(ConstRawPtr rawv) const override  {return myM->myIsZero(rawv);} // v == 0
//???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.
    bool myIsEqual(ConstRawPtr rawx, ConstRawPtr rawy) const override  {return myM->myIsEqual(rawx, rawy);}

  private: // data members
    const FreeModule myM;
    std::vector<ModuleElem> myGensValue;
    // mutable member fields
    mutable bool myTidyGensIsValid;
    mutable std::vector<ModuleElem> myMinGensValue;
    mutable std::vector<ModuleElem> myTidyGensValue;
//???    std::vector<ModuleElem>& ComputeTidyGens() const;
  };



  SubmoduleImpl::SubmoduleImpl(const module& M, const std::vector<ModuleElem>& gens):
      myM(M),
      myGensValue(gens),
      myTidyGensIsValid(false)
  {
    CoCoA_ASSERT(IsFreeModule(M));
    for (long i=0; i < len(gens); ++i)
      if (owner(gens[i]) != M)  CoCoA_THROW_ERROR1(ERR::MixedModules);
    myRefCountZero();
  }


//   ideal::ideal(const std::vector<RingElem>& gens)
//   {
//     if (gens.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
//     if (!HasUniqueOwner(gens)) CoCoA_THROW_ERROR1(ERR::MixedRings);
//     ideal tmp = owner(gens[0])->myIdealCtor(gens);
//     myPtr = tmp.myPtr;
//     myPtr->myRefCountInc();
//   }


  const std::vector<ModuleElem>& SubmoduleImpl::myGBasis(const CpuTimeLimit& CheckForTimeout) const
  {
    CoCoA_ASSERT(IsSparsePolyRing(myRing()));
    if (!IsField(CoeffRing(myRing())))  CoCoA_THROW_ERROR2(ERR::NYI, "coeffs not in a field");

    //    if (IhaveMonomialGens()) return myGBasisMonId();
    if (myTidyGensIsValid) return myTidyGensValue;
    CoCoA_ASSERT(myTidyGensValue.empty());
    //    if (IamZero()) return myTidyGensValue;
    ComputeGBasis(myTidyGensValue, myMinGensValue, myGens(), CheckForTimeout);
    myTidyGensIsValid = true;
    return myTidyGensValue;
  }


  const std::vector<ModuleElem>& SubmoduleImpl::myMinGens(const CpuTimeLimit& CheckForTimeout) const
  {
    if (!myMinGensValue.empty()) return myMinGensValue;
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    if (!HasPositiveGrading(myRing()))  CoCoA_THROW_ERROR1(ERR::ReqPositiveGrading);
    myGBasis(CheckForTimeout);
    return myMinGensValue;
  }


  const std::vector<ModuleElem>& SubmoduleImpl::myTidyGens(const CpuTimeLimit& CheckForTimeout) const
  {
    //    if (!myTidyGensIsValid)
    if (!IsSparsePolyRing(myRing()))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of module");
    return myGBasis(CheckForTimeout);
  }


  ConstRefRingElem SubmoduleImpl::myCompt(const RawPtr& rawv, long pos) const
  {
    CoCoA_ASSERT(0 <= pos && pos < myNumCompts());
    return myM->myCompt(rawv, pos);
  }


  namespace{  // anonymous
  //??? the following functions to compute NR will be replaced by GBMill

    int FindReducerIndex(ConstRefPPMonoidElem pp, long posn, const vector<ModuleElem>& g)
  {
    const long nelems = len(g);
    long posn_gi;
    for (long i=0; i < nelems; ++i)
      if (posn == (posn_gi=LPosn(g[i])))
        if (IsDivisible(pp, LPP(g[i][posn_gi])))
          return i;
    return -1;
  }


  inline int FindReducerIndex(const ModuleElem& F, const vector<ModuleElem>& g)
  {
    if ( IsZero(F) ) return -1;
    return FindReducerIndex(LPP(F), LPosn(F), g);
  }


  void ReduceLM(ModuleElem& F, const vector<ModuleElem>& g)
  {
    long i;
    while ( (i = FindReducerIndex(F, g) ) != -1)
    {
      const long pg = LPosn(g[i]);
      const long pF = pg; /// should be equal to LPosn(F);  assert???
      const auto q = monomial(RingOf(owner(F)), 
                              LC(F[pF])/LC(g[i][pg]), LPP(F[pF])/LPP(g[i][pg]));
      F -= q*g[i]; // mult on LEFT!!
    }
  }


  void reduce(ModuleElem& F, const vector<ModuleElem>& g)
  {
    ReduceLM(F, g);
//     while ( !IsActiveZero(F) )
//     {
//       F->myMoveToNextLM();
//       ReduceActiveLM(F, v);
//     }
  }

  ModuleElem NR(const ModuleElem& f, const vector<ModuleElem>& g)
  {
    if (!IsField(CoeffRing(RingOf(owner(f)))))
      CoCoA_THROW_ERROR2(ERR::NYI, "coeffs not in a field");//???
    if ( IsZero(f) ) return f;
    ModuleElem F=f;
    reduce(F, g);
    return F;
  }

  } // anonymous namespace


  bool IsElem(const ModuleElem& v, const module& M)
  {
    CoCoA_ASSERT(IsFGModule(M));
    if (owner(v) != AmbientFreeModule(M))  CoCoA_THROW_ERROR1(ERR::MixedModules);
    //    return I->IhaveElem(raw(r));
    //??? for FGmodule only 
    return IsZero(NR(v, TidyGens(M)));
  }


  void SubmoduleImpl::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams

    out << "submodule(" << myM << ", [";
    if (!myGensValue.empty()) out << myGensValue[0];
    for (long i=1; i < len(myGensValue); ++i)
    {
      out << ", " << myGensValue[i];
    }
    out << "])";
  }


  void SubmoduleImpl::myOutput_OM(OpenMathOutput& OMOut, ConstRawPtr rawv) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "ModuleElement"); // BUG: what should this OMSymbol be???
    OMOut << myM;
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "list"); // BUG: what should this OMSymbol be???
    OMOut << myNumCompts();
    myM->myOutput_OM(OMOut, rawv); // BUG: this should be a "naked" output???
    OMOut->mySendApplyEnd();
    OMOut->mySendApplyEnd();
  }


  void SubmoduleImpl::myOutputSelf_OM(OpenMathOutput& OMOut) const
  {
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "submodule"); // BUG: what should this OMSymbol be???
    OMOut << myM;
    OMOut->mySendApplyStart();
    OMOut << OpenMathSymbol("???", "list"); // BUG: what should this OMSymbol be???
    OMOut << len(myGensValue);
    for (long i=0; i < len(myGensValue); ++i)
      OMOut << myGensValue[i];  // BUG: this should be a "naked" output???
    OMOut->mySendApplyEnd();
    OMOut->mySendApplyEnd();
  }


 //???    bool IsZeroAddMul(RawPtr& rawlhs, RingElemConstRawPtr rawy, ConstRawPtr rawz) const;  // lhs += y*z, result says whether lhs == 0.

  //-- pseudo-ctors

  FGModule submodule(const FGModule& M, const std::vector<ModuleElem>& gens)
  { return FGModule(new SubmoduleImpl(M, gens)); }


  FGModule submodule(const std::vector<ModuleElem>& gens)
  {
    if (gens.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(gens))  CoCoA_THROW_ERROR1(ERR::MixedModules);
    return submodule(owner(gens[0]), gens);
  }


  FGModule submodule(const ModuleElem& v1)
  {
    vector<ModuleElem> gens;
    gens.push_back(v1);
    return submodule(gens);
  }
  

  FGModule submodule(const ModuleElem& v1, const ModuleElem& v2)
  {
    vector<ModuleElem> gens;
    gens.push_back(v1);gens.push_back(v2);
    return submodule(gens);
  }
  

  FGModule submodule(const ModuleElem& v1, const ModuleElem& v2, const ModuleElem& v3)
  {
    vector<ModuleElem> gens;
    gens.push_back(v1);gens.push_back(v2);gens.push_back(v3);
    return submodule(gens);
  }
  

  FGModule submodule(const ModuleElem& v1, const ModuleElem& v2, const ModuleElem& v3, const ModuleElem& v4)
  {
    vector<ModuleElem> gens;
    gens.push_back(v1);gens.push_back(v2);gens.push_back(v3);gens.push_back(v4);
    return submodule(gens);
  }
  

  FGModule SubmoduleCols(const FGModule& F, ConstMatrixView M)
  {
    const std::vector<ModuleElem>& e = gens(F);
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    if (len(e) != nrows)  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    std::vector<ModuleElem> g(ncols, zero(F));
    for (long i=0; i<nrows; ++i)
      for (long j=0; j<ncols; ++j)
        g[j] += M(i,j) * e[i];
    return FGModule(new SubmoduleImpl(F, g));
  }
  

  FGModule SubmoduleRows(const FGModule& F, ConstMatrixView M)
  {
    const std::vector<ModuleElem>& e = gens(F);
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    if (len(e) != ncols)  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    std::vector<ModuleElem> g(nrows, zero(F));
    for (long i=0; i<nrows; ++i)
      for (long j=0; j<ncols; ++j)
        g[i] += M(i,j) * e[j];
    return submodule(F, g);
  }


  FGModule SubmoduleOfMinGens(const FGModule& F)
  {
    if (IsFreeModule(F)) return F;
    return submodule(AmbientFreeModule(F), MinGens(F));
  }
  

  matrix GensAsRows(const FGModule& Mod)
  {
    const std::vector<ModuleElem>& g = gens(Mod);
    matrix M = NewDenseMat(RingOf(Mod), len(g), NumCompts(Mod));
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i=0; i<nrows; ++i)
      for (long j=0; j<ncols; ++j)
        SetEntry(M,i,j, g[i][j]);
    return M;
  }
  

  matrix GensAsCols(const FGModule& Mod)
  {
    const std::vector<ModuleElem>& g = gens(Mod);
    matrix M = NewDenseMat(RingOf(Mod), NumCompts(Mod), len(g));
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i=0; i<nrows; ++i)
      for (long j=0; j<ncols; ++j)
        SetEntry(M,i,j, g[j][i]);
    return M;
  }
  

  namespace // anonymous
  {
    // template? ???

    ModuleElem InsertZeroes(const FreeModule& F, const vector<RingElem>& L, const ModuleElem& v)
    {
      const std::vector<ModuleElem>& e = gens(F);
      ModuleElem w(F);
      long j=0;
      for (long i=0; i<len(L); ++i)
        if (IsZero(L[i]))   ++j;  else  w += v[i-j]*e[i];
      return w;
    }

    ModuleElem InsertZeroes(const FreeModule& F, const vector<ModuleElem>& L, const ModuleElem& v)
    {
      const std::vector<ModuleElem>& e = gens(F);
      ModuleElem w(F);
      long j=0;
      for (long i=0; i<len(L); ++i)
        if (IsZero(L[i]))   ++j;  else  w += v[i-j]*e[i];
      return w;
    }

  } // anonymous
  

  FGModule syz(const std::vector<RingElem>& g)
  { return syz(NewFreeModuleForSyz(g, CoCoA_ERROR_CONTEXT), g); }


  FGModule syz(const FreeModule& F, const std::vector<RingElem>& g)
  {
    if (g.empty())  CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!IsField(CoeffRing(RingOf(F))))  CoCoA_THROW_ERROR2(ERR::NYI, "coeffs not in a field");//???
    if (NumCompts(F)!=len(g))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    vector<RingElem> g_non0;
    //-- remove zero gens
    for (long i=0; i<len(g); ++i)
    {
      if ( (!g_non0.empty()) && owner(g[i])!=owner(g_non0[0]))  CoCoA_THROW_ERROR1(ERR::MixedRings);
      if (!IsZero(g[i]))
        g_non0.push_back(g[i]);
    }
    vector<ModuleElem> SyzVec_non0;
    ComputeSyz(SyzVec_non0, F, g_non0);
    if (len(g_non0)==len(g))  return submodule(F, SyzVec_non0);
    const std::vector<ModuleElem>& e = gens(F);
    //-- syzygies of zero gens
    vector<ModuleElem> SyzVec;
    for (long i=0; i<len(g); ++i)
      if (IsZero(g[i]))  SyzVec.push_back(e[i]);
    //-- interweave 0 in SyzVec_non0
    for (long i=0; i<len(SyzVec_non0); ++i)
      SyzVec.push_back(InsertZeroes(F, g, SyzVec_non0[i]));
    return submodule(SyzVec);
  }


  FGModule SyzOfGens(const FreeModule& F, const ideal& I)
  { return syz(F, gens(I)); }


  FGModule SyzOfGens(const FreeModule& F, const FGModule& N)
  {
    if (!IsField(CoeffRing(RingOf(F))))  CoCoA_THROW_ERROR2(ERR::NYI, "coeffs not in a field");//???
    if (NumCompts(F)!=len(gens(N)))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    const vector<ModuleElem>& g = gens(N);
    vector<ModuleElem> g_non0;
    //-- remove zero gens
    for (long i=0; i<len(g); ++i)
      if (!IsZero(g[i]))  g_non0.push_back(g[i]);
    vector<ModuleElem> SyzVec_non0;
    ComputeSyz(SyzVec_non0, F, g_non0);
    if (len(g_non0)==len(g))  return submodule(F, SyzVec_non0);
    const std::vector<ModuleElem>& e = gens(F);
    //-- syzygies of zero gens
    vector<ModuleElem> SyzVec;
    for (long i=0; i<len(g); ++i)
      if (IsZero(g[i]))  SyzVec.push_back(e[i]);
    //-- interweave 0 in SyzVec_non0
    for (long i=0; i<len(SyzVec_non0); ++i)
      SyzVec.push_back(InsertZeroes(F, g, SyzVec_non0[i]));
    return submodule(SyzVec);
  }
  

  FGModule SyzOfGens(const ideal& I)
  {
    if (!IsField(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR2(ERR::NYI, "coeffs not in a field");//???
    const vector<RingElem>& g = gens(I);
    for (long i=0; i<len(g); ++i)
      if (IsZero(g[i]))
        CoCoA_THROW_ERROR2(ERR::BadArg, "with zero generator(s) first arg must be the output free module");
    vector<ModuleElem> syzvec;
    FreeModule F = NewFreeModuleForSyz(gens(I), CoCoA_ERROR_CONTEXT);
    return SyzOfGens(F, I);
  }
  

  FGModule SyzOfGens(const FGModule& N)
  {
    if (!IsField(CoeffRing(RingOf(N))))
      CoCoA_THROW_ERROR2(ERR::NYI, "coeffs not in a field");//???
    const vector<ModuleElem>& g = gens(N);
    for (long i=0; i<len(g); ++i)
      if (IsZero(g[i]))
        CoCoA_THROW_ERROR2(ERR::BadArg, "with zero generator(s) first arg must be the output free module");
    vector<ModuleElem> syzvec;
    FreeModule F = NewFreeModuleForSyz(gens(N), CoCoA_ERROR_CONTEXT);
    return SyzOfGens(F, N);
  }
  

  bool IsContained(const module& M, const module& N)
  {
    if (!IsSparsePolyRing(RingOf(M)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    const FGModule M_fg(M, CoCoA_ERROR_CONTEXT);
    const vector<ModuleElem>& g = gens(M_fg);
    for (long i=0; i < len(g); ++i)
      if (!IsElem(g[i], N)) return false;
    return true;
  }


  bool IsHomog(const module& M)
  {
    const FGModule M_fg(M, CoCoA_ERROR_CONTEXT);
    const SparsePolyRing P(RingOf(M_fg), CoCoA_ERROR_CONTEXT);
////    if (!IsSparsePolyRing(RingOf(M)))
////      CoCoA_THROW_ERROR1(ERR::NotSparsePolyRing); 
    if (GradingDim(P)==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (IsZero(M_fg))  return true;
    // Now we know I is non-trivial.
    const vector<ModuleElem>& GB = TidyGens(M);
    const long GBsize = len(GB); // MUST BE A REDUCED GBASIS !!!
    for (long i=0; i < GBsize; ++i)
      if (!IsHomog(GB[i]))  return false;
    return true;
  }


  // intersection


  FGModule LT(const module& M)
  {
    const FGModule M_fg(M, CoCoA_ERROR_CONTEXT);
    const SparsePolyRing P(RingOf(M_fg), CoCoA_ERROR_CONTEXT);
//// if (!IsSparsePolyRing(RingOf(M)))  CoCoA_THROW_ERROR1(ERR::NotSparsePolyRing);
////    const SparsePolyRing P = RingOf(M); 
    if (GradingDim(P)==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    const vector<ModuleElem>& e = gens(AmbientFreeModule(M));
    const vector<ModuleElem>& GB = TidyGens(M);
    vector<ModuleElem> LTs;
    const long GBsize = len(GB);
    for (long i=0; i < GBsize; ++i)
    {
      long j = LPosn(GB[i]);
      LTs.push_back(monomial(P, LPP(GB[i][j]))*e[j]);
    }
    return submodule(AmbientFreeModule(M), LTs);
  }


//   FGModule SubmoduleOfGBasis(const module& M)
//   {
//     ideal J(RingOf(I), GBasis(I));
//     SetGBasisAsGens(J);
//     return J;
//   }


}  // end of namespace CoCoA
