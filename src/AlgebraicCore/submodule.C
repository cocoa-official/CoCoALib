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
      if (owner(gens[i]) != M)
        CoCoA_THROW_ERROR(ERR::MixedModules, "SubmoduleImpl(M, gens)");
    myRefCountZero();
  }


//   ideal::ideal(const std::vector<RingElem>& gens)
//   {
//     if (gens.empty()) CoCoA_THROW_ERROR("Empty list of generators: need ring", "ideal(gens)");
//     if (!HasUniqueOwner(gens)) CoCoA_THROW_ERROR(ERR::MixedRings, "ideal(gens)");
//     ideal tmp = owner(gens[0])->myIdealCtor(gens);
//     myPtr = tmp.myPtr;
//     myPtr->myRefCountInc();
//   }


  const std::vector<ModuleElem>& SubmoduleImpl::myGBasis(const CpuTimeLimit& CheckForTimeout) const
  {
    CoCoA_ASSERT(IsSparsePolyRing(myRing()));
    if (!IsField(CoeffRing(myRing())))
      CoCoA_THROW_ERROR("ERR:NYI coeffs not in a field", "SubmoduleImpl::myGBasis");//???

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
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))
      CoCoA_THROW_ERROR("Input is not homogeneous", "SubmoduleImpl::myMinGens");
    if (!HasPositiveGrading(myRing()))
      CoCoA_THROW_ERROR("Require positive grading", "SubmoduleImpl::myMinGens()");
    myGBasis(CheckForTimeout);
    return myMinGensValue;
  }


  const std::vector<ModuleElem>& SubmoduleImpl::myTidyGens(const CpuTimeLimit& CheckForTimeout) const
  {
    //    if (!myTidyGensIsValid)
    if (!IsSparsePolyRing(myRing()))
      CoCoA_THROW_ERROR(ERR::NYI, "SubmoduleImpl::myTidyGens");
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
      CoCoA_THROW_ERROR("ERR:NYI coeffs not in a field", "NR");//???
    if ( IsZero(f) ) return f;
    ModuleElem F=f;
    reduce(F, g);
    return F;
  }

  } // anonymous namespace


  bool IsElem(const ModuleElem& v, const module& M)
  {
    CoCoA_ASSERT(IsFGModule(M));
    if (owner(v) != AmbientFreeModule(M))
      CoCoA_THROW_ERROR(ERR::MixedModules, "IsElem(v, M)");
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
    if (gens.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    if (!HasUniqueOwner(gens)) CoCoA_THROW_ERROR(ERR::MixedModules, "submodule(gens)");
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
    if (len(e) != nrows)
      CoCoA_THROW_ERROR(ERR::IncompatDims,"SubmoduleCols");
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
    if (len(e) != ncols)
      CoCoA_THROW_ERROR(ERR::IncompatDims,"SubmoduleRows");
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
    if (g.empty())
      CoCoA_THROW_ERROR("ERR:empty vector", "syz");
    if (!IsField(CoeffRing(RingOf(F))))
      CoCoA_THROW_ERROR("ERR:NYI coeffs not in a field", "syz");//???
    if (NumCompts(F)!=len(g))
      CoCoA_THROW_ERROR("wrong number of components in free module", "syz");
    vector<RingElem> g_non0;
    //-- remove zero gens
    for (long i=0; i<len(g); ++i)
    {
      if ( (!g_non0.empty()) && owner(g[i])!=owner(g_non0[0]))
        CoCoA_THROW_ERROR(ERR::MixedRings, "syz");
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
    if (!IsField(CoeffRing(RingOf(F))))
      CoCoA_THROW_ERROR("ERR:NYI coeffs not in a field", "SyzOfGens");//???
    if (NumCompts(F)!=len(gens(N)))
      CoCoA_THROW_ERROR("wrong number of components in free module", "SyzOfGens");
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
      CoCoA_THROW_ERROR("ERR:NYI coeffs not in a field", "SyzOfGens");//???
    const vector<RingElem>& g = gens(I);
    for (long i=0; i<len(g); ++i)
      if (IsZero(g[i]))
        CoCoA_THROW_ERROR("Zero generator(s): first arg must be the output free module",
                    "SyzOfGens");
    vector<ModuleElem> syzvec;
    FreeModule F = NewFreeModuleForSyz(gens(I), CoCoA_ERROR_CONTEXT);
    return SyzOfGens(F, I);
  }
  

  FGModule SyzOfGens(const FGModule& N)
  {
    if (!IsField(CoeffRing(RingOf(N))))
      CoCoA_THROW_ERROR("ERR:NYI coeffs not in a field", "SyzOfGens");//???
    const vector<ModuleElem>& g = gens(N);
    for (long i=0; i<len(g); ++i)
      if (IsZero(g[i]))
        CoCoA_THROW_ERROR("Zero generator(s): first arg must be the output free module",
                    "SyzOfGens");
    vector<ModuleElem> syzvec;
    FreeModule F = NewFreeModuleForSyz(gens(N), CoCoA_ERROR_CONTEXT);
    return SyzOfGens(F, N);
  }
  

  bool IsContained(const module& M, const module& N)
  {
    if (!IsSparsePolyRing(RingOf(M)))
      CoCoA_THROW_ERROR(ERR::NYI, "SubmoduleImpl::IsContained");
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
////      CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, "IsHomog(submodule)"); 
    if (GradingDim(P)==0)
      CoCoA_THROW_ERROR(ERR::ReqNonZeroGradingDim, "IsHomog(submodule)");
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
////    if (!IsSparsePolyRing(RingOf(M)))
////      CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, "LT(submodule)"); 
////    const SparsePolyRing P = RingOf(M); 
    if (GradingDim(P)==0)
      CoCoA_THROW_ERROR(ERR::ReqNonZeroGradingDim, "LT(submodule)");
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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/submodule.C,v 1.53 2024/07/15 17:12:42 bigatti Exp $
// $Log: submodule.C,v $
// Revision 1.53  2024/07/15 17:12:42  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.52  2024/07/02 09:57:58  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.51  2024/03/26 12:59:32  bigatti
// Summary: improved/rearranged syz  (and using ErrorContext)
//
// Revision 1.50  2024/03/15 18:57:33  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.49  2023/07/04 09:20:13  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.48  2023/06/27 17:45:34  abbott
// Summary: Added check for positive grading to myMinGens
//
// Revision 1.47  2022/08/05 20:48:50  abbott
// Summary: Took NumRows/NumCols out of for-loop exit conds
//
// Revision 1.46  2022/03/07 11:00:46  abbott
// Summary: Added new fn syz(L)
//
// Revision 1.45  2022/03/07 09:47:10  bigatti
// Summary: added func syz(const FreeModule& F, const std::vector<RingElem>& g);
//
// Revision 1.44  2022/02/18 14:12:02  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.43  2022/02/08 20:18:55  abbott
// Summary: Renamed OpenMath output fns (added suffix _OM) (redmine 1528)
//
// Revision 1.42  2021/10/30 16:51:19  abbott
// Summary: Used keyword override (redmine 1625)
//
// Revision 1.41  2020/06/22 15:45:17  abbott
// Summary: ReduceLM now multiplies on the left
//
// Revision 1.40  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.39  2020/02/11 16:56:43  abbott
// Summary: Corrected last update (see redmine 969)
//
// Revision 1.38  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.37  2018/08/02 16:48:43  bigatti
// -- just cleaned up some commented code
//
// Revision 1.36  2018/08/01 15:12:15  bigatti
// -- fixed SyzOfGens in case of generators = 0
//
// Revision 1.35  2018/05/25 09:24:47  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.34  2018/05/18 16:42:11  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.33  2018/05/18 12:25:54  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.32  2018/05/17 15:58:40  bigatti
// -- renamed VectorOperations --> VectorOps
//
// Revision 1.31  2017/12/21 12:19:41  bigatti
// -- added check for coeffs in a field on some operations
//
// Revision 1.30  2017/11/20 20:10:26  bigatti
// -- mimimalized --> SubmoduleOfMinGens
//
// Revision 1.29  2015/06/11 16:57:24  bigatti
// -- using new functions monomial(ring, pp) and monomial(ring, expv)
//
// Revision 1.28  2015/04/13 15:44:46  abbott
// Summary: Changed error code in IsIn and submodule (pseudo-ctor)
// Author: JAA
//
// Revision 1.27  2014/07/31 14:45:19  abbott
// Summary: Merged io.H and UtilsTemplate.H into new header VectorOperations.H
// Author: JAA
//
// Revision 1.26  2014/07/30 14:13:56  abbott
// Summary: Changed BaseRing into RingOf; myBaseRing --> myRing
// Author: JAA
//
// Revision 1.25  2014/07/14 15:10:05  abbott
// Summary: Changed include of tmp.H into UtilsTemplate.H
// Author: JAA
//
// Revision 1.24  2014/07/09 14:27:53  abbott
// Summary: Removed AsFreeModule and AsFGModule
// Author: JAA
//
// Revision 1.23  2014/07/07 13:29:54  abbott
// Summary: Removed AsSparsePolyRing
// Author: JAA
//
// Revision 1.22  2014/04/17 13:39:54  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.21  2014/04/10 13:04:55  bigatti
// -- minimalized(FGModule)
//
// Revision 1.20  2014/04/09 13:11:06  bigatti
// -- added submodule(vector<ModuleElem>)
//
// Revision 1.19  2014/03/26 15:23:44  bigatti
// -- added MinGens for submodules
//
// Revision 1.18  2013/08/02 14:41:41  bigatti
// -- changed LPos --> LPosn
//
// Revision 1.17  2013/07/31 14:47:17  bigatti
// added LT(module)
//
// Revision 1.16  2013/07/31 10:20:33  bigatti
// -- added IsZero, IsHomog
// -- move inliners into class definition
//
// Revision 1.15  2013/06/06 05:47:34  bigatti
// -- pseudoctor is now called "submodule" instead of "NewSubmodule"
// -- added pseudoctor with 1,2,3,4 generators
//
// Revision 1.14  2013/06/03 10:05:59  bigatti
// -- added IsElem, IsContained for modules
//
// Revision 1.13  2013/03/26 14:57:18  abbott
// Corrected silly typo.
//
// Revision 1.12  2013/03/25 17:29:59  abbott
// Changed formal parameter names (M for matrix, N for submodule).
//
// Revision 1.11  2013/03/15 14:58:48  bigatti
// -- added SyzOfGens for modules (something inside is not working though)
//
// Revision 1.10  2013/02/21 16:56:16  bigatti
// -- added GensAsRows, GensAsRows, SyzOfGens
//
// Revision 1.9  2013/01/30 15:46:11  bigatti
// -- added NewSubmoduleCols/Rows
//
// Revision 1.8  2013/01/25 16:15:06  bigatti
// -- added function myGBasis (for SparsePolyRing)
// -- modified myTidyGens
// -- changed myGensArray/myTidyGensArray --> myGensValue/myTidyGensValue
// -- fixed myTidyGensIsValid in constructor
//
// Revision 1.7  2012/10/12 12:38:18  abbott
// Removed element accessor (via operator[]) and non-const mem fn  ModuleBase::myCompt.
//
// Revision 1.6  2012/05/28 09:18:20  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.5  2011/03/10 16:39:33  abbott
// Replaced (very many) size_t by long in function interfaces (for rings,
// PPMonoids and modules).  Also replaced most size_t inside fn defns.
//
// Revision 1.4  2009/12/03 17:26:34  abbott
// Renamed EFGModule to FGModule.
// Renamed ModuleBase member fns  myInit -> myNew, myKill -> myDelete.
// Removed some cruft (old code that was not used by anyone).
//
// Revision 1.3  2008/05/30 12:50:48  abbott
// Aligned some comments.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.5  2007/03/08 18:22:28  cocoa
// Just whitespace cleaning.
//
// Revision 1.4  2007/01/15 14:39:12  bigatti
// -- added prefix "raw" to RawPtr arguments names
//
// Revision 1.3  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 10:15:52  cocoa
// In response to Susan's bug: a fiasco when compiling with CoCoA_MEMPOOL_DEBUG
// set wrongly.  Moved several implementation classes out of their header files
// into the implementation files.  Several functions had to be uninlined.
// Also corrected position of #include, etc.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.7  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.6  2006/04/21 14:56:33  cocoa
// Changed return type of myCompt member function: now it returns a
// ConstRefRingElem instead of a RingElem (i.e. a copy).
//
// Revision 1.5  2006/03/15 18:09:31  cocoa
// Changed names of member functions which print out their object
// into myOutputSelf -- hope this will appease the Intel C++ compiler.
//
// Revision 1.4  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
// Revision 1.3  2005/11/29 13:04:47  cocoa
// -- added "const" to myCompt argument
//
// Revision 1.2  2005/11/24 16:09:38  cocoa
// -- added operator[] for ModuleElem
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.4  2005/04/20 15:40:47  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.3  2005/04/19 14:06:03  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/02/11 14:15:20  cocoa
// New style ring elements and references to ring elements;
// I hope I have finally got it right!
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.6  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.5  2004/11/11 14:06:20  cocoa
// -- moved CVS log to the bottom
// -- minor changes for doxygen
//
// Revision 1.4  2004/06/29 17:10:22  cocoa
// Partially tidied use of "protected" and "private" in various
// base classes.  Checking in at the end of the day -- it works,
// and I wouldn't want it to be lost next time point's disk
// misbehaves.
//
// Revision 1.3  2004/05/24 15:52:13  cocoa
// Major update:
//   new error mechanism
//   many fixes
//   RingHoms almost work now
//   RingFloat much improved
//
// Revision 1.2  2004/01/28 15:37:58  cocoa
// Fairly major update: resuscitated "old style" code which didn't compile
// under the current organization.
//
// Revision 1.1.1.1  2003/09/24 12:55:43  cocoa
// Imported files
//
