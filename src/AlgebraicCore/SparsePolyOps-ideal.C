//   Copyright (c)  2018  John Abbott and Anna M. Bigatti
//   Authors:  2018  John Abbott and Anna M. Bigatti
//             2017  Alice Moallemy (translation CoCoA5-->CoCoALib)

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


// Source code for ideals in SparsePolyRing (functions and member functions)

#include "CoCoA/SparsePolyOps-ideal.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/CpuTimeLimit.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/DenseMatrix.H" // for MultiplicationMat/myDiv
#include "CoCoA/MatrixOps.H" // for LinSolve
#include "CoCoA/MatrixView.H" // for ZeroMat
#include "CoCoA/PPMonoidHom.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H" // for IsQuotientRing
#include "CoCoA/RingFp.H" // for IsRingFp
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H" // for IsQQ
#include "CoCoA/RingZZ.H" // for IsZZ
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-MinPoly.H" // for MinPolyDef
#include "CoCoA/SparsePolyOps-ideal-monomial.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/TmpGOperations.H"  // for myIntersect, my Elim..
#include "CoCoA/TmpPPVector.H"  // for interreducing in GBasisByHomog
#include "CoCoA/TmpUniversalInvolutiveBasisContainer.H" // for ideal ctor
#include "CoCoA/assert.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"  // for myGcd
#include "CoCoA/random.H" // for RandomLongStream
#include "CoCoA/time.H"
#include "CoCoA/verbose.H"
//#include "CoCoA/ideal.H"     // for myGcd
//#include "CoCoA/matrix.H" // for OrdMat, myDivMod

#include <algorithm>
//using std::max;     // for MaxExponent, StdDeg
using std::remove;  // for myColon
using std::sort;    // for AreGoodIndetNames
//#include <functional>
//using std::not1;    // for AreLPPSqFree
//using std::ptr_fun; // for AreLPPSqFree
#include <iostream>
// using std::ostream in SparsePolyRingBase::myOutput
//#include <iterator>
//using std::back_inserter;
#include <list>
//#include <vector>
using std::vector;

namespace CoCoA
{


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGens() const
  { return myGensValue; }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myTidyGens(const CpuTimeLimit& CheckForTimeout) const
  {
    if (!IsField(CoeffRing(myRing())))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);
    return myGBasis(CheckForTimeout);
  }


  const std::vector<RingElem>& GBasis(const ideal& I,
                                      const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
//    if (!IsField(CoeffRing(RingOf(I))))
//      CoCoA_THROW_ERROR(ERR::NotField, "GBasis(I)");
    return TidyGens(I, CheckForTimeout);
  }


  std::vector<RingElem> GBasisTrunc(const ideal& I,
                                    const long TruncDeg,
                                    const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisTrunc(TruncDeg, CheckForTimeout);
  }


  const std::vector<RingElem>& GBasisByHomog(const ideal& I,
                                             const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisByHomog(CheckForTimeout);
  }


  // const std::vector<RingElem>& ReducedGBasisByHomog(const ideal& I,
  //                                                   const CpuTimeLimit& CheckForTimeout)
  // {
  //   const SparsePolyRingBase::IdealImpl* const ptrI = 
  //     SparsePolyRingBase::IdealImpl::ourGetPtr(I);
  //   return ptrI->myGBasisByHomog(CheckForTimeout, );
  // }


  std::vector<RingElem> GBasisSelfSatCore(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisSelfSatCore();
  }


  std::vector<RingElem> GBasisRealSolve(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myGBasisRealSolve();
  }


  const std::vector<RingElem>& ReducedGBasis(const ideal& I)
  {
    if (IsCommutative(RingOf(I))) return GBasis(I); // the same (2017)
    CoCoA_THROW_ERROR2(ERR::NYI, "non-commutative");  // ExtAlg e Weyl?
    return GBasis(I); // just to keep the compiler quiet
  }


  const std::vector<RingElem>& MinGens(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (len(gens(I))>1 && !HasPositiveGrading(RingOf(I)))
      CoCoA_THROW_ERROR1(ERR::ReqPositiveGrading);
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->myMinGens();
  }


  std::vector<ideal> SparsePolyRingBase::IdealImpl::myPrimaryDecomposition() const
  {
    const bool OverField = IsField(CoeffRing(myRing()));
    //if (IhaveMonomialGens()) return myPrimaryDecomposition_MonId();
    if (OverField && IhaveSqFreeMonomialGens()) return myPrimaryDecomposition_MonId();
    if (IamZeroDim()) return myPrimaryDecomposition_0dim();
    CoCoA_THROW_ERROR2(ERR::NYI, "general ideal");
    return vector<ideal>(); // just to keep compiler quiet
  }


  ideal LT(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
  CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const std::vector<RingElem>& GB = TidyGens(I);
    std::vector<RingElem> v;
    const SparsePolyRing P = RingOf(I);
    for (long i=0; i<len(GB); ++i)
      v.push_back(monomial(P, LPP(GB[i])));
    return ideal(P, v);
  }


  ideal LF(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRing P = RingOf(I);
    if ( GradingDim(P)==0 ) CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    const std::vector<RingElem>& GB = TidyGens(I);
    std::vector<RingElem> v;
    for (long i=0; i<len(GB); ++i)  v.push_back(LF(GB[i]));
    return ideal(P, v);
  }


  ideal homog(const ideal& I, ConstRefRingElem x)
  {
    if (!IsSparsePolyRing(RingOf(I)))  CoCoA_THROW_ERROR1(ERR::ReqSparsePolyRing);
    if (AreGensMonomial(I)) return I;
    if (IsZero(I)) return I;
    std::vector<RingElem> HomogIdealGens;
    std::vector<RingElem> v(1,x);
    ComputeHomogenization(HomogIdealGens, gens(I), v);
    return ideal(HomogIdealGens);
  }


  ideal IdealOfGBasis(const ideal& I)
  {
    ideal J(RingOf(I), GBasis(I));
    SetGBasisAsGens(J);
//     if (!uncertain3(IamPrime3Flag)) ...
//     if (!uncertain3(IamMaximal3Flag
    return J;
  }


  ideal IdealOfMinGens(const ideal& I)
  {
    ideal J = I;
    MinGens(J); // so GBasis and such are stored in original ideal
    MakeUnique(J)->myMinimalize();
    return J;
  }


  ideal elim(const ideal& I, const std::vector<RingElem>& ElimIndets) //const CpuTimeLimit& CheckForTimeout)
  {
    if (!IsSparsePolyRing(RingOf(I)))
  CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    //    if (!IsField(CoeffRing(RingOf(I)))) CoCoA_THROW_ERROR(ERR::NotField, "GBasis(I)");
    if (IsZero(I)) return I;
    ideal J(zero(RingOf(I)));
    MakeUnique(J)->myAssignElim(I, ElimIndets);
    return J;
  }


  // Anna: must be friend for ourGetPtr
  bool HasGBasis(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
  CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveGBasis();
  }
  

  // Anna: must be friend for ourGetPtr
  bool AreGensMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
  CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveMonomialGens();
  }
  

  // Anna: must be friend for ourGetPtr
  bool AreGensSqFreeMonomial(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
  CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    return ptrI->IhaveSqFreeMonomialGens();
  }
  

  // Anna: must be friend for ourGetPtr
  void SetGBasisAsGens(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
  CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    const SparsePolyRingBase::IdealImpl* const ptrI = 
      SparsePolyRingBase::IdealImpl::ourGetPtr(I);
    ptrI->mySetGBasisAsGens();
  }


  //-- IdealImpl ----------------------------------------

  ideal SparsePolyRingBase::myIdealCtor(const std::vector<RingElem>& gens) const
  {
    return ideal(new IdealImpl(SparsePolyRing(this), gens)); //??? ugly ???
  }


  SparsePolyRingBase::IdealImpl::IdealImpl(const SparsePolyRing& P, const std::vector<RingElem>& gens):
      myP(P),
//      myGensValue(gens),
      myInvBasisContainerPtr(new Involutive::UniversalInvolutiveBasisContainer(gens)),
      IhaveGBasisValue(false)
  {
    for (const auto& f:gens)
    {
      if (IsZero(f)) continue;
      if (IsOne(f) || (IsField(CoeffRing(P)) && IsConstant(f) ))
      {
        myGensValue.clear();
        myGensValue.push_back(myRing()->myOne());
        break;
      }
      myGensValue.push_back(f);
    }
    // IhaveMonomialGens3Value = uncertain3; // default for bool3
    // IhaveSqFreeMonomial3Gens = uncertain3; // default for bool3
/////    if (!IsField(CoeffRing(P)))
/////      CoCoA_THROW_ERROR("ERR:NYI ideal of polynomials with coeffs not in a field", "ideal(SparsePolyRing, gens)");//???
  }


  IdealBase* SparsePolyRingBase::IdealImpl::myClone() const
  {
    return new IdealImpl(*this);
  }


  const SparsePolyRing& SparsePolyRingBase::IdealImpl::myRing() const
  {
    return myP;
  }


  bool SparsePolyRingBase::IdealImpl::IamZero() const
  {
    for (long i=0; i<len(myGens()); ++i)
      if (!IsZero(myGens()[i])) return false;
    return true;
  }


  
  namespace // anonymous --------------------------------
  {
    bool AreLPPSqFree(const std::vector<RingElem>& v)
    {
      const long n = len(v);
      for (long i=0; i < n; ++i)
        if (!IsSqFree(LPP(v[i]))) return false;
      return true;
//   We *DO NOT USE* STL algorithm because std::ptr_fun does not work if the fn has formal params which are of reference type
//       return find_if(v.begin(), v.end(),
//                      not1(ptr_fun(CoCoA::IsSqFreeLPP)))
// 	//                     not1(ptr_fun(static_cast<bool(*)(ConstRefRingElem)>(CoCoA::IsRadLPP))))
//         == v.end();
    }

    } // anonymous end ----------------------------------


  bool SparsePolyRingBase::IdealImpl::IhaveGBasis() const
  { return IhaveGBasisValue; }


  bool SparsePolyRingBase::IdealImpl::IhaveMonomialGens() const
  {
    if (IsUncertain3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = AreMonomials(myGensValue);
    return IsTrue3(IhaveMonomialGens3Value);
  }


  bool SparsePolyRingBase::IdealImpl::IhaveSqFreeMonomialGens() const
  {
    if (IsUncertain3(IhaveSqFreeMonomialGens3Value))
    {
      if (!IhaveMonomialGens()) IhaveSqFreeMonomialGens3Value = false3;
      else IhaveSqFreeMonomialGens3Value = AreLPPSqFree(myGensValue);
    }
    return IsTrue3(IhaveSqFreeMonomialGens3Value);
  }


  void SparsePolyRingBase::IdealImpl::mySetGBasisAsGens() const
  {
    myGBasisValue = myGensValue;
    IhaveGBasisValue = true;
  }


  void SparsePolyRingBase::IdealImpl::myReset() const
  {
    IamMaximal3Flag = uncertain3;
    IamPrimary3Flag = uncertain3;
    IamPrime3Flag = uncertain3;
    IamRadical3Flag = uncertain3;
    IhaveSqFreeMonomialGens3Value = uncertain3;
    IhaveMonomialGens3Value = uncertain3;
    IhaveGBasisValue = false;
    myGBasisValue.clear();
    myMinGensValue.clear();
  }


  void SparsePolyRingBase::IdealImpl::myTestIsMaximal() const
  {
    if (IamZero()) { myAssignMaximalFlag(false); return; }// <0> in K[X]
    if (!IsField(CoeffRing(myP)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);//???
    if (NumIndets(myP) == 1)  // univariate poly ring 
    { myAssignMaximalFlag(IsIrred(myGBasis(NoCpuTimeLimit())[0])); return; }
    if (IamZeroDim())
      myTestIsMaximal_0dim(); // assigns flags
    else
      myAssignMaximalFlag(false);
  }


  void SparsePolyRingBase::IdealImpl::myTestIsPrimary() const
  {
    //    std::cout << "gens" << myGensValue << std::endl;
    if (!IsField(CoeffRing(myP)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);//???
    if (NumIndets(myP) == 1)  // univariate poly ring 
    {
      const factorization<RingElem> F = factor(myGBasis(NoCpuTimeLimit())[0]); // PID
      myAssignPrimaryFlag(len(F.myFactors()) == 1);
      if (IamPrimary())
        myAssignMaximalFlag(F.myMultiplicities()[0] == 1); // since we have the info
      return;
    }
    if (!IamZeroDim())
      CoCoA_THROW_ERROR("Not Yet Implemented for general ideals", "myTestIsPrimary");
    myTestIsPrimary_0dim();
  }


  void SparsePolyRingBase::IdealImpl::myTestIsPrime() const
  {
    if (NumIndets(myP) == 1 && IsField(CoeffRing(myP)))
    { myTestIsMaximal(); return; }
    CoCoA_THROW_ERROR(ERR::NYI, "SparsePolyRingBase::IdealImpl::myTestIsPrime()");//???
  }


  void SparsePolyRingBase::IdealImpl::myTestIsRadical() const
  {
    if (!IsField(CoeffRing(myP)))  CoCoA_THROW_ERROR1(ERR::ReqCoeffsInField);//???
    if (IhaveMonomialGens()) { myTestIsRadical_MonId(); return; }
    if (NumIndets(myP) == 1)  // univariate poly ring 
    {
      const factorization<RingElem> F = factor(myGBasis(NoCpuTimeLimit())[0]); // PID
      if (len(F.myFactors()) != 1) myAssignPrimeFlag(false);
      if (F.myMultiplicities()[0] == 1) myAssignMaximalFlag(true);
      return;
    }
    if (!IamZeroDim())
      CoCoA_THROW_ERROR("Not Yet Implemented for general ideals", "myTestIsRadical");
    myTestIsRadical_0dim();
  }


  void SparsePolyRingBase::IdealImpl::myReduceMod(RingElemRawPtr rawf) const
  {
    //??? very basic default implementation
    RingElem tmp = NR(RingElemAlias(myP, rawf), myGBasis(NoCpuTimeLimit()));
    myP->mySwap(rawf, raw(tmp));
  }


  bool SparsePolyRingBase::IdealImpl::IhaveElem(RingElemConstRawPtr rawf) const
  {
    RingElem g = RingElemAlias(myP, rawf);
    myReduceMod(raw(g));
    return IsZero(g);
  }


  bool SparsePolyRingBase::IdealImpl::IamZeroDim() const
  {
    const vector<RingElem>& GB = myTidyGens(NoCpuTimeLimit());
    const long GBlen = len(GB);
    const int nvars = NumIndets(myRing());
    vector<bool> AlreadySeen(nvars);
    int NumIndetPowers = 0;
    long index; BigInt IgnoreExp; // for rtn vals from IsIndetPosPower
    for (long i=0; i < GBlen; ++i)
    {
      if (IsIndetPosPower(index, IgnoreExp, LPP(GB[i])) && !AlreadySeen[index])
      {
        AlreadySeen[index] = true;
        if (++NumIndetPowers == nvars) return true;
      }
    }
    return false;
  }
  

  const SparsePolyRingBase::IdealImpl* SparsePolyRingBase::IdealImpl::ourGetPtr(const ideal& I)
  {
    return dynamic_cast<const SparsePolyRingBase::IdealImpl*>(I.myIdealPtr());
  }


  void SparsePolyRingBase::IdealImpl::myAdd(const ideal& Jin)
  {
    if (IsZero(Jin)) return;
    if (!myGensValue.empty() && IsOne(myGensValue[0])) return;
    // ANNA: clever insert (skip 0s?, deal with 1s)
    if ( (IhaveGBasis() && IamOne()) || (HasGBasis(Jin) && IsOne(Jin)) )
    {
      myGensValue.clear();
      myGensValue.push_back(myRing()->myOne());
      myReset();
      return;
    }
    if (IsZero(Jin)) return;
//    myGensValue.insert(myGensValue.end(), gens(Jin).begin(), gens(Jin).end());
    for (const auto& f:gens(Jin))
    {
      if (IsZero(f)) continue;
      if (IsOne(f) || (IsField(CoeffRing(myRing())) && IsConstant(f) ))
      {
        myGensValue.clear();
        myGensValue.push_back(myRing()->myOne());
        break;
      }
      myGensValue.push_back(f);
    }
    bool3 IhaveMG_old = IhaveMonomialGens3Value;
    bool3 IhaveSFMG_old = IhaveSqFreeMonomialGens3Value;
    myReset();
    // we can recover some info about monomial gens:
    IhaveMonomialGens3Value = IhaveMG_old;
    IhaveSqFreeMonomialGens3Value = IhaveSFMG_old;
    const IdealImpl* const J = ourGetPtr(Jin);
    if (IsTrue3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = J->IhaveMonomialGens3Value;
    if (IsTrue3(IhaveSqFreeMonomialGens3Value))
      IhaveSqFreeMonomialGens3Value = J->IhaveSqFreeMonomialGens3Value;
  }


  void SparsePolyRingBase::IdealImpl::myMul(const ideal& Jin)
  {
    const bool OverField = IsField(CoeffRing(RingOf(Jin)));
    if (OverField && IhaveMonomialGens() && AreGensMonomial(Jin))
    {
      // Special handling for monomial ideals over a field
      myMul_MonId(Jin);
      return;
    }
    vector<RingElem> tmpV;
    const SparsePolyRingBase::IdealImpl* const J = ourGetPtr(Jin);
    for (const auto& gI: myGensValue)
      for (const auto& gJ: J->myGensValue)  tmpV.push_back(gI*gJ);
    swap(tmpV, myGensValue);  // ANNA does this make copies???  2010-02-03
    myReset(); // can we recover some info?
  }


  void SparsePolyRingBase::IdealImpl::myIntersect(const ideal& J)
  {
    CoCoA_ASSERT(myRing() == RingOf(J));
    if (IamZero()) return;
    CoCoA_ASSERT(!IsZero(J));
    const bool OverField = IsField(CoeffRing(RingOf(J)));
    if (OverField && IhaveMonomialGens() && AreGensMonomial(J))
    {
      myIntersect_MonId(J);
      return;
    }
    if (!OverField) CoCoA_THROW_ERROR(ERR::NYI, "intersection over non field");
    ComputeIntersection(myGensValue, myGensValue, gens(J));
    myReset(); // can we recover some info?
  }


  void SparsePolyRingBase::IdealImpl::myColon(const ideal& J)
  {
    CoCoA_ASSERT(myRing() == RingOf(J));
    if (IsZero(J))
      myGensValue = vector<RingElem>{one(myRing())};
    else
    {
      const bool OverField = IsField(CoeffRing(RingOf(J)));
      if (OverField && IhaveMonomialGens() && AreGensMonomial(J))
      {
        myColon_MonId(J);
        return;
      }
      const RingElem Z(zero(myRing()));
      myGensValue.erase(remove(myGensValue.begin(), myGensValue.end(),Z),
                        myGensValue.end());
      ComputeColon(myGensValue, myGensValue, gens(J));
    }
    myReset(); // can we recover some info?
  }


  void SparsePolyRingBase::IdealImpl::mySaturate(const ideal& J)
  {
    CoCoA_ASSERT(myRing() == RingOf(J));
    if (IsZero(J))
      myGensValue = vector<RingElem>{one(myRing())};
    else
      ComputeSaturation(myGensValue, myGensValue, gens(J));
    myReset();
  }


  void SparsePolyRingBase::IdealImpl::myMinimalize()
  {
    if (GradingDim(myRing())==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (!IsHomog(myGensValue))  CoCoA_THROW_ERROR1(ERR::ReqHomog);
    myGBasis(NoCpuTimeLimit()); // this sets GBasis and MinGens
    myGensValue = myMinGens(); // if monomial ideal min gens are in GBasis
    if (IsFalse3(IhaveMonomialGens3Value))
      IhaveMonomialGens3Value = uncertain3;
    if (IsFalse3(IhaveSqFreeMonomialGens3Value))
      IhaveSqFreeMonomialGens3Value = uncertain3;
  }


  void SparsePolyRingBase::IdealImpl::myAssignElim(const ideal& I, const std::vector<RingElem>& ElimIndets)
  {
    //    if (HasGBasis(I)) and elim ordering ...
    const bool OverField = IsField(CoeffRing(myRing()));
    if (OverField && AreGensMonomial(I))
    {
      myAssignElim_MonId(I, ElimIndets);
      return;
    }
    const RingElem Z(zero(myRing()));
    PPMonoidElem ElimIndetsProd(myRing()->myPPM());
    for (const auto& x: ElimIndets)
    {
      if (!IsIndet(x))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
      ElimIndetsProd *= LPP(x);
    }
    std::vector<RingElem> ElimGens;
    ComputeElim(ElimGens, gens(I), ElimIndetsProd);
    myReset(); // can we recover some info from I?
    swap(myGensValue, ElimGens);
  }


  namespace {  // anonymous ------------------------------
    RingElem CoeffOfTermSparse(ConstRefRingElem f, ConstRefPPMonoidElem pp)
    {
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        if (PP(itf) == pp) return coeff(itf);
      return zero(CoeffRing(owner(f)));
    }

    matrix MultiplicationMat(ConstRefRingElem f, const ideal& I)
    {
      std::vector<PPMonoidElem> QB = QuotientBasis(I);
      SparsePolyRing P = owner(f);
      matrix Mf = NewDenseMat(CoeffRing(P), len(QB), len(QB));
      for (long j=0; j<len(QB); ++j)
      {
        RingElem tmpf = f;
        P->myMulByPP(raw(tmpf), raw(QB[j]));
        tmpf = NF(tmpf, I);
        for (long i=0; i<len(QB); ++i)
          SetEntry(Mf,i,j, CoeffOfTermSparse(tmpf,QB[i]));
      }
      return Mf;
    }

  } // anonymous end -------------------------------------------


  bool SparsePolyRingBase::IdealImpl::myDivMod(RingElemRawPtr rawlhs, RingElemConstRawPtr rawnum, RingElemConstRawPtr rawden) const
  {
// check for IsZeroDivisor(den) is done by operator/
    const SparsePolyRing P = myRing();
    //    if (IsField(CoeffRing(P)) && P->myIsConstant(rawden))
    if (P->myIsConstant(rawden))
    {
      RingElem GDivF = RingElemAlias(P, rawnum);
      P->myDivByCoeff(raw(GDivF), raw(P->myLC(rawden))); // exc safe?
      P->mySwap(rawlhs, raw(GDivF));
      return true;
    }
    if (IamZeroDim())
    {
      ideal I = ideal(const_cast<IdealImpl*>(this)); // Tappullus Horribilis
      std::vector<PPMonoidElem> QB = QuotientBasis(I);
      long LenQB = len(QB);
      RingElem G = RingElemAlias(P,rawnum);
      RingElem F = RingElemAlias(P,rawden);
      matrix coeffsG = NewDenseMat(ZeroMat(CoeffRing(P), len(QB), 1));
      for (long i=0; i<LenQB; ++i)
        SetEntry(coeffsG,i,0, CoeffOfTermSparse(G, QB[i]));
      matrix coeffsGDivF = LinSolve(MultiplicationMat(F,I), coeffsG);
      if (NumRows(coeffsGDivF) == 0) return false; // no quotient in ring.
      RingElem GDivF(P);
      for (long i=0; i<LenQB; ++i)
        GDivF += monomial(P, coeffsGDivF(i,0), QB[i]);
      P->mySwap(rawlhs, raw(GDivF));
      return true;
    }

// The following should be a general solution...
///    auto tmp = GenRepr(rawnum, ideal(rawden)+this);
///    return tmp[0];

    CoCoA_THROW_ERROR2(ERR::NYI, "general case");
    return false; // just to keep the compiler quiet!!
  }
  

//////   GBasis //////////////////////////////

  void SparsePolyRingBase::IdealImpl::myGBasis_EasyCases() const
  {
    if (IhaveGBasis()) return;
    const bool OverField = IsField(CoeffRing(myRing()));
    if (OverField && IhaveMonomialGens())      myGBasis_MonId();
    else if (IamZero())
    {
      CoCoA_ASSERT(myGBasisValue.empty());
      CoCoA_ASSERT(myMinGensValue.empty());
      IhaveGBasisValue = true;
    }
    else if (len(myGensValue)==1)
    {
      myMinGensValue = myGensValue;
      myGBasisValue.push_back(monic(myGensValue[0]));
      IhaveGBasisValue = true;
    }
//  general case:  IhaveGBasisValue = false;
  }


  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasis(const CpuTimeLimit& CheckForTimeout) const
  {
    myGBasis_EasyCases();
    if (IhaveGBasis()) return myGBasisValue;
    vector<RingElem> MinGens;
    ComputeGBasis(myGBasisValue, MinGens, myGensValue, CheckForTimeout);
    if (!MinGens.empty()) // non-trivial MinGens is non-empty only if ideal is homog
      myMinGensValue = MinGens;
    IhaveGBasisValue = true;
    return myGBasisValue;
  }


  std::vector<RingElem> SparsePolyRingBase::IdealImpl::myGBasisTrunc(long TruncDeg, const CpuTimeLimit& CheckForTimeout) const
  {
    myGBasis_EasyCases();
    if (IhaveGBasis()) return myGBasisValue;
    return myGBasisTrunc_compute(TruncDeg, CheckForTimeout);
  }
  
  std::vector<RingElem> SparsePolyRingBase::IdealImpl::myGBasisTrunc_compute(long TruncDeg, const CpuTimeLimit& CheckForTimeout) const
  {
    if (TruncDeg <= 0) CoCoA_THROW_ERROR1(ERR::ReqPositive); 
    vector<RingElem> MinGens;
    vector<RingElem> GB;
    ComputeGBasisTrunc(GB, MinGens, TruncDeg, myGensValue, CheckForTimeout);
    if (!MinGens.empty()) // MinGens is non-empty only if ideal is homog
      myMinGensValue = MinGens;
    if (TruncDeg != 0) // GBasis is actually truncated
      return GB;
    IhaveGBasisValue = true;
    swap(myGBasisValue, GB);
    return myGBasisValue;
  }

  
  namespace { // anonymous
        
    PPOrdering HomogOrd(const PPOrdering& PPO)
    {
      ConstMatrixView M(OrdMat(PPO));
      long n = NumCols(M);
      bool IsM0NonNeg = true;
      for (long j=0; j<n; ++j) if (M(0,j)<0) IsM0NonNeg = false;   
      RingElem I = one(RingZZ());
      RingElem Z = zero(RingZZ());
      if (IsM0NonNeg)
        return NewMatrixOrdering(MakeTermOrdMat(BlockMat2x2(
                                                            RowMat(OrdMat(PPO),0),
                                                            RowMat(vector<RingElem>(1,I)),
                                                            OrdMat(PPO),
                                                            ColMat(vector<RingElem>(n,Z)))),
                                 1);
      return NewMatrixOrdering(MakeTermOrdMat(ConcatVer(
                                                        RowMat(vector<RingElem>(n+1,I)),
                                                        BlockMat2x2(
                                                                    RowMat(vector<RingElem>(n,I)),
                                                                    RowMat(vector<RingElem>(1,Z)),
                                                                    OrdMat(PPO),
                                                                    ColMat(vector<RingElem>(n,Z))))),
                               1);
    }

    
    SparsePolyRing NewPolyRingForHomog(SparsePolyRing P)
    {
      long n = NumIndets(P);
      if (IsStdDegRevLex(ordering(PPM(P))))
        return NewPolyRing(CoeffRing(P), NewSymbols(n+1));
      if (IsLex(ordering(PPM(P))))
        return NewPolyRing(CoeffRing(P), NewSymbols(n+1), StdDegLex);
      return NewPolyRing(CoeffRing(P), NewSymbols(n+1), HomogOrd(ordering(PPM(P))));
    }

  }
  

  namespace { // anonymous
    bool LPPLessThan(const RingElem& f, const RingElem& g)
    { return LPP(f) < LPP(g); }
  } // anonymous namespace

  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myGBasisByHomog(const CpuTimeLimit& CheckForTimeout) const
  {
    const bool OverField = IsField(CoeffRing(myRing()));
    if (IhaveGBasis()) return myGBasisValue;
    if (OverField && IhaveMonomialGens()) return myGBasis_MonId();
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    // check GradingDim
    // make correct ordering
    const SparsePolyRing&  P = myRing();
    SparsePolyRing Ph = NewPolyRingForHomog(P);
    const long n = NumIndets(P);
    const RingElem& h = indet(Ph, n);
    const std::vector<RingElem>& X = indets(Ph);
    const RingHom phi =PolyAlgebraHom(P,Ph,vector<RingElem>(X.begin(),X.begin()+n));
    std::vector<RingElem> v = indets(P);  v.push_back(one(P));
    const RingHom dehom = PolyAlgebraHom(Ph,P,v);
    const std::vector<RingElem>& g = myGens();
    std::vector<RingElem> gh; // use transform here???
    for (long i=0; i<len(g); ++i)
      if (!IsZero(g[i]))
        gh.push_back(homog(phi(g[i]), h));
    std::vector<RingElem> GB = dehom( GBasis(ideal(gh),CheckForTimeout) );
    gh.clear();
    // interreduce GB
    std::vector<RingElem> GBaux;
    sort(GB.begin(), GB.end(), LPPLessThan);
    //    PPVector LT_GB(PPM(P), NewDivMaskEvenPowers());
    PPWithMask LT_GBi(one(PPM(P)), NewDivMaskEvenPowers());
    PPVector LT_GB_min(PPM(P), NewDivMaskEvenPowers());
    for (long i=0; i<len(GB); ++i)
    {
      //      PushBack(LT_GB, LPP(GB[i]));
      LT_GBi = LPP(GB[i]);
      if (!IsDivisible(LT_GBi, LT_GB_min))
      {
        PushBack(LT_GB_min, LT_GBi);
        GBaux.push_back(monic(NR(GB[i], GBaux)));
      }
    }
    swap(myGBasisValue, GBaux);
    IhaveGBasisValue = true;
    return myGBasisValue;
  }


  std::vector<RingElem> SparsePolyRingBase::IdealImpl::myGBasisSelfSatCore() const
  {
    if (IhaveGBasis()) return myGBasisValue;
    const bool OverField = IsField(CoeffRing(myRing()));
    if (OverField && IhaveMonomialGens()) return myGBasis_MonId();
    CoCoA_ASSERT(myGBasisValue.empty());
    if (IamZero()) return myGBasisValue;
    vector<RingElem> GBSelfSat;
    ComputeGBasisSelfSatCore(GBSelfSat, myGensValue);
    return GBSelfSat;
  }


  std::vector<RingElem> SparsePolyRingBase::IdealImpl::myGBasisRealSolve() const
  {
    if (IamZero()) return myGBasisValue;
    vector<RingElem> GBRealSolve;
    ComputeGBasisRealSolve(GBRealSolve, myGensValue);
    return GBRealSolve;
  }


  namespace 
  {
    long MaxDeg(const std::vector<RingElem>& F)
    {
      long D = 0, d;
      for (const auto& f:F)  if ( (d=ConvertTo<long>(wdeg(f)[0])) > D)  D = d;
      return D;
    }
  }
  
  
  const std::vector<RingElem>& SparsePolyRingBase::IdealImpl::myMinGens() const
  {
    myGBasis_EasyCases();
    if (!myMinGensValue.empty()) return myMinGensValue;
    const bool OverField = IsField(CoeffRing(myRing()));
    if (OverField && IhaveMonomialGens()) return myMinGens_MonId();
    if (GradingDim(myRing())==0 || !IsHomog(myGensValue))
      CoCoA_THROW_ERROR2(ERR::ReqHomog, "try MinSubsetOfGens instead");
    //// WARNING: positive grading checked in fn "MinGens(I)", not here
    myGBasisTrunc_compute(MaxDeg(myGensValue), NoCpuTimeLimit());
    return myMinGensValue;
  }


  bool IsZeroDim(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    if (IsZero(I) || IsOne(I)) return false;
    // Now we know I is non-trivial.
//     const SparsePolyRing P = RingOf(I);
//     const vector<RingElem>& GB = TidyGens(I);
//     const long GBlen = len(GB); // MUST BE A REDUCED GBASIS !!!
//     long NumIndetPowers = 0;
//     for (long i=0; i < GBlen; ++i)
//       if (IsIndetPosPower(LPP(GB[i])))
//         ++NumIndetPowers;
//     return (NumIndetPowers == NumIndets(P));
    return SparsePolyRingBase::IdealImpl::ourGetPtr(I)->IamZeroDim();
  }


  bool IsHomog(const ideal& I)
  {
    VerboseLog VERBOSE("IsHomog(ideal)");
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    if (GradingDim(RingOf(I))==0)  CoCoA_THROW_ERROR1(ERR::ReqNonZeroGradingDim);
    if (IsZero(I)) return true;
    if (AreGensMonomial(I))  // Slug #1739
    {
      VERBOSE(90) << "AreGensMonomial(I) is true" << std::endl;
      return true;
    }
    if (IsHomog(gens(I)))  // Slug #1739
    {
      VERBOSE(90) << "IsHomog(gens(I))) is true" << std::endl;      
      return true;
    }
    // Now we know I is non-trivial.
    return IsHomog(TidyGens(I)); // MUST BE A REDUCED GBASIS !!!
  }


  RingElem DenSigma(const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    if (!IsFractionField(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR2(ERR::NotFracField, "CoeffRing");
    return CommonDenom(ReducedGBasis(I));
  }

  // template???
  bool IsSigmaGoodPrime(const BigInt& p, const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    if (!IsQQ(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR2(ERR::BadRing, "CoeffRing must be RingQQ");
    if (!IsPrime(p))
      CoCoA_THROW_ERROR2(ERR::BadArg, "first argument must be a prime number");
    return !IsDivisible(CommonDenom(ReducedGBasis(I)), p);
  }

  // template???
  bool IsSigmaGoodPrime(const long p, const ideal& I)
  {
    if (!IsSparsePolyRing(RingOf(I)))
      CoCoA_THROW_ERROR2(ERR::ReqSparsePolyRing, "ring of I");
    if (!IsQQ(CoeffRing(RingOf(I))))
      CoCoA_THROW_ERROR2(ERR::BadRing, "Coeff ring must be RingQQ");
    if (!IsPrime(p))
      CoCoA_THROW_ERROR2(ERR::BadArg, "first argument must be a prime number");
    return !IsDivisible(CommonDenom(ReducedGBasis(I)), p);
  }
  

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-ideal.C,v 1.76 2024/09/04 15:06:03 bigatti Exp $
// $Log: SparsePolyOps-ideal.C,v $
// Revision 1.76  2024/09/04 15:06:03  bigatti
// Summary: Error Codes: now using ReqCoeffsInField
//
// Revision 1.75  2024/08/01 16:50:25  bigatti
// Summary: error codes: Req(Elem)SparsePolyRing;
// now using CoCoA_THROW_ERROR1/2
//
// Revision 1.74  2024/07/23 13:40:37  bigatti
// Summary: moved myGBasis_EasyCases from myGBasisTrunc into myGBasisTrunc_compute/myMinGens
//
// Revision 1.73  2024/07/23 13:24:05  bigatti
// Summary: added myGBasisTrunc_compute and myMinGens_MonId (Redmine #1824)
//
// Revision 1.72  2024/07/22 16:21:39  abbott
// Summary: Cleaned myGBasis_EasyCases: if ideal is 0, it sets myGBValue & myMinGensValue
//
// Revision 1.71  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.70  2024/07/02 15:36:15  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.69  2024/07/02 09:57:58  bigatti
// Summary: error codes, first batch:
// ReqUnivariate, ReqNonZero, ReqNonZeroGradingDim, ReqNonZeroModulus,  ReqNonZeroRingElem
//
// Revision 1.68  2024/05/30 09:42:16  bigatti
// Summary: removed MinGens_almost (redmine #1740)
//
// Revision 1.67  2024/05/28 12:56:04  bigatti
// Summary: added GBasisTrunc, myGBasisTrunc (used by myMinGens);
// separated common code myGBasis_EasyCases
//
// Revision 1.66  2024/05/24 14:37:44  bigatti
// Summary: added MinGens_almost
//
// Revision 1.65  2024/03/22 13:48:20  bigatti
// Summary: GBasisByHomog for all orderings - Redmine #1212
//
// Revision 1.64  2024/03/15 18:55:25  abbott
// Summary: Corrected impl of mySaturate (redmine 1790)
//
// Revision 1.63  2024/03/15 15:40:39  bigatti
// Summary: changed myElim and myElim_MonId into myAssignElim and myAssignElim_MonId.
// Added function elim.
// Extended GBasisByHomog to any ordering.
//
// Revision 1.62  2024/02/15 15:30:50  bigatti
// Summary: clever management of 0 and 1 in ideal ctor and ideal sum redmine #1647
//
// Revision 1.61  2024/02/07 10:25:10  bigatti
// Summary: modified  SparsePolyRingBase::IdealImpl::IdealImpl(P, gens)
// detecting 1 and 0.  All other ctors call this one
//
// Revision 1.60  2024/02/06 10:58:14  bigatti
// Summary: just removed commented out code (radical)
//
// Revision 1.59  2024/02/05 11:25:35  bigatti
// Summary: improved IsHomog;  moved radical code to radical
//
// Revision 1.58  2023/07/04 09:20:13  abbott
// Summary: Changed TimeOut to Timeout
//
// Revision 1.57  2023/06/27 17:56:44  abbott
// Summary: Disabled Positive grading check if ideal is created principal
//
// Revision 1.56  2023/06/27 17:45:00  abbott
// Summary: Added check for positive grading to MinGens
//
// Revision 1.55  2023/03/27 13:59:17  abbott
// Summary: Improved err mesg
//
// Revision 1.54  2022/11/26 13:58:54  abbott
// Summary: Fixed redmine 1714
//
// Revision 1.53  2022/02/18 14:11:59  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.52  2022/02/04 21:41:06  abbott
// Summary: Changed name MakeTermOrd to MakeTermOrdMat (redine 854)
//
// Revision 1.51  2022/01/21 08:56:44  abbott
// Summary: Corrected previous check-in comment
//
// Revision 1.50  2022/01/20 19:16:14  abbott
// Summary: Added check for zero poly in myGBasisByHomog (redmine 1646)
//
// Revision 1.49  2021/10/04 08:58:57  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.48  2021/04/28 12:43:20  abbott
// Summary: Short-cuts to monomial ideals only if coeffs are in a field (redmine 1381)
//
// Revision 1.47  2021/04/28 09:04:52  abbott
// Summary: Removed check that coeffring is field in pseudo-ctor for ideal; added checks that coeffring is field in myMul "monomial short-cut", and in myTidyGens (redmine 1381)
//
// Revision 1.46  2021/03/03 22:10:02  abbott
// Summary: Changed arg names
//
// Revision 1.45  2020/06/17 15:49:28  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.44  2020/03/06 18:40:07  abbott
// Summary: Added some const
//
// Revision 1.43  2020/03/06 16:38:48  bigatti
// -- added GBasisByHomog with DegLex
//
// Revision 1.42  2020/02/26 15:29:15  bigatti
// -- commented out ASSERT for GBasis (called by MinGens redmine#1420)
//
// Revision 1.41  2020/02/14 09:01:10  bigatti
// -- removed typo
//
// Revision 1.40  2020/02/14 08:58:52  bigatti
// -- mingens: for gbasis computation with Buchberger,
//    if myMinGensValue not set (IdealOfProjectivePoints)
//
// Revision 1.39  2020/02/13 14:01:23  bigatti
// -- now using auto in for loops
//
// Revision 1.38  2020/02/12 10:27:38  bigatti
// -- improved myAssignBLAHFlag
//
// Revision 1.37  2020/02/12 09:42:17  abbott
// Summary: Corrected defn of myTestIsPrimary
//
// Revision 1.36  2020/02/12 09:01:47  bigatti
// -- changed myTestIsMaximal etc to return void (and consequences)
//
// Revision 1.35  2020/01/26 14:37:24  abbott
// Summary: Removed useless include (of NumTheory)
//
// Revision 1.34  2020/01/13 17:42:38  bigatti
// -- minor change in GBasisByHomog
//
// Revision 1.33  2019/12/28 17:54:50  abbott
// Summary: New defn of IamZeroDim now works with any GB (previous reqd minGB)
//
// Revision 1.32  2019/12/27 17:42:31  bigatti
// -- fixed interreduction in GBasisByHomog
// -- IsZeroDim --> IamZeroDim in member functions
//
// Revision 1.31  2019/10/21 16:40:46  bigatti
// -- added VERBOSE(1000) for some functions (to see if they are called)
//
// Revision 1.30  2019/10/15 12:57:55  bigatti
// -- renamed files for ideals
//
// Revision 1.29  2019/10/03 14:54:37  bigatti
// -- added radical_MonId
//
// Revision 1.28  2019/10/02 10:40:47  bigatti
// -- minor change in view of adding optimized radical for MonomialIdeals
//
// Revision 1.27  2019/09/27 14:56:55  abbott
// Summary: Fixed redmine 1322; some cleaning
//
// Revision 1.26  2019/09/27 14:44:06  abbott
// Summary: Corrected use of empty inside CoCoA_ASSERT
//
// Revision 1.25  2019/09/25 16:04:19  bigatti
// -- gbasis for principal ideals now just computes monic(f)
// -- added radical_tmp, checking the implemented cases (will become radical)
//
// Revision 1.24  2019/03/04 11:14:04  abbott
// Summary: Added "just to keep compiler quiet"
//
// Revision 1.23  2018/12/17 15:21:41  bigatti
// -- timeout in GBasisByHomog
//
// Revision 1.22  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.21  2018/08/06 09:38:29  bigatti
// -- renamed GBasisViaHomog --> GBasisByHomog
//
// Revision 1.20  2018/08/06 08:57:48  bigatti
// -- added timeout for GBasisByHomog
// -- now using GBasisByHomog in IsPrimary and radical
//
// Revision 1.19  2018/08/05 16:31:25  bigatti
// -- added GBasisByHomog
//
// Revision 1.18  2018/06/27 12:15:19  abbott
// Summary: Renamed RealSolveCore to RealSolve
//
// Revision 1.17  2018/06/27 08:50:39  abbott
// Summary: Revised to work with new CpuTimeLimit
//
// Revision 1.16  2018/05/25 09:24:46  abbott
// Summary: Major redesign of CpuTimeLimit (many consequences)
//
// Revision 1.15  2018/05/18 12:23:50  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.14  2018/05/17 15:49:00  bigatti
// -- sorted #includes
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.13  2018/04/18 15:39:34  abbott
// Summary: Added missing return in NYI fn.
//
// Revision 1.12  2018/04/16 21:49:26  bigatti
// -- added IamZeroDim
// -- added myPrimaryDecomposition_0dim
//
// Revision 1.11  2018/04/10 14:59:32  bigatti
// -- now member function for PrimaryDecomposition
//
// Revision 1.10  2018/04/10 14:51:43  bigatti
// -- added virtual myPrimaryDecomposition (with default implementation)
//
// Revision 1.9  2018/04/10 14:20:47  bigatti
// -- fixed includes
//
// Revision 1.8  2018/04/09 16:25:53  bigatti
// -- functions for zer-dim ideals moved into SparsePolyOps-IdealZeroDim
//
// Revision 1.7  2018/04/06 17:28:22  bigatti
// -- fixed includes
//
// Revision 1.6  2018/04/04 12:36:51  bigatti
// -- radical: returning I itself, in IsRadical3 = true
//
// Revision 1.5  2018/03/29 13:09:04  bigatti
// -- improved isprimary (with gbasis timeout)
//
// Revision 1.4  2018/03/29 09:36:40  bigatti
// -- added member functions myTestIsRadical, myTestIsPrimary and flags
//
// Revision 1.3  2018/03/20 11:48:12  bigatti
// -- new code for 0dim ideals (radical, primary...)
//
// Revision 1.2  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.1  2018/03/12 14:38:41  bigatti
// -- renamed SparsePoly-ideal/involutive into SparsePolyOps-ideal/involutive
//
// Revision 1.3  2018/03/09 14:54:01  bigatti
// -- removed useless includes
//
// Revision 1.2  2018/02/27 17:26:07  bigatti
// -- split off involutive part
//
// Revision 1.1  2018/02/27 17:12:37  bigatti
// -- Renamed SparsePolyRing_ideal.C --> SparsePoly-ideal.C
//
