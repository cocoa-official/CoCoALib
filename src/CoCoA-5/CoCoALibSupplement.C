//   Copyright (c) 2010-2017,2022  John Abbott,  Anna M Bigatti
//   Orig authors: 2010-2017 Anna Maria Bigatti, 2010 Giovanni Lagorio (lagorio@disi.unige.it)
//
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

#include "CoCoALibSupplement.H"
#include "CoCoA/library.H"

#include <sstream>  // for ErrorMessage
#include <string>
//#include <vector> // already included by library.H
using std::vector;

namespace CoCoA {


  namespace // anonymous namespace
  {
    // for NBM_forC5 and QuotientBasis_forC5
    std::vector<RingElem> VectorRingElem(const SparsePolyRing& P, const std::vector<PPMonoidElem>& v)
    {
      const long N = len(v);
      vector<RingElem> ans; ans.reserve(N);
      for (long i=0; i < N; ++i)
        ans.push_back(monomial(P, 1, v[i]));
      return ans;
    }
  } // anonymous namespace


  std::vector<long> VectorLong(const std::vector<BigInt>& BigIntVec, const std::string& FuncName)
  {
    const ErrorInfo ErrMesg(ERR::BadSymbolRange, FuncName);
    vector<long> v;
    for (long i=0; i<len(BigIntVec); ++i)
      v.push_back(ConvertTo<long>(BigIntVec[i],  ErrMesg));
    return v;
  }


  std::vector<long> VectorLongDecr1(const std::vector<BigInt>& BigIntVec, const ERR::ID& ErrID, const std::string& FnName)
  {
    const ErrorInfo ErrMesg(ErrID, FnName);
    vector<long> v;
    for (long i=0; i<len(BigIntVec); ++i)
      v.push_back(-1+ConvertTo<long>(BigIntVec[i],  ErrMesg));
    return v;
  }


  std::string fold_forC5(const std::string& str, const BigInt& MaxLineLen)
  {
    const ErrorInfo ErrMesg(ERR::OutOfRange, "fold");
    const long MaxWidth = ConvertTo<long>(MaxLineLen, ErrMesg);
    if (MaxWidth < 1)  throw ErrMesg; /// ??? CoCoA_THROW_ERROR1(ERR::OutOfRange); ????
    return fold(str, MaxWidth);
  }


  bool IsTerm_forC5(ConstRefRingElem f)
  {
    return IsMonomial(f)&&IsOne(LC(f));
  }


  long UnivariateIndetIndex_forC5(ConstRefRingElem f)
  {
    return 1 + UnivariateIndetIndex(f);
  }


  RingElem LT_forC5(ConstRefRingElem f)
  {
    return LPP_forC5(f);
  }  ///< NB result belongs to ring (owner(f))


  ModuleElem LT_forC5(const ModuleElem& f)
  {
    //    return LPP(f)*gens(owner(f))[LPos(f)];
    return LPP_forC5(f)*gens(owner(f))[LPosn(f)];
  }  ///< NB result belongs to module


  long LPosn_forC5(const ModuleElem& f)
  { return LPosn(f)+1; }


  long FirstNonZeroPosn_forC5(const ModuleElem& f)
  { return FirstNonZeroPosn(f)+1; }


  RingElem LPP_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsSparsePolyRing(R)) return monomial(R, LPP(f));
    return IndetPower(R, 0, deg(f)); // univariate case
  }


  RingElem LPP_forC5(const ModuleElem& f)
  {
    const ring R = RingOf(owner(f));
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsSparsePolyRing(R)) return monomial(R, LPP(f));
    //    PolyRing Rx = AsPolyRing(R); // univariate case
    //    return IndetPower(Rx, 0, deg(f));
    CoCoA_THROW_ERROR2(ERR::NYI, "LPP_forC5 univariate module");
    return IndetPower(R, 0, 1);
  }


  RingElem LM_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsSparsePolyRing(R)) return monomial(R, LC(f), LPP(f));
    return LC(f)*IndetPower(R, 0, deg(f));
  }


  ModuleElem LM_forC5(const ModuleElem& f)
  {
    long LPosnf = LPosn(f);
    return LM_forC5(f[LPosnf])*gens(owner(f))[LPosnf];
  }  ///< NB result belongs to module


  // bool IsCoprime_forC5(ConstRefRingElem f, ConstRefRingElem g)
  // {
  //   if (owner(f) != owner(g))
  //     CoCoA_THROW_ERROR1(ERR::MixedRings);
  //   if (IsZZ(owner(f))) return IsCoprime(f,g); // integers
  //   if (!IsPolyRing(owner(f)))
  //     CoCoA_THROW_ERROR2(ERR::BadArg, "Args must be PP, or elem of ZZ or in poly ring over field");
  //   if (IsMonomial(f) && IsOne(LC(f)) && IsMonomial(g) && IsOne(LC(g)))
  //     return IsCoprime(LPP(f), LPP(g)); // power products; NB result belongs to ring (owner(f))

  //   if (!IsField(CoeffRing(owner(f))))
  //     CoCoA_THROW_ERROR1("Args must be PP, or elem of ZZ or in poly ring over field");
  //   return IsCoprime(f,g); // polys in k[x,y,...] with k a field
  // }


  namespace
  {

    RingElem CoeffOfTermSparse(ConstRefRingElem f, ConstRefPPMonoidElem pp)
    {
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
//??? [redmine 1545]     {  const int sign = cmp(PP(itf), pp); if (sign < 0) break; if (sgn == 0) return coeff(itf); }
        if (PP(itf) == pp) return coeff(itf);
      return zero(CoeffRing(owner(f)));
    }

    std::vector<RingElem> CoefficientsDense(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const DenseUPolyRing P = owner(f);
      for (long i=P->myDegPlus1(raw(f))-1; i>=0; --i)
        if (!IsZero(coeff(f,i))) v.push_back(coeff(f,i));
      return v;
    }

    std::vector<RingElem> CoefficientsSparse(ConstRefRingElem f)
    {
      vector<RingElem> v;
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        v.push_back(coeff(itf));
      return v;
    }

    std::vector<RingElem> MonomialsDense(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const DenseUPolyRing P = owner(f);
      for (long i=P->myDegPlus1(raw(f))-1; i>=0; --i)
        if (!IsZero(coeff(f,i))) v.push_back(monomial(P, coeff(f,i), i));
      return v;
    }

    std::vector<RingElem> MonomialsSparse(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const SparsePolyRing P = owner(f);
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        v.push_back(monomial(P, coeff(itf), PP(itf)));
      return v;
    }

    std::vector<RingElem> SupportDense(ConstRefRingElem f)
    {
      CoCoA_THROW_ERROR2(ERR::NYI, "SupportDense is disabled");
      vector<RingElem> v;
      const DenseUPolyRing P = owner(f);
      for (long i=P->myDegPlus1(raw(f))-1; i>=0; --i)
        if (!IsZero(coeff(f,i))) v.push_back(monomial(P, one(CoeffRing(P)), i));
      return v;
    }

    std::vector<RingElem> SupportSparse(ConstRefRingElem f)
    {
      vector<RingElem> v;
      const SparsePolyRing P = owner(f);
      for (SparsePolyIter itf=BeginIter(f); !IsEnded(itf); ++itf)
        v.push_back(monomial(P, PP(itf)));
      return v;
    }
  }


  RingElem CoeffOfTerm_forC5(ConstRefRingElem f, ConstRefRingElem t)
  {
    const ring& R = owner(f);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (owner(t) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsMonomial(t) || !IsOne(LC(t)))  CoCoA_THROW_ERROR2(ERR::BadArg, "arg 2 should be a term");

    if (IsSparsePolyRing(R))
      return CoeffOfTermSparse(f, LPP(t));
    return coeff(f, deg(t));
  }


  std::vector<RingElem> coefficients_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsSparsePolyRing(R))  return CoefficientsSparse(f);
    return CoefficientsDense(f);
  }


  std::vector<RingElem> monomials_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsSparsePolyRing(R))  return MonomialsSparse(f);
    return MonomialsDense(f);
  }


  std::vector<RingElem> support_forC5(ConstRefRingElem f)
  {
    ring R = owner(f);
    if (!IsPolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    if (IsSparsePolyRing(R))  return SupportSparse(f);
    return SupportDense(f);
  }


  bool IsFactorClosed_forC5(const std::vector<RingElem>& v)
  {
    if (v.empty()) CoCoA_THROW_ERROR1(ERR::ReqNonEmpty);
    const long n = len(v);
    if (n == 1 && IsOne(v[0])) return true;
    const ring R = owner(v[0]);
    for (long i=1; i < n; ++i)
      if (owner(v[i]) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (!IsSparsePolyRing(R))  CoCoA_THROW_ERROR1(ERR::ReqElemPolyRing);
    for (long i=0; i < n; ++i)
      if (!IsMonomial(v[i]))  CoCoA_THROW_ERROR1(ERR::ReqMonomialGens);
    vector<PPMonoidElem> SetPP; SetPP.reserve(n);
    for (long i=0; i < n; ++i)
      SetPP.push_back(LPP(v[i]));
    return IsFactorClosed(SetPP);
  }


  namespace
  {
    // needs a decent implementation, this one isn't
    RingElem DensePolyRec(const PolyRing& P, long n, long d)
    {
      if (n==1) return IndetPower(P,0,d);
      if (d==1)
      {
        RingElem s=indet(P,0);
        for (long i=1; i<n; ++i)  s+=indet(P,i);
        return s;
      }
      return indet(P,n-1)*DensePolyRec(P,n,d-1) + DensePolyRec(P,n-1,d);
    }
  }


  RingElem DensePoly_forC5(const ring& R, const BigInt& D)
  {
    CoCoA_ASSERT(IsPolyRing(R));
    const long d = ConvertTo<long>(D, ErrorInfo(ERR::ArgTooBig, "DensePoly_forC5"));
    if (d<0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    if (d==0) return one(R);
    return DensePolyRec(R, NumIndets(R), d);
  }


  std::vector<BigInt> DegreeToVec(degree d)
  {
    vector<BigInt> v;
    for (long i=0; i<GradingDim(d); ++i)  v.push_back(d[i]);
    return v;
  }


  std::vector<BigInt> wdeg_forC5(const ModuleElem& f)
  {
    degree d = wdeg(f);
    vector<BigInt> v;
    for (long i=0; i<GradingDim(d); ++i)  v.push_back(d[i]);
    return v;
  }


  BigInt NextPrime_forC5(const BigInt& N)
  {
    if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "NextPrime_forC5");
    const long ans = NextPrime(ConvertTo<long>(N, ErrMesg));
    if (ans == 0)  CoCoA_THROW_ERROR2(ERR::ArgTooBig, "see manual for NextProbPrime");
    return BigInt(ans);
  }


  BigInt PrevPrime_forC5(const BigInt& N)
  {
    if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "PrevPrime_forC5");
    const long ans = PrevPrime(ConvertTo<long>(N, ErrMesg));
    return BigInt(ans);
  }


  factorization<BigInt> FactorINT_forC5(const BigInt& N)
  {
    // 2024-03-18 CoCoALib does arg checking anyway
//    if (N == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    return factor(N);
  }


//  2023-12-22  Commented out since FactorINT automatically times out.
//   factorization<BigInt> FactorINT_timeout_forC5(const BigInt& N, const BigInt& TimeLimit)
//   {
//     if (N == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
// //???    const ErrorInfo ErrMesg(ERR::ArgTooBig, "FactorINT_TimeOut_forC5");
//     long secs; if (!IsConvertible(secs, TimeLimit))  CoCoA_THROW_ERROR2(ERR::ArgTooBig, "FactorINT_timeout");
//     if (secs < 0 || secs > 10000)  CoCoA_THROW_ERROR2("Time limit must be between 0 and 10000", "FactorINT_timeout");
//     return factor(N, CpuTimeLimit(secs));
//   }


  factorization<BigInt> FactorINT_TrialDiv_forC5(const BigInt& N, const BigInt& TrialLimit)
  {
    if (N == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "FactorINT_TrialDiv_forC5");
    return factor_TrialDiv(N, ConvertTo<long>(TrialLimit, ErrMesg));
  }


  factorization<BigInt> FactorINT_PollardRho_forC5(const BigInt& N, const BigInt& TrialLimit)
  {
    if (N == 0)  CoCoA_THROW_ERROR1(ERR::ReqNonZero);
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "FactorINT_PollardRho_forC5");
    return factor_PollardRho(N, ConvertTo<long>(TrialLimit, ErrMesg));
  }


  std::vector<long> AllFactors_forC5(const BigInt& N)
  {
    if (N < 0)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "AllFactors_forC5");
    return AllFactors(ConvertTo<long>(N, ErrMesg));
  }


  BigRat CpuTime_forC5()
  { return ConvertTo<BigRat>(CpuTime()); }


  BigRat ElapsedTime_forC5()
  { return ConvertTo<BigRat>(ElapsedTime()); }


  //  std::string date_forC5()
  BigInt date_forC5()
  {
    long date, time;
    DateTime(date, time);
    return BigInt(date);  //*1000000 + time;
  }


  BigInt TimeOfDay_forC5()
  {
    long date, time;
    DateTime(date, time);
    return BigInt(time);
  }


  void reseed_forC5(const BigInt& seed)
  { reseed(GlobalRandomSource(), seed); }


  const BigInt random_forC5(const BigInt& lo, const BigInt& hi)
  { return RandomBigInt(GlobalRandomSource(), lo, hi); }


  BigInt CertifyNotPower_forC5(BigInt N, BigInt UPBexp)
  {  return BigInt(CertifyNotPower(N, ConvertTo<long>(UPBexp)));  }


  std::vector<long> InvTotient_forC5(const BigInt& N) // N must fit into long!!
  {  return InvTotient(ConvertTo<long>(N));  }

  BigInt lcm_forC5(const std::vector<BigInt>& v)
  {  // v is not empty
    BigInt res = v[0];
    for (long i=1; i<len(v); ++i)  res = lcm(res, v[i]);
    return res;
  }


  BigInt gcd_forC5(const std::vector<BigInt>& v)
  {  // v is not empty
    BigInt res = v[0];
    for (long i=1; i<len(v); ++i)  res = gcd(res, v[i]);
    return res;
  }


  RingElem lcm_forC5(const std::vector<RingElem>& v)
  {  // v is not empty
    RingElem res = v[0];
    for (long i=1; i<len(v); ++i)  res = lcm(res, v[i]);
    return res;
  }


  std::vector<BigInt> CoprimeFactorBasis_forC5(const std::vector<BigInt>& v)
  {
    CoprimeFactorBasis_BigInt gfb;
    gfb.myAddInfo(v);
    return FactorBase(gfb);
  }


  std::vector<RingElem> CoprimeFactorBasis_forC5(const std::vector<RingElem>& v)
  {
    CoprimeFactorBasis_RingElem gfb;
    gfb.myAddInfo(v);
    return FactorBase(gfb);
  }


  RingElem gcd_forC5(const std::vector<RingElem>& v)
  {  // v is not empty
    RingElem res = v[0];
    for (long i=1; i<len(v); ++i)  res = gcd(res, v[i]);
    return res;
  }


  RingElem ContentWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    vector<long> indices(len(v));
    for (long i=0; i<len(v); ++i)
      if (!IsIndet(indices[i], v[i]))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
    return ContentWRT(f, indices);
  }


  std::vector<CoeffPP> CoefficientsWRT_forC5(ConstRefRingElem f, const std::vector<RingElem>& v)
  {
    vector<long> indices(len(v));
    for (long i=0; i<len(v); ++i)
      if (!IsIndet(indices[i], v[i]))  CoCoA_THROW_ERROR1(ERR::ReqIndet);
    return CoefficientsWRT(f, indices);
  }


  std::vector<RingElem> CoeffListWRT_forC5(ConstRefRingElem f, ConstRefRingElem x)
  { return CoeffVecWRT(f, x); }


  std::vector<RingElem> CoeffListWRTSupport_forC5(ConstRefRingElem f, ConstRefRingElem basis)
  { return CoeffVecWRTSupport(f, basis); }



  RingElem HomogCompt_forC5(ConstRefRingElem f, const BigInt& D)
  {
    long d;
    if (!IsConvertible(d, D)) CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return HomogCompt(f,d);
  }


///???  RingElem resultant_forC5(ConstRefRingElem f, ConstRefRingElem g)

  RingElem resultant_forC5(ConstRefRingElem f, ConstRefRingElem g, ConstRefRingElem x)
  {
    long x_index; // will be assigned by call to IsIndet
    if (!IsIndet(x_index, x))  CoCoA_THROW_ERROR2(ERR::ReqIndet, "arg 3");
    return resultant(f,g, x_index);
  }


///???  RingElem discriminant_forC5(ConstRefRingElem f);
  RingElem discriminant_forC5(ConstRefRingElem f, ConstRefRingElem x)
  {
    long x_index; // will be assigned by call to IsIndet
    if (!IsIndet(x_index, x))  CoCoA_THROW_ERROR2(ERR::ReqIndet, "arg 2");
    return discriminant(f, x_index);
  }


  std::vector<RingElem> QuotientBasis_forC5(const ideal& I)
  { return VectorRingElem(RingOf(I), QuotientBasis(I)); }


  std::vector<RingElem> QuotientBasisSorted_forC5(const ideal& I)
  { return VectorRingElem(RingOf(I), QuotientBasisSorted(I)); }


  std::vector<RingElem> JanetBasis_forC5(const ideal& I)
  { return Involutive::JanetBasis(I); }


  void ApproxPointsNBM_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV,
                 const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance)
  {
    std::vector<vector<RingElem> > pts(GetRows(CanonicalHom(RingQQ(), CoeffRing(P))( OrigPts )));
    if (NumRows(OrigTolerance) != 1)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    vector<RingElem> tol(GetRow(OrigTolerance, 0));
    vector<PPMonoidElem> tmpQB0;
    vector<RingElem> BB0;
    vector<RingElem> AV0;
    ApproxPts::NBM(tmpQB0, BB0, AV0, P, pts, tol);
    vector<RingElem> QB0(VectorRingElem(P, tmpQB0));
    swap(QB, QB0);
    swap(BB, BB0);
    swap(AV, AV0);
  }


  void ApproxPointsSOI_forC5(std::vector<RingElem>& QB, std::vector<RingElem>& BB, std::vector<RingElem>& AV,
                             const SparsePolyRing& P, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance)
  {
    std::vector<vector<RingElem> > pts(GetRows(CanonicalHom(RingQQ(), CoeffRing(P))( OrigPts )));
    if (NumRows(OrigTolerance) != 1)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    vector<RingElem> tol(GetRow(OrigTolerance, 0));
    vector<PPMonoidElem> tmpQB0;
    vector<RingElem> BB0;
    vector<RingElem> AV0;
    RingElem gamma(RingQQ(), BigRat(1,10));
    ApproxPts::SOI(tmpQB0, BB0, AV0, P, pts, tol, gamma);
    vector<RingElem> QB0(VectorRingElem(P, tmpQB0));
    swap(QB, QB0);
    swap(BB, BB0);
    swap(AV, AV0);
  }


  RingElem ClosePassingPoly_forC5(const ring& R, ConstMatrixView OrigPts, ConstMatrixView OrigTolerance)
  {
    std::vector<vector<RingElem> > pts(GetRows(CanonicalHom(RingQQ(), CoeffRing(R))( OrigPts )));
    if (NumRows(OrigTolerance) != 1)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    vector<RingElem> tol(GetRow(OrigTolerance, 0));
    vector<RingElem> poly0;
    vector<ApproxPts::PointR> NewPointsVec0;
    RingElem MaxTol=zero(RingQQ());
    for (long i=0; i<len(tol); ++i)
      if (MaxTol<tol[i]) MaxTol = tol[i];
    ApproxPts::VanishPoly(NewPointsVec0, poly0, R, pts, tol, MaxTol);
    return poly0[0];
  }


  void PreprocessPts_forC5(const std::string& WhichAlgm,
                           std::vector< std::vector<RingElem> >& NewPts,
                           std::vector<long>& weights,
                           ConstMatrixView OrigPts,
                           ConstMatrixView epsilon)
  {
    if (!IsOrderedDomain(RingOf(OrigPts)))  CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    if (RingOf(epsilon) != RingOf(OrigPts))  CoCoA_THROW_ERROR1(ERR::MixedRings);
    const long NumPts = NumRows(OrigPts);
    const long dim = NumCols(OrigPts);
    if (NumCols(epsilon) != dim || NumRows(epsilon) != 1)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    for (long j=0; j < dim; ++j)
      if (epsilon(0,j) <= 0) CoCoA_THROW_ERROR2(ERR::ReqPositive, "epsilon");
    vector< vector<RingElem> > pts(NumPts, vector<RingElem>(dim));
    for (long i=0; i < NumPts; ++i)
      for (long j=0; j < dim; ++j)
        pts[i][j] = OrigPts(i,j);
    vector<RingElem>  eps(dim);
    for (long j=0; j < dim; ++j)
      eps[j] = epsilon(0,j);

    if (WhichAlgm == "grid")
      PreprocessPtsGrid(NewPts, weights, pts, eps);
    else if (WhichAlgm == "aggr")
      PreprocessPtsAggr(NewPts, weights, pts, eps);
    else if (WhichAlgm == "subdiv")
      PreprocessPtsSubdiv(NewPts, weights, pts, eps);
    else if (WhichAlgm == "auto")
      PreprocessPts(NewPts, weights, pts, eps);
    else CoCoA_THROW_ERROR2(ERR::ShouldNeverGetHere, "Unknown preprocessing algm");
  }


  RingElem HilbertNumQuot_forC5(const ideal& I)
  {
    PPOrdering PPO(ordering(PPM(RingOf(I))));
    if (GradingDim(PPO) ==0)  // aff -- see hp.cpkg5
      CoCoA_THROW_ERROR2(ERR::NYI, "GradingDim(PPO) is 0");
    //    if (GradingDim(PPO) >1)  // multigraded
    if (GradingDim(PPO) ==1)  // standard case
    {
      if (IsStdGraded(PPO))
        return HilbertNumQuot(I);
      else
        return MGHilbertNumQuot(I);
    }
    return MGHilbertNumQuot(I);
  }


  std::vector<BigInt> ContFrac_forC5(const BigRat& q)
  {
    vector<BigInt> ans;
    if (IsZero(q)) { ans.push_back(BigInt(0)); return ans; }
    for (ContFracIter iter(q); !IsEnded(iter); ++iter)
      ans.push_back(*iter);
    return ans;
  }


  std::vector<BigRat> CFApproximants_forC5(const BigRat& q)
  {
    vector<BigRat> ans;
    if (IsZero(q)) { ans.push_back(BigRat(0)); return ans; }
    for (CFApproximantsIter iter(q); !IsEnded(iter); ++iter)
      ans.push_back(*iter);
    return ans;
  }


  std::vector<BigInt> BinomialRepr_forC5(const BigInt& N, const BigInt& r)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "BinomialRepr");
    std::vector<BigInt> ans = BinomialRepr(N, ConvertTo<long>(r,  ErrMesg));
    ans.erase(ans.begin()); // erase the first elem (which is always 0)
    return ans;
  }

  BigInt BinomialReprShift_forC5(const BigInt& N, const BigInt& r, const BigInt& shift1, const BigInt& shift2)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "BinomialReprShift");
    return BinomialReprShift(N,
                             ConvertTo<long>(r,      ErrMesg),
                             ConvertTo<long>(shift1, ErrMesg),
                             ConvertTo<long>(shift2, ErrMesg));
  }


  BigInt NumPartitions_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "NumPartitions");
    return NumPartitions(ConvertTo<long>(N,  ErrMesg));
  }


  BigInt CatalanNumber_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "CatalanNumber");
    return CatalanNumber(ConvertTo<long>(N,  ErrMesg));
  }


  ring RingQQt_forC5(const BigInt& NumIndets)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, "RingQQt");
    return RingQQt(ConvertTo<long>(NumIndets, ErrMesg));
  }


  //----- matrix --------------------------------------------------

  matrix NewDenseMat_forC5(const ring& R, const BigInt& NR, const BigInt& NC)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return NewDenseMat(R,
                       ConvertTo<long>(NR,  ErrMesg),
                       ConvertTo<long>(NC,  ErrMesg));
  }


  matrix RevLexMat_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return NewDenseMat(RevLexMat(ConvertTo<long>(N,  ErrMesg)));
  }


  matrix StdDegLexMat_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return NewDenseMat(StdDegLexMat(ConvertTo<long>(N,  ErrMesg)));
  }


  matrix StdDegRevLexMat_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return NewDenseMat(StdDegRevLexMat(ConvertTo<long>(N,  ErrMesg)));
  }


  matrix XelMat_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return NewDenseMat(XelMat(ConvertTo<long>(N,  ErrMesg)));
  }


  matrix LexMat_forC5(const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return NewDenseMat(LexMat(ConvertTo<long>(N,  ErrMesg)));
  }


  void SetEntry_forC5(MatrixView& M, const BigInt& I, const BigInt& J, ConstRefRingElem x)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long i = ConvertTo<long>(I,  ErrMesg);
    const long j = ConvertTo<long>(J,  ErrMesg);
    SetEntry(M, i-1, j-1, x);
  }


  std::vector<RingElem> GetRow_forC5(ConstMatrixView M, const BigInt& I)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long i = ConvertTo<long>(I,  ErrMesg);
    return GetRow(M, i-1);
  }


  std::vector<RingElem> GetCol_forC5(ConstMatrixView M, const BigInt& J)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long j = ConvertTo<long>(J,  ErrMesg);
    return GetCol(M, j-1);
  }


  void SetRow_forC5(MatrixView& M, const BigInt& I, const std::vector<RingElem>& v)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long i = -1 + ConvertTo<long>(I,  ErrMesg);
    if (len(v) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    for (long j=0; j<NumCols(M); ++j)  SetEntry(M, i, j, v[j]);
  }


  void SetCol_forC5(MatrixView& M, const BigInt& J, const std::vector<RingElem>& v)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long j = -1 + ConvertTo<long>(J,  ErrMesg);
    if (len(v) != NumRows(M))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    for (long i=0; i<NumRows(M); ++i)  SetEntry(M, i, j, v[i]);
  }


  void SwapRows_forC5(matrix& M, const BigInt& row1, const BigInt& row2)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long r1 = -1 + ConvertTo<long>(row1,  ErrMesg);
    const long r2 = -1 + ConvertTo<long>(row2,  ErrMesg);
    SwapRows(M, r1, r2);
  }


  void SwapCols_forC5(matrix& M, const BigInt& col1, const BigInt& col2)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    const long c1 = -1 + ConvertTo<long>(col1,  ErrMesg);
    const long c2 = -1 + ConvertTo<long>(col2,  ErrMesg);
    SwapCols(M, c1, c2);
  }


  ConstMatrix IdentityMat_forC5(const ring& R, const BigInt& N)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return IdentityMat(R, ConvertTo<long>(N,  ErrMesg));
  }


  ConstMatrix ZeroMat_forC5(const ring& R, const BigInt& NRows, const BigInt& NCols)
  {
    const ErrorInfo ErrMesg(ERR::ArgTooBig, CoCoA_ERROR_CONTEXT);
    return ZeroMat(R,
                   ConvertTo<long>(NRows,  ErrMesg),
                   ConvertTo<long>(NCols,  ErrMesg));
  }


  ConstMatrixView transposed_forC5(ConstMatrixView M)
  { return transpose(M); }


//   std::vector<RingElem> minors_forC5(const BigInt& N, ConstMatrixView M)
//   {
//     const char* const FnName = "minors";
//     vector<RingElem> v;
//     size_t n;
//     if (!IsConvertible(n, N))
//       CoCoA_THROW_ERROR2(ERR::BadRowIndex, FnName);
//     if (n>NumRows(M))  CoCoA_THROW_ERROR2(ERR::BadRowIndex, FnName);
//     if (n>NumCols(M))  CoCoA_THROW_ERROR2(ERR::BadColIndex, FnName);
//     vector<size_t> rows(n);
//     vector<size_t> cols(n);
// ---- we need a tuple-generator ----
//         v.push_back(det(submat(M,rows,cols)));
//     return v;
//   }
//   }




  //-------- ideal --------------------------------------------------



  //-------- modules --------------------------------------------------

  module NewFreeModule_forC5(const ring& R, ConstMatrixView M)
  {
    std::vector<degree> shifts;
    long n = NumCols(M);
    degree d(n);
    for (long i=0; i<NumRows(M); ++i)
    {
      for (long j=0; j<n; ++j)
        SetComponent(d, j, ConvertTo<BigInt>(M(i,j)));
      shifts.push_back(d);
    }
    return NewFreeModule(R, shifts);
  }


  long NumCompts_forC5(const module& M)
  {
    if (!IsFinGenModule(M))  CoCoA_THROW_ERROR1(ERR::ReqFinGenModule);
    return NumCompts(M);
  }


  ModuleElem NewFreeModuleElem(const module& M, const vector<RingElem>& v)
  {
    if (!IsFreeModule(M))  CoCoA_THROW_ERROR1(ERR::ReqFinGenModule);
    const vector<ModuleElem>& e = gens(M);
    if (len(e) != len(v))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    ModuleElem res(M);
    for (long i=0; i<len(e); ++i)  res += v[i]*e[i];
    return res;
  }


  const std::vector<ModuleElem>& CanonicalBasis_forC5(const module& F)
  {
    if (!IsFreeModule(F))  CoCoA_THROW_ERROR1(ERR::ReqFinGenModule);
    return CanonicalBasis(FreeModule(F));
  }
  

  FinGenModule SubmoduleCols_forC5(const module& F, ConstMatrixView M)
  {
    if (!IsFreeModule(F))  CoCoA_THROW_ERROR1(ERR::ReqFinGenModule);
    return SubmoduleCols(F, M);
  }


  FinGenModule SubmoduleRows_forC5(const module& F, ConstMatrixView M)
  {
    if (!IsFreeModule(F))  CoCoA_THROW_ERROR1(ERR::ReqFinGenModule);
    return SubmoduleRows(F, M);
  }

  
  // Delete me if you want to (by Nico Mexis)
  RingElem GraeffeN_forC5(ConstRefRingElem f, const BigInt& N)
  {
    long n;
    if (!IsConvertible(n, N)) CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return GraeffeN(f, n);
  }

  
  RingElem cyclotomic_forC5(const BigInt& N, ConstRefRingElem x)
  {
    long n;
    if (!IsConvertible(n, N))  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return cyclotomic(n, x);
  }

  
  BigInt CyclotomicIndex_forC5(ConstRefRingElem f)
  { return BigInt(CyclotomicIndex(f)); }

  
  BigInt CyclotomicTest_forC5(ConstRefRingElem f)
  { return BigInt(CyclotomicTest(f)); }

  
  BigInt LRSDegeneracyOrder_forC5(ConstRefRingElem f)
  { return BigInt(LRSDegeneracyOrder(f)); }

  
  std::vector<BigInt> LRSDegeneracyOrders_forC5(ConstRefRingElem f)
  {
// ALMOST    return VectorLong(LRSDegeneracyOrders(f));
    const vector<unsigned long> ords = LRSDegeneracyOrders(f);
    vector<BigInt> ans(len(ords));
    for (long i=0; i < len(ords); ++i)
      ans[i] = BigInt(ords[i]);
    return ans;
  }

  
  std::vector<BigInt> LRSDegeneracyOrdersUnchecked_forC5(ConstRefRingElem f, const BigInt& VerifLevel)
  {
    long Vlev;
    if (VerifLevel < 1)  CoCoA_THROW_ERROR2(ERR::BadArg, "arg2: VerifLevel must be >= 1");
    if (!IsConvertible(Vlev, VerifLevel))
      CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    
// ALMOST    return VectorLong(LRSDegeneracyOrders(f,VerificationLevel(Vlev)));
    const vector<unsigned long> ords = LRSDegeneracyOrders(f, VerificationLevel(Vlev));
    vector<BigInt> ans(len(ords));
    for (long i=0; i < len(ords); ++i)
      ans[i] = BigInt(ords[i]);
    return ans;
  }
  

  bool IsLRSDegenerateOrder_forC5(ConstRefRingElem f, const BigInt& N)
  {
    unsigned long n;
    if (!IsConvertible(n, N)) CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    return IsLRSDegenerateOrder(f, n);
  }


  bool IsLRSDegenerate_forC5(ConstRefRingElem f)
  { return IsLRSDegenerate(f); }

  
  BigInt MoebiusFn_forC5(const BigInt &N)
  { return BigInt(MoebiusFn(ConvertTo<long>(N))); }

}

