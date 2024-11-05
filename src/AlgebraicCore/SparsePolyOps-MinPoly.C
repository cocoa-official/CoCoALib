//   Copyright (c)  2016-2017 John Abbott, Anna M. Bigatti
//   Author: 2016 John Abbott, Anna M. Bigatti, Elisa Palezzato

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


#include "CoCoA/SparsePolyOps-MinPoly.H"

#include "CoCoA/BigIntOps.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/LinDepMill.H"
#include "CoCoA/MatrixOps.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/NumTheory-prime.H" // for NextPrime
#include "CoCoA/PPMonoidHom.H"  // for PPMonoidHom, GeneralHom
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H" // for CanonicalRepr / DefiningIdeal
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"  // for IsQQ
#include "CoCoA/RingZZ.H"
#include "CoCoA/SmallFpImpl.H" // for SmallFpImpl::ourMaxModulus
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-ideal.H" // for QuotientBasisSorted
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H
//#include "CoCoA/SparsePolyRing.H" // for HasGBasis
#include "CoCoA/VectorOps.H" // for HasUniqueOwner
#include "CoCoA/VerificationLevel.H"  // for MinPolyModular
#include "CoCoA/assert.H"
#include "CoCoA/bool3.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H" // for CheckForInterrupt
#include "CoCoA/time.H"
#include "CoCoA/utils.H" // for len
#include "CoCoA/verbose.H"  // for VerboseLog

#include <algorithm>
using std::min;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;
//#include <map>  // for QBPositions, but seems worse
//using std::insert;

namespace CoCoA
{

  namespace // anonymous
  {

    // >>>ASSUMES<<< that QB is in increasing term-order!
    void coefficients(vector<RingElem>& coeffs, ConstRefRingElem f, const vector<PPMonoidElem>& QB)
    {
      if (len(coeffs)!=len(QB))
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "len(coeffs)!=len(QB)");
      SparsePolyIter itf=BeginIter(f);
      long i=len(QB)-1;
      for (; i>=0 && !IsEnded(itf); --i)
      {
        CoCoA_ASSERT(PP(itf) <= QB[i]);
        if (PP(itf) == QB[i])
        {
          coeffs[i] = coeff(itf);
          ++itf;
        }
        else
          coeffs[i] = 0;
      }
      for (; i>=0; --i)  coeffs[i] = 0;
    }


    // >>>ASSUMES<<< that QB is in increasing term-order!
    void CoefficientsInMatCol(MatrixView& M, long j,
                              ConstRefRingElem f, const vector<PPMonoidElem>& QB)
    {
      if (NumRows(M)!=len(QB))
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "NumCols(M)!=len(QB)");
      //for (long i=0; i<NumRows(M); ++i) SetEntry(M,i,j, CoeffOfTerm(p,QB[i]));
      SparsePolyIter itf=BeginIter(f);
      long i=len(QB)-1;
      for (; i>=0 && !IsEnded(itf); --i)
      {
        CoCoA_ASSERT_ALWAYS(PP(itf) <= QB[i]);
        if (PP(itf) == QB[i])
        {
          SetEntry(M,i,j, coeff(itf));
          ++itf;
        }
        else
          SetEntry(M,i,j, 0);
      }
      for (; i>=0; --i)  SetEntry(M,i,j, 0);
    }


//const std::map<PPMonoidElem, long> QBPositions(const vector<PPMonoidElem>& QB)
//{
//  std::map<PPMonoidElem, long> res;
//  for (long i=len(QB)-1; i>=0; --i)  res.insert(std::make_pair(QB[i], i));
//  return res;
//}


    long position(ConstRefPPMonoidElem elem, vector<PPMonoidElem> QB)
    {
      //for (long i=1; i<len(QB); ++i)  // wrong way!
      for (long i=len(QB)-1; i>=0; --i)  if (elem==QB[i]) return i;
      CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere); // only for elem in QB
      return 0;  // just to keep compiler quiet
    }


    //long position(ConstRefPPMonoidElem elem, const std::map<PPMonoidElem,long>& QBPosn)
    //{ return QBPosn[elem]; }


    matrix MultiplicationMat(ConstRefRingElem f_orig, const ideal& I)
    {
      VerboseLog VERBOSE("MultiplicationMat");
      const PolyRing P(RingOf(I), CoCoA_ERROR_CONTEXT);
      if (P != owner(f_orig))  CoCoA_THROW_ERROR1(ERR::MixedRings);
      if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
      const vector<PPMonoidElem>& x_PPM = indets(PPM(P));
      const vector<RingElem>& x_P = indets(P);
      const vector<PPMonoidElem> QB = QuotientBasisSorted(I);
      const long lenQB = len(QB);
      //      std::map<PPMonoidElem, long> QBPosn = QBPositions(QB);
      RingElem f = NF(f_orig, I);
      RingElem p;
      long MinDegI = deg(GBasis(I)[0]);
      VERBOSE(90) << lenQB<<" x "<< lenQB<< std::endl;
      matrix M = NewDenseMat(CoeffRing(P), lenQB, lenQB);
      for (long j=0; j<lenQB; ++j)
      {
        VERBOSE(99) << j<< std::endl;
        //        if (deg(QB[j])<12)
        if (deg(QB[j]) < 3*MinDegI)
          p = NF(f*monomial(P, QB[j]), I);
        else
        {
          //std::cout << "*";
          long k = NumIndets(P);
          for (--k; k>=0; --k)  if (IsDivisible(QB[j], x_PPM[k])) break;
          //long quot_pos = position(QB[j]/x_PPM[k], QB);
          long quot_pos = position(QB[j]/x_PPM[k], QB);
          RingElem NFquot = zero(P); // scalar product:
          for (long d=0; d<lenQB; ++d)
            NFquot += monomial(P, M(d,quot_pos), QB[d]);
          p = NF(NFquot*x_P[k], I);
        }
        //for (long i=0; i<len(QB); ++i)  SetEntry(M,i,j, CoeffOfTerm(p,QB[i]));
        CoefficientsInMatCol(M,j, p, QB);
      }
      return M;
    }


    RingElem NFPower(ConstRefRingElem t, long exp, const ideal& I)
    {
      //// ring compatible?
      if (exp<0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
      if (exp==0) return one(owner(t));
      if (exp==1) return t;
      if (IsEven(exp))
        return NFPower(NF(power(t,2), I), exp/2, I);
      return NF(NFPower(NF(power(t,2), I), exp/2, I)*t, I);
    }


    namespace // anonymous
    {

      matrix ChangeBasis(const ideal& I, const std::vector<PPMonoidElem>& QB2)
      {
        const std::vector<PPMonoidElem> QB = QuotientBasisSorted(I);
        const ring P = RingOf(I);
        matrix M = NewDenseMat(CoeffRing(P), len(QB), len(QB));
        for (long j=0; j<len(QB2); ++j)
        {
          RingElem p = NF(monomial(P, QB2[j]), I);
          //for (long i=0; <len(QB);++i)  SetEntry(M,i,j,CoeffOfTerm(p,QB[i]));
          CoefficientsInMatCol(M, j, p, QB);
        }
        if (IsZeroDet(M))  CoCoA_THROW_ERROR2(ERR::BadArg, "QB2 is not a basis");
        return M;
      }

    }

  } // anonymous namespace





  matrix FrobeniusMat(const ideal& I)
  {
    VerboseLog VERBOSE("FrobeniusMat");
    ring P = RingOf(I);
    long p = ConvertTo<long>(characteristic(P));
    if (p==0) CoCoA_THROW_ERROR2(ERR::BadRing, "Characteristic must be finite");
    long e = LogCardinality(CoeffRing(P));
    if (p==0)  CoCoA_THROW_ERROR2(ERR::BadRing, "Field must be finite");

    const vector<PPMonoidElem> QB = QuotientBasisSorted(I);
    const long lenQB = len(QB);
    //    std::map<PPMonoidElem, long> QBPosn = QBPositions(QB);
    matrix M = NewDenseMat(ZeroMat(CoeffRing(P), lenQB, lenQB));
    SetEntry(M, 0,0, one(CoeffRing(P)));

    RingElem NFj(P);
    RingElem NFk(P);
    RingElem NFquot(P);
    double U = CpuTime();
    for (long j=1; j<lenQB; ++j)
    {
      long k = j-1;
      for (; k>0; --k)  if (IsDivisible(QB[j], QB[k])) break;
      if (k==0)
        NFj = NFPower(monomial(P, QB[j]), p, I);
      else
      {
        const long quot_pos = position(QB[j]/QB[k], QB);
        NFk = zero(P);
        NFquot = zero(P);
        for (long d=0; d<lenQB; ++d)
        { // scalar products:
          NFk    += monomial(P, M(d,k),        QB[d]);
          NFquot += monomial(P, M(d,quot_pos), QB[d]);
        }
        NFj = NF(NFk*NFquot, I);
      }
      // coefficients
      //for (long i=0;i<NumRows(M);++i) SetEntry(M,i,j,CoeffOfTerm(NFj,QB[i]));
      CoefficientsInMatCol(M, j, NFj, QB);
    }
    VERBOSE(90) << "cycle time --> " << CpuTime() - U << std::endl;
    //    VERBOSE(99) << "M = " << M << std::endl;
    return power(M, e);
  }


  matrix FrobeniusMat(const ideal& I, const std::vector<PPMonoidElem>& QB2)
  {
    matrix FM = FrobeniusMat(I);
    matrix BC = ChangeBasis(I,QB2);
    return inverse(BC)*FM*BC;
  }


  RingElem MinPoly(ConstMatrixView M, ConstRefRingElem x)
  {
    if (NumRows(M)!=NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (!IsIndet(x))  CoCoA_THROW_ERROR2(ERR::ReqIndet, "2nd arg");

    const ring& RM = RingOf(M); // we shall check that this is CoeffRing(Rx)
    const PolyRing Rx(owner(x), CoCoA_ERROR_CONTEXT);

    // TODO: convert M into matrix of coefficients (if possible)
    if (RM==Rx)  CoCoA_THROW_ERROR2(ERR::NYI, "RingOf(M)==Rx");
    const ring& K = CoeffRing(Rx);
    //   If (IsZZ(R) or IsQQ(R)) Then
    //     M := CanonicalHom(R, K)( M );
    //     R := K;
    //   EndIf;
    if (K!=RM)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    long NRows = NumRows(M);

    LinDepMill ILD(K, NRows*NRows);
///    vector<RingElem> flattenv(NRows*NRows, zero(K));

    matrix MD = NewDenseMat(IdentityMat(K, NRows));
    vector<RingElem> flattenv = FlattenByRows(MD);
    ILD.myAppendVec(flattenv);

    /*C++17 [[maybe_unused]]*/    double S1=0, S2=0, S5=0, V;
    while ( (ILD.myLinReln()).empty() )
    {
      V = CpuTime(); MD = M*MD;                    S1 += (CpuTime()-V);
      V = CpuTime(); flattenv = FlattenByRows(MD); S2 += (CpuTime()-V);
      V = CpuTime(); ILD.myAppendVec(flattenv);    S5 += (CpuTime()-V);
    }

    double Y = CpuTime();
    const RingHom& phi = CoeffEmbeddingHom(Rx);
    RingElem mi = zero(owner(x));
    const vector<RingElem>& linreln = ILD.myLinReln();
    for (long d=0; d<len(linreln); ++d)  mi += phi(linreln[d])*power(x,d);
    Y = CpuTime() - Y;

    return mi;
  }


//   RingElem MinPolyQuotMat2(ConstRefRingElem L, const ideal& I, ConstRefRingElem x)
//   {
//     return MinPoly(MultiplicationMat(L,I), x);
//   }


  // RingElem CharPolyQuot(ConstRefRingElem f_orig, const ideal& I, ConstRefRingElem x)
  // {
  //   VerboseLog VERBOSE("CharPolyQuot");
  //   //VERBOSE(99) << f_orig << std::endl;
  //   double U;
  //   RingElem f = NF(f_orig, I);
  //   U = CpuTime();
  //   matrix M = MultiplicationMat(f,I);
  //   VERBOSE(90) << "Time MultiplicationMat: " << (CpuTime() - U) << std::endl;

  //   U = CpuTime();

  //   const long LenQB = len(QuotientBasis(I)); // multiplicity
  //   M = CanonicalHom(RingOf(M),owner(x))(M);
  //   M -= DiagMat(std::vector<RingElem>(LenQB, x));
  //   return det(M);
  // }


  RingElem MinPolyQuotMat(ConstRefRingElem f_orig, const ideal& I, ConstRefRingElem x)
  {
    VerboseLog VERBOSE("MinPolyQuotMat");
    //VERBOSE(99) << f_orig << std::endl;
    double U;
    RingElem f = NF(f_orig, I);
    U = CpuTime();
    matrix M = MultiplicationMat(f,I);
    VERBOSE(90) << "Time MultiplicationMat: " << (CpuTime() - U) << std::endl;

    U = CpuTime();

    const PolyRing P(owner(x), CoCoA_ERROR_CONTEXT);
    const ring& K = CoeffRing(P);
    //    const vector<PPMonoidElem> QB = QuotientBasisSorted(I);
    const long LenQB = len(QuotientBasis(I)); // multiplicity
    LinDepMill ILD(K, LenQB);
////    vector<RingElem> FlatVec(LenQB, zero(K));
    matrix fPow = NewDenseMat(ZeroMat(K, LenQB, 1));
    SetEntry(fPow, 0,0, one(K));

    vector<RingElem> FlatVec = FlattenByRows(fPow);
    ILD.myAppendVec(FlatVec);

    U = CpuTime();
    double S1 = 0;
    double S2 = 0;
    double S5 = 0;
    double V;
    while ( (ILD.myLinReln()).empty() )
    {
      CheckForInterrupt("MinPolyQuotMat");
      V = CpuTime(); fPow = M*fPow;                 S1 += (CpuTime()-V);
      V = CpuTime(); FlatVec = FlattenByRows(fPow); S2 += (CpuTime()-V);
      V = CpuTime(); ILD.myAppendVec(FlatVec);      S5 += (CpuTime()-V);
    }
    VERBOSE(90) << "time cycle: --> " << (CpuTime() - U) << std::endl;

    double Y = CpuTime();
    RingHom phi = CoeffEmbeddingHom(owner(x));
    //ScalarProduct(col, [x^d | d in 0..i-1]);
    RingElem mp = zero(owner(x));
    const vector<RingElem>& linreln = ILD.myLinReln();
    for (long d=0; d<len(linreln); ++d)  mp += phi(linreln[d])*power(x,d);
    Y = CpuTime() - Y;

    VERBOSE(90) << "  time: M*fPow --> " << S1 << std::endl
                << "  time: flatten --> " << S2 << std::endl
                << "  time: LinDepMill --> " << S5 << std::endl
                << "Post-cycle --> " << Y << std::endl;
    return mp;
  }


  RingElem MinPolyQuotElim(ConstRefRingElem L, const ideal& I, ConstRefRingElem x)
  {
    VerboseLog VERBOSE("MinPolyQuotElim");
    //VERBOSE(99) << L << std::endl;
    // controllare anelli
    ring P = RingOf(I);
    long n = NumIndets(P);
    //  myIndets = concat(IndetSymbols(P), [Record [head ="aux", indices =[]]]);
    ring S = NewPolyRing(CoeffRing(P), NewSymbols(n+1));
    vector<RingElem> indetS_1;
    for (long i=0; i<n; ++i) indetS_1.push_back(indet(S,i));
    RingHom phi = PolyAlgebraHom(P, S, indetS_1);
    vector<RingElem> zero_x;
    for (long i=0; i<n; ++i) zero_x.push_back(zero(P));
    zero_x.push_back(x);
    RingHom psi = PolyAlgebraHom(S, P, zero_x);
    RingElem s = indet(S, n); // auxiliary variable
    ideal J = elim(ideal(phi(gens(I))) + ideal(s-phi(L)), indetS_1);
    return psi(monic(gens(J)[0]));
  }


  RingElem MinPolyQuotDef(ConstRefRingElem f_orig, const ideal& I, ConstRefRingElem x)
  {
    VerboseLog VERBOSE("MinPolyQuotDef");
    //  const double T = CpuTime();
    const PolyRing P(owner(f_orig), CoCoA_ERROR_CONTEXT);
    const PolyRing P_x(owner(x), CoCoA_ERROR_CONTEXT);
    const ring& K = CoeffRing(P);
    // fare tutti i controlli sugli anelli
    if (K != CoeffRing(P_x))  CoCoA_THROW_ERROR2(ERR::BadRing, "incompatible coeff rings");
    if (IsOne(I))  CoCoA_THROW_ERROR2(ERR::BadArg, "Ideal is (1)");
    //    std::cout << x << std::endl;
    //    std::cout << f_orig << std::endl;
    //    std::cout << (owner(f_orig) == RingOf(I)) << std::endl;
    RingElem f = NF(f_orig, I);
    //    std::cout << f << std::endl;
    vector<PPMonoidElem> QB = QuotientBasisSorted(I); // 2018-01 AMB
    //    std::cout << "len(QB) " << len(QB) << std::endl;
    LinDepMill ILD(K, len(QB));
    //    std::cout << "after ILD " << std::endl;
    vector<RingElem> coeffs(len(QB), zero(K));

    RingElem p = one(P);
    coeffs[0] = one(K);
    ILD.myAppendVec(coeffs);

    double W = CpuTime();
    double S1 = 0;
    double S2 = 0;
    double S3 = 0;
    double S4 = 0;
    double U;
    while ( (ILD.myLinReln()).empty() )
    {
      CheckForInterrupt("MinPolyQuotDef");
      U = CpuTime();  p *= f;  S1 = S1+(CpuTime()-U);
      U = CpuTime();  p = NF(p, I);  S2 = S2+(CpuTime()-U);
      U = CpuTime();  coefficients(coeffs, p, QB);  S3 = S3+(CpuTime()-U);
      U = CpuTime();  ILD.myAppendVec(coeffs);  S4 = S4+(CpuTime()-U);
    }
    W = CpuTime() - W;
    VERBOSE(90) << "time cycle: --> " << (W) << std::endl;

    double Y = CpuTime();
    RingHom phi = CoeffEmbeddingHom(P_x);
    //ScalarProduct(col, [x^d | d in 0..i-1]);
    RingElem mp = zero(P_x);
    const vector<RingElem>& linreln = ILD.myLinReln();
    for (long d=0; d<len(linreln); ++d)  mp += phi(linreln[d])*power(x,d);
    Y = CpuTime() - Y;

    VERBOSE(90)  << "Total time --> " << S1+S2+S3+S4 << std::endl
                 << "  time: p*f --> " << S1 << std::endl
                 << "  time: NF --> " << S2 << std::endl
                 << "  time: coefficients --> " << S3 << std::endl
                 << "  time: LinDepMill --> " << S4 << std::endl
                 << "Post-cycle --> " << Y << std::endl;
    return monic(mp);
  }


  RingElem MinPolyQuotDefLin(ConstRefRingElem f_orig, const ideal& I, ConstRefRingElem x)
  {
    VerboseLog VERBOSE("MinPolyQuotDefLin");
    //  const double T = CpuTime();
    const PolyRing P(owner(f_orig), CoCoA_ERROR_CONTEXT);
    const PolyRing P_x(owner(x), CoCoA_ERROR_CONTEXT);
    const ring& K = CoeffRing(P);
    RingHom psi = CoeffEmbeddingHom(P);
    // fare tutti i controlli sugli anelli
    if (K != CoeffRing(P_x))  CoCoA_THROW_ERROR2(ERR::BadRing, "incompatible coeff rings");
    if (IsOne(I))  CoCoA_THROW_ERROR2(ERR::BadArg, "Ideal is (1)");
    RingElem f = NF(f_orig, I);
    vector<PPMonoidElem> QB = QuotientBasisSorted(I); // 2018-01 AMB
    //  std::map<PPMonoidElem, long> QBPosn = QBPositions(QB);
    LinDepMill ILD(K, len(QB));

    //VERBOSE(99) << f_orig << std::flush;
    //---- "MultiplicationMat"  fxQB
    double U = CpuTime();
    const vector<RingElem>& x_P = indets(P);
    const vector<PPMonoidElem>& x_PPM = indets(PPM(P));
    vector<RingElem> fxQB(len(QB));
    for (long j=0; j<len(QB); ++j)
    {
      if (deg(QB[j])<12)
        fxQB[j] = NF(f*monomial(P, QB[j]), I);
      else
      {
        long k = NumIndets(P);
        for (--k; k>=0; --k)  if (IsDivisible(QB[j], x_PPM[k])) break;
        fxQB[j] = NF(fxQB[position(QB[j]/x_PPM[k], QB)]*x_P[k], I);
      }
    }
    VERBOSE(90) << "Time MultiplicationMat: " << (CpuTime() - U) << std::endl;
    //---- MultiplicationMat - end

    vector<RingElem> coeffs(len(QB), zero(K));

    RingElem p = one(P);
    coeffs[0] = one(K);
    ILD.myAppendVec(coeffs);
    double W = CpuTime();
    double S1 = 0;
    double S2 = 0;
    double S3 = 0;
    double S4 = 0;
    RingElem partial = zero(P);
    while ( (ILD.myLinReln()).empty() )
    {
      CheckForInterrupt("MinPolyQuotDefLin");
      U = CpuTime();  //p = p*f;
      partial = zero(P);
      for (SparsePolyIter itp=BeginIter(p); !IsEnded(itp); ++itp)
        partial += fxQB[position(PP(itp), QB)]*psi(coeff(itp));
      p = partial;
      S1 = S1+(CpuTime()-U);
      //    U = CpuTime();  p = NF(p, I);  S2 = S2+(CpuTime()-U);
      U = CpuTime();  coefficients(coeffs, p, QB);  S3 = S3+(CpuTime()-U);
      U = CpuTime();  ILD.myAppendVec(coeffs);  S4 = S4+(CpuTime()-U);
    }
    W = CpuTime() - W;
    VERBOSE(90) << "  ... MinPolyDefQuot: cycle " << W << std::endl;

    double Y = CpuTime();
    RingHom phi = CoeffEmbeddingHom(P_x);
    //ScalarProduct(col, [x^d | d in 0..i-1]);
    RingElem mp = zero(P_x);
    const vector<RingElem>& linreln = ILD.myLinReln();
    for (long d=0; d<len(linreln); ++d)  mp += phi(linreln[d])*power(x,d);
    Y = CpuTime() - Y;

    VERBOSE(90)  << "Time in cycle --> " << S1+S2+S3+S4 << std::endl
                 << "  p*f --> " << S1 << std::endl
                 << "  NF --> " << S2 << std::endl
                 << "  coefficients --> " << S3 << std::endl
                 << "  LinDepMill --> " << S4 << std::endl
                 << "After cycle --> " << Y << std::endl;
    return monic(mp);
  }

  
  RingElem MinPolyModular(ConstRefRingElem f, const ideal& I,
                          ConstRefRingElem x, VerificationLevel VerLev);
  

  RingElem MinPolyQuot(ConstRefRingElem f, const ideal& I, ConstRefRingElem x, VerificationLevel VerLev)
  {
    if (!IsIndet(x))  CoCoA_THROW_ERROR2(ERR::ReqIndet, "3rd arg");
    const SparsePolyRing& P = RingOf(I);
    if (owner(f) != P)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (CoeffRing(owner(x)) != CoeffRing(P))  CoCoA_THROW_ERROR1(ERR::MixedCoeffRings);
    ideal J(zero(P));
    if (HasGBasis(I) || IsGuaranteed(VerLev))
    {
      vector<RingElem> RGB = ReducedGBasis(I);
      vector<long> IndetsUsed = IndetsIn(RGB);
      const int nvars = NumIndets(P);
      vector<bool> seen(nvars);
      for (long i: IndetsUsed)
      {
        seen[i] = true;
      }
      for (long i: IndetsIn(f))
      {
        seen[i] = true;
      }
      vector<long> IndetsNotUsed;
      for (int i=0; i < nvars; ++i)
        if (!seen[i])
          RGB.push_back(indet(P,i));
      J = ideal(RGB);      
    }
    else // !HasGBasis(I)
    {
      vector<RingElem> G = gens(I);
      vector<long> IndetsUsed = IndetsIn(G);
      const int nvars = NumIndets(P);
      vector<bool> seen(nvars);
      for (long i: IndetsUsed)
      {
        seen[i] = true;
      }
      for (long i: IndetsIn(f))
      {
        seen[i] = true;
      }
      vector<long> IndetsNotUsed;
      for (int i=0; i < nvars; ++i)
        if (!seen[i])
          G.push_back(indet(P,i));
      J = ideal(G);
    }
    if (IsQQ(CoeffRing(owner(f))))  return MinPolyModular(f, J, x, VerLev);
    return MinPolyQuotDef(f, J, x);
  }


  //--------- MinPoly in QuotientRing -----------------------

  RingElem MinPolyMat(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (!IsQuotientRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::NotElemQuotientRing);
    return MinPolyQuotMat(CanonicalRepr(f), DefiningIdeal(owner(f)), x);
  }

  RingElem MinPolyDef(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (!IsQuotientRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::NotElemQuotientRing);
    return MinPolyQuotDef(CanonicalRepr(f), DefiningIdeal(owner(f)), x);
  }

  RingElem MinPolyElim(ConstRefRingElem f, ConstRefRingElem x)
  {
    if (!IsQuotientRing(owner(f)))  CoCoA_THROW_ERROR1(ERR::NotElemQuotientRing);
    return MinPolyQuotElim(CanonicalRepr(f), DefiningIdeal(owner(f)), x);
  }


  RingElem MinPoly(ConstRefRingElem f, ConstRefRingElem x)
  { return MinPolyDef(f, x); }


//---------- MinPoly modular -----------------------------------------

  RingElem LiftPolyFromFpxToQQx(const PolyRing& QQx, ConstRefRingElem f)
  {
    const SparsePolyRing& Fpx(owner(f));
    const PPMonoidHom psi(GeneralHom(PPM(Fpx), indets(PPM(QQx))));
    RingElem ans(QQx);
    for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
      ans += monomial(QQx, ConvertTo<BigInt>(coeff(it)), psi(PP(it)));
    return ans;
  }


  namespace //  anonymous
  {
    RingElem UniPolyToQQ(ConstRefRingElem f, ConstRefRingElem x)
    {
      const SparsePolyRing& QQx = owner(x);
      const SparsePolyRing& Fpx(owner(f));
      const PPMonoidHom psi(GeneralHom(PPM(Fpx), std::vector<PPMonoidElem>(NumIndets(Fpx),LPP(x))));
      RingElem ans(QQx);
      for (SparsePolyIter it=BeginIter(f) ; !IsEnded(it) ; ++it )
        ans += monomial(QQx, ConvertTo<BigInt>(coeff(it)), psi(PP(it)));
      return ans;
    }

    vector<RingElem> PolyVectorToQQ(const PolyRing QQx, const vector<RingElem>& v_p)
    {
      vector<RingElem> v_QQ;
      for (long i=0; i<len(v_p); ++i)
        v_QQ.push_back(LiftPolyFromFpxToQQx(QQx, v_p[i]));
      return v_QQ;
    }

  }


  namespace // anonymous
  {

    // NF(f(g), I) == 0?  (using Horner)
    bool IsZeroEvalUniPolyMod(ConstRefRingElem f, ConstRefRingElem g, const ideal& I)
    {
      if (!HasGBasis(I))  CoCoA_THROW_ERROR2(ERR::BadArg, "No precomputed GBasis: cannot verify");
      const ring& Pf = owner(f);
      const ring& Pg = owner(g);
      if (CoeffRing(Pf) != CoeffRing(Pg))  CoCoA_THROW_ERROR2(ERR::MixedCoeffRings, "arg 1 and 2");
      if (Pg != RingOf(I))  CoCoA_THROW_ERROR2(ERR::MixedRings, "arg 2 and 3");
      const RingElem z = indet(Pf, UnivariateIndetIndex(f));
      //  const RingHom phi = CoeffEmbeddingHom(Pg);
      //      VERBOSE(50) << "horner-"<<i << std::endl;
      const std::vector<RingElem> c = CoeffVecWRT(f, z);
      const RingHom phi = CoeffEmbeddingHom(Pg);
      //      const long NF_freq = (NumIndets(Pg)>2)?1:2;
      RingElem eval_f = zero(Pg);
      for (long d=len(c)-1; d>=0; --d)
      {
        eval_f *= g;
        if (!IsZero(c[d]))  eval_f += phi(LC(c[d]));
        //    if (d%NF_freq == 0)  // NF_freq = 1 or 1?  Or better not?
        eval_f = NF(eval_f, I);
      }
      return IsZero(eval_f);
    }


    // RingElem EvalUniPolyModTF(f, g, I)    // removed - not always good

    bool IsZeroEvalUniPolyModFp(ConstRefRingElem f, ConstRefRingElem g,
                                const ideal& I, VerificationLevel VerLev)
    {
      VerboseLog VERBOSE("EvalUniPolyModTF");
      const ring& Pf = owner(f);
      const ring& Pg = owner(g);
      if (CoeffRing(Pf) != CoeffRing(Pg))  CoCoA_THROW_ERROR2(ERR::MixedCoeffRings, "arg 1 and 2");
      if (Pg != RingOf(I))  CoCoA_THROW_ERROR2(ERR::MixedRings, "arg 2 and 3");
      const RingElem z = indet(Pf, UnivariateIndetIndex(f));
      const std::vector<RingElem> c = CoeffVecWRT(f, z);
      long p = 0;
      RingElem eval_f;
      for (long i=0; i < level(VerLev); ++i)
      {
        p = RandomSmallPrime((1UL<<31)-1);
        if (HasGBasis(I)) 
          while (p < 1222333444 and !IsSigmaGoodPrime(p,I))
            p = RandomSmallPrime((1UL<<31)-1);
        else
          while (p < 1222333444 and IsDivisible(CommonDenom(gens(I)), p))
            p = RandomSmallPrime((1UL<<31)-1);
        VERBOSE(80) << "Testing prime "<< p << std::endl;
        SparsePolyRing Pp = NewPolyRing(NewZZmod(p), symbols(PPM(Pg)), ordering(PPM(Pg)));
        RingHom pi = PolyRingHom(Pg, Pp, QQEmbeddingHom(Pp), indets(Pp));
        vector<RingElem> Gp;
        if (HasGBasis(I)) Gp = pi(GBasis(I));
        else Gp = GBasis(ideal(pi(gens(I))));
        RingElem gp = pi(g);
        //      VERBOSE(50) << "horner-"<<i << std::endl;
        eval_f = zero(Pp);
        for (long d=len(c)-1; d>=0; --d)
        {
          eval_f *= gp;
          if (not(IsZero(c[d])))  eval_f += QQEmbeddingHom(Pp)(LC(c[d]));
          // if (d%N == 0)  // NF every N steps not good in general
          eval_f = NR(eval_f, Gp);
        }
        if (!IsZero(eval_f))  return false;
      }
      return true;
    }

  } // anonymous namespace


  RingElem MinPolyModular(ConstRefRingElem f, const ideal& I,
                          ConstRefRingElem x, VerificationLevel VerLev)
  {
    VerboseLog VERBOSE("MinPolyModular");
    const SparsePolyRing QQx = RingOf(I);
    if (!IsQQ(CoeffRing(QQx)))  CoCoA_THROW_ERROR2(ERR::BadRing, "CoeffRing for I must be QQ");
    if (!IsQQ(CoeffRing(owner(x))))  CoCoA_THROW_ERROR2(ERR::BadRing, "CoeffRing for indet must be QQ");
    std::vector<RingElem> G;
    if (HasGBasis(I) || IsGuaranteed(VerLev))  // GB needed for verification
    {
      if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
      G = ReducedGBasis(I);
    }
    else
    {
      VERBOSE(20) << "WARNING: no precomputed GBasis, using gens!" << std::endl;
      G = gens(I);
    }
    long PrimeCount = 0;
    PrimeSeqForCRT primes;

    SmallPrime p = *primes; // first prime in list (will be discarded)
    const vector<symbol> syms = NewSymbols(NumIndets(QQx));
    RingElem mpCRT(owner(x));
    RingElem GDenom = CommonDenom(G[0]);
    for (long i=1; i<len(G); ++i) GDenom = lcm(GDenom, CommonDenom(G[i]));
    while (IsZero(mpCRT))
    {
      CheckForInterrupt("MinPolyModular");
      p = NextPrime(primes);  ++PrimeCount;
      VERBOSE(80) << PrimeCount << ": prime is " << p << std::endl;
      if (IsDivisible(GDenom, p))
      {
        VERBOSE(80) << "   UGLY PRIME: going to another prime" << std::endl;
        continue;
      }
      SparsePolyRing Fpx = NewPolyRing(NewZZmod(p), syms);
      RingHom pi = PolyRingHom(QQx, Fpx, CanonicalHom(RingQQ(), Fpx), indets(Fpx));
      ideal Ip = ideal(pi(G));
      if ((!HasGBasis(I)) && (!IsZeroDim(Ip)))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
      RingElem xFp = indet(Fpx,0);
      mpCRT = UniPolyToQQ(MinPolyQuot(pi(f), Ip, xFp), x);
      //VERBOSE(99) << "mpCRT = " << mpCRT << std::endl;
    }
    RingElem mp;
    BigInt mpModulus(p);
    while (true)
    {
      CheckForInterrupt("MinPolyModular");
      p = NextPrime(primes);  ++PrimeCount;
      VERBOSE(80) << PrimeCount << ": prime is " << p << std::endl;
      if (IsDivisible(GDenom, p))
      {
        VERBOSE(80) << "   UGLY PRIME: going to another prime" << std::endl;
        continue;
      }
      SparsePolyRing Fpx = NewPolyRing(NewZZmod(p), syms);
      RingHom pi = PolyRingHom(QQx, Fpx, CanonicalHom(RingQQ(), Fpx), indets(Fpx));
      ideal Jp = ideal(pi(G));
      if ((!HasGBasis(I)) && (!IsZeroDim(Jp)))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
      mp = MinPolyQuot(pi(f), Jp, indet(Fpx,0));
      if (!HasGBasis(I))
        CRTPoly(mpCRT,mpModulus,  mpCRT,mpModulus,  UniPolyToQQ(mp,x),BigInt(p));
      else
      {
        if (deg(mp) == deg(mpCRT))
          CRTPoly(mpCRT,mpModulus, mpCRT,mpModulus, UniPolyToQQ(mp,x),BigInt(p));
        else
        {
          if (deg(mp) < deg(mpCRT))
            VERBOSE(80) << "DISCARDED: deg(mp) < deg(mpCRT)" << std::endl;
          else
          {
            VERBOSE(80) << "SWAPPED: deg(mp) > deg(mpCRT)" << std::endl;
            mpModulus = p;
            mpCRT = UniPolyToQQ(mp,x);
          }
        }
      }
      //VERBOSE(99) << "mpCRT = " << mpCRT << std::endl;
      try
      {
        const double StartReconstr = CpuTime();
        const RingElem Rat_mp = RatReconstructPoly(mpCRT, mpModulus); // may throw
        VERBOSE(90) << "Reconstruction time: " << CpuTime()-StartReconstr << std::endl;

        bool IsVerified;
        double t = CpuTime();
        if (IsGuaranteed(VerLev))
        {
          if (HasGBasis(I))
          {
            IsVerified = IsZeroEvalUniPolyMod(Rat_mp, f, I);
            VERBOSE(90) << "Verification-QQ time:" << CpuTime()-t << std::endl;
          }
          else
          {
            CoCoA_THROW_ERROR2(ERR::ShouldNeverGetHere, "No precomputed GBasis");
            IsVerified = true; // just to keep the compiler quiet
//           VERBOSE(20) << "WARNING: no precomputed GBasis, no exact verification" << std::endl;
          }
        }
        else // >= 0
        {
          IsVerified = IsZeroEvalUniPolyModFp(Rat_mp, f, I, VerLev);
          VERBOSE(90) << "Verification-"<< VerLev <<"Fp time:"
                      << CpuTime()-t << std::endl;
        }
        if (IsVerified) return Rat_mp;
        VERBOSE(80) << "verification FAILED (degree="<< deg(Rat_mp) << ")"
                    << std::endl;
      }
      catch (const CoCoA::ErrorInfo& err)
      {
        if (err != ERR::CannotReconstruct) throw;
      }
    } // while true
  }



  //--------- ShapeLemma -----------------------

  std::vector<RingElem> ShapeLemmaCore(const ideal& I);

/// seems to work
  std::vector<RingElem> ShapeLemma(const ideal& I)
  {
    VerboseLog VERBOSE("ShapeLemma");
    if (!IsQQ(CoeffRing(RingOf(I)))) return ShapeLemmaCore(I);
    //VERBOSE(99) << I << std::endl;
    const SparsePolyRing QQx = RingOf(I);
    //  std::vector<RingElem> SL_p;
    std::vector<RingElem> SL_pQQ;
    std::vector<RingElem> SL_CRT;
    std::vector<RingElem> G;
    if (HasGBasis(I))
    {
      if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
      //    G = ReducedGBasis(I);
      G = GBasis(I);
    }
    else
    {
      VERBOSE(20) << "WARNING: no precomputed GBasis, using gens!" << std::endl;
      G = gens(I);
    }
    long PrimeCount = 1;
    PrimeSeqForCRT primes;
    SmallPrime p = NextPrime(primes);
    VERBOSE(30) << PrimeCount << ": prime is " << p << std::endl;
    SparsePolyRing Fpx = NewPolyRing(NewZZmod(p), NewSymbols(NumIndets(QQx)));
    RingHom pi = PolyRingHom(QQx, Fpx, CanonicalHom(RingQQ(), Fpx), indets(Fpx));
    ideal Ip = ideal(pi(G));
    if ((!HasGBasis(I)) && (!IsZeroDim(Ip)))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
    RingElem xFp = indet(Fpx,0);
    SL_CRT = PolyVectorToQQ(QQx, ShapeLemmaCore(Ip));
    //VERBOSE(99) << "SL_CRT = " << SL_CRT << std::endl;
    BigInt CRTModulus(p);
    while (true)
    {
      CheckForInterrupt("ShapeLemma");
      p = NextPrime(primes);  ++PrimeCount;
      VERBOSE(30) << PrimeCount << ": prime is " << p << std::endl;
      Fpx = NewPolyRing(NewZZmod(p), NewSymbols(NumIndets(QQx)));
      pi = PolyRingHom(QQx, Fpx, CanonicalHom(RingQQ(), Fpx), indets(Fpx));
      Ip = ideal(pi(G));
      if ((!HasGBasis(I)) && (!IsZeroDim(Ip)))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
      vector<RingElem> SL_pQQ = PolyVectorToQQ(QQx, ShapeLemmaCore(Ip));
      if (len(SL_pQQ)!=len(SL_CRT))
        CoCoA_THROW_ERROR2(ERR::IncompatDims, "SL_pQQ & SL_CRT");
      //    VERBOSE(20) << "before CRT " << std::endl;
      BigInt OldCRTModulus(CRTModulus);
      for (long i=0; i<len(SL_pQQ); ++i)
        CRTPoly(SL_CRT[i],CRTModulus, SL_CRT[i],OldCRTModulus, SL_pQQ[i],BigInt(p));
      //    VERBOSE(20) << "CRTModulus = " << CRTModulus << std::endl;
      //    VERBOSE(20) << "SL_CRT = " << SL_CRT << std::endl;
      try
      {
        vector<RingElem> result;
        VERBOSE(50) << "RatReconstructPoly " << std::endl;
        for (long i=0; i<len(SL_pQQ); ++i)
          result.push_back(RatReconstructPoly(SL_CRT[i], CRTModulus));
        return result;
      }
      catch (const CoCoA::ErrorInfo& e)
      {
//        if (message(e) != "cannot reconstruct rational")
//          std::cout << "---- "<< message(e) << std::endl;
      }
    } // while true


  }


  std::vector<RingElem> ShapeLemmaCore(const ideal& I) // FGLM
  {
    VerboseLog VERBOSE("ShapeLemmaCore");
    ring P = RingOf(I);
    ring K = CoeffRing(P);
    RingHom phi = CoeffEmbeddingHom(P);
    // fare tutti i controlli sugli anelli
    if (IsOne(I))  CoCoA_THROW_ERROR2(ERR::BadArg, "Ideal is (1)");
    if (!IsZeroDim(I))  CoCoA_THROW_ERROR1(ERR::ReqZeroDim);
    vector<PPMonoidElem> QB = QuotientBasisSorted(I); // 2018-01 AMB
    vector<RingElem> coeffs(len(QB), zero(K));
    LinDepMill ILD(K, len(QB));

    RingElem z = indet(P, NumIndets(P)-1);
    RingElem p = one(P);
    coeffs[0] = one(K);
    ILD.myAppendVec(coeffs);
    while ( (ILD.myLinReln()).empty() )
    {
      p = NF(p*z, I);
      coefficients(coeffs, p, QB);
      ILD.myAppendVec(coeffs);
    }
    if (len(ILD.myLinReln()) != len(QB)+1)
      CoCoA_THROW_ERROR2(ERR::BadArg, "Last indet not good: len(linreln)!=len(QB)");
    RingElem MinPoly_z(P);
    const vector<RingElem>& linreln = ILD.myLinReln();
    long d;
    for (d=0; d<len(linreln); ++d)  MinPoly_z += phi(linreln[d])*power(z,d); // AMB 2024 remove phi?
    VERBOSE(40) << "Found univariate polynomial" << std::endl;
    //VERBOSE(99) << MinPoly_z << std::endl;
    RingElem prev_f = MinPoly_z -phi(linreln[d-1])*power(z,d-1); /// non serve???
    std::vector<RingElem> FirstIndetsExpr;
    for (long i=0; i<NumIndets(P)-1; ++i)
    {
      p = NF(indet(P, i), I);
      coefficients(coeffs, p, QB);
      ILD.myAppendVec(coeffs);
      RingElem f(P);
      const vector<RingElem>& linreln = ILD.myLinReln();
      for (long d=0; d<len(linreln)-1; ++d)  f += phi(linreln[d])*power(z,d);
      //    FirstIndetsExpr.push_back(f-prev_f+phi(linreln[len(linreln)-1])*indet(P,i));
      FirstIndetsExpr.push_back(indet(P,i)+f); /// + oppure -?
      prev_f = f;
      VERBOSE(40) << "found polynomial for indet " << i << std::endl;
      //VERBOSE(99) << f << std::endl;
    }
    FirstIndetsExpr.push_back(MinPoly_z);
    return FirstIndetsExpr;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-MinPoly.C,v 1.29 2024/09/04 15:05:23 bigatti Exp $
// $Log: SparsePolyOps-MinPoly.C,v $
// Revision 1.29  2024/09/04 15:05:23  bigatti
// Summary: Error codes: MixedSizes renamed IncompatDims
//
// Revision 1.28  2024/08/28 15:21:56  bigatti
// Summary: Error renamed: Not ==> ReqSquareMatrix, ReqFullRank, ReqTermOrdering
//
// Revision 1.27  2024/08/28 13:30:52  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2
//
// Revision 1.26  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.25  2024/07/02 15:36:15  bigatti
// Summary:Changed error codes: LogZero into ReqNonZero
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.24  2024/03/15 18:57:32  abbott
// Summary: Improved use of PolyRing inside fns; employed CoCoA_ERROR_CONTEXT
//
// Revision 1.23  2024/03/15 15:43:47  bigatti
// Summary: using elim;  tested ShapeLemma/FGLM;  using IsZeroDet;  added verbosity
//
// Revision 1.22  2022/11/16 19:46:04  abbott
// Summary: Added commented-out [[maybe unused]]
//
// Revision 1.21  2022/08/05 20:46:46  abbott
// Summary: Took NumRows/NumCols out of for-loop exit conds
//
// Revision 1.20  2022/02/18 14:11:58  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.19  2021/10/04 08:57:01  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.18  2021/02/10 19:41:44  abbott
// Summary: Increased level required to see primes during MinPoly
//
// Revision 1.17  2020/10/02 19:05:26  abbott
// Summary: Replaced MatrixFlatten by FlattenByRows
//
// Revision 1.16  2020/06/17 15:49:27  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.15  2020/03/04 20:06:02  abbott
// Summary: Added 2 comments
//
// Revision 1.14  2020/02/12 09:12:14  bigatti
// -- just a verbosity text
//
// Revision 1.13  2019/12/21 16:45:47  abbott
// Summary: Corrected indentation (I think -- ediff gets confused)
//
// Revision 1.12  2019/03/06 16:01:46  bigatti
// -- fixed: in MinPolyModular guaranteed actually computes GBasis, if missing
//
// Revision 1.11  2019/03/04 16:37:48  abbott
// Summary: Changed comment
//
// Revision 1.10  2019/03/01 13:31:02  bigatti
// -- cout for warning missing verification
//
// Revision 1.9  2018/12/20 15:08:04  bigatti
// -- verbosity to 90 for pinting times
//
// Revision 1.8  2018/12/17 15:16:31  bigatti
// -- code for verification without gbasis
//
// Revision 1.7  2018/09/28 15:54:04  abbott
// Summary: Removed pseudo-ctors NewPolyRing which took just num of indets; now must specify their names
//
// Revision 1.6  2018/05/22 09:53:14  bigatti
// -- fixed MinPoly heuristic verification (mod p)
//
// Revision 1.5  2018/05/22 09:19:16  bigatti
// -- fixed bug for QQ verification
//
// Revision 1.4  2018/05/18 16:38:52  bigatti
// -- added include SparsePolyOps-RingElem.H
//
// Revision 1.3  2018/05/18 12:22:30  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.2  2018/05/17 16:01:49  bigatti
// -- renamed MatrixOperations --> MatrixOps
// -- renamed VectorOperations --> VectorOps
// -- sorted #includes
// -- added #include SparsePolyIter
//
// Revision 1.1  2018/04/06 15:13:26  bigatti
// -- renamed MinPoly.C
//
// Revision 1.54  2018/04/04 12:37:33  bigatti
// -- minor
//
// Revision 1.53  2018/03/20 13:57:35  bigatti
// -- minor adjustment
//
// Revision 1.52  2018/03/15 14:18:06  bigatti
// -- added files SparsePolyOps-ideal.H and SparsePolyOps-involutive.H
//
// Revision 1.51  2018/03/15 10:47:01  abbott
// Summary: Added new fns, IsGuaranteed & level (for VerificationLevel)
//
// Revision 1.50  2018/03/13 18:03:52  bigatti
// -- MinPolyModular: now using VerificationLevel class
//
// Revision 1.49  2018/03/13 17:43:38  bigatti
// -- now MinPolyQuot takes a verification level
//    instead of being called MinPolyQuotHeuristic
//
// Revision 1.48  2018/03/13 17:36:35  abbott
// Summary: Removed old fn CoeffOfTerm
//
// Revision 1.47  2018/03/12 12:50:46  abbott
// Summary: Removed trailing spaces
//
// Revision 1.46  2018/03/09 15:18:42  bigatti
// -- removed map (worse than before)
//
// Revision 1.45  2018/03/09 14:53:33  bigatti
// -- now using coefficients and CoefficientsInMatCol instead of CoeffOfTerm
//    (and must use QuotientBasisSorted)
//
// Revision 1.44  2018/03/02 15:34:48  bigatti
// -- removed pointless verbosity (printing polys)
//
// Revision 1.43  2018/03/02 14:13:15  bigatti
// -- verbosity to 99 for printing polys
//
// Revision 1.42  2018/03/02 14:00:29  bigatti
// -- removed verbosed printing
//
// Revision 1.41  2018/03/01 20:18:20  bigatti
// -- fixed
//
// Revision 1.40  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.39  2018/02/27 16:15:06  abbott
// Summary: Now uses PrimeSeqForCRT instead of NextPrime
//
// Revision 1.38  2018/02/27 15:35:57  abbott
// Summary: Now using PrimeSeqForCRT
//
// Revision 1.37  2018/02/27 10:57:09  abbott
// Summary: Added include NumTheory_prime
//
// Revision 1.36  2018/02/22 17:08:01  abbott
// Summary: Changed initial choice of prime
//
// Revision 1.35  2018/02/22 16:52:21  bigatti
// -- added MinPolyQuotHeuristic
//
// Revision 1.34  2018/01/17 12:40:51  abbott
// Summary: Start prime in MinPolyModular now depends on SmallFpImpl::ourMaxModulus
//
// Revision 1.33  2018/01/17 10:53:00  abbott
// Summary: Many minor changes (e.g. const ring& instead of just ring)
//
// Revision 1.32  2018/01/12 18:13:52  bigatti
// -- improved coefficients
//
// Revision 1.31  2017/09/12 06:01:25  bigatti
// -- Fixed immature reconstraction in MinPolyModular
//
// Revision 1.30  2017/09/06 11:56:28  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.29  2017/07/14 19:23:17  abbott
// Summary: Changed NewRingFp into NewZZmod -- safer, more general, and just as fast
//
// Revision 1.28  2017/07/07 10:23:26  bigatti
// -- printing of polynomials raised to level 85
//
// Revision 1.27  2017/07/07 09:48:40  bigatti
// -- changed NextPrime --> PrevPrime
//
// Revision 1.26  2017/06/29 07:20:12  bigatti
// -- input check for 3rd indet
//
// Revision 1.25  2017/06/26 13:17:24  bigatti
// -- exported utility function LiftPolyFromFpxToQQx
//
// Revision 1.24  2017/05/15 10:24:37  bigatti
// -- detecting UGLY PRIME using divisibility (instead of try/catch)
//
// Revision 1.23  2017/05/11 08:46:46  bigatti
// -- added error ERR::CannotReconstruct
// -- cleaned up code accordingly
//
// Revision 1.22  2017/05/11 07:20:28  bigatti
// -- fixed bug about first prime being an ugly prime
//
// Revision 1.21  2017/05/09 13:52:23  bigatti
// -- fixed bug about uncaught ugly prime
// -- checked in buggy buggy ShapeLemma function
//
// Revision 1.20  2017/04/07 14:19:54  bigatti
// -- increased verbosity level to 80/90
//
// Revision 1.19  2017/02/16 12:41:40  bigatti
// -- highlighted code for CoefficientsInMatCol
//
// Revision 1.18  2017/02/06 16:08:34  bigatti
// -- just some verbosity
//
// Revision 1.17  2017/02/01 09:48:04  bigatti
// -- added verbosity printing
//
// Revision 1.16  2016/12/06 16:08:23  bigatti
// -- improved MinPolyModular
//
// Revision 1.15  2016/11/25 17:06:53  bigatti
// -- increased verbosity for some debugging entries
//
// Revision 1.14  2016/11/23 13:33:11  bigatti
// -- MinPolyModular called automatically by MinPolyQuot
// -- MinPolyModular with verbose and interrupt
//
// Revision 1.13  2016/10/27 13:07:24  bigatti
// -- added MinPolyQuotDefLin
// -- cleaned up code, unified notation
//
// Revision 1.12  2016/10/24 12:10:11  bigatti
// -- added ShapeLemma (first draft)
//
// Revision 1.11  2016/10/10 16:22:15  bigatti
// -- renamed myColTbl --> myColIndices
// -- aligned code for timings in MinPoly functions
// -- added check for ideal(1)
//
// Revision 1.10  2016/09/22 14:13:46  bigatti
// -- added safety checks in FrobeniusMat
//
// Revision 1.9  2016/09/08 14:13:42  bigatti
// -- using now HasUniqueOwner
//
// Revision 1.8  2016/06/27 14:50:28  bigatti
// -- now FrobeniusMat may take two args
//
// Revision 1.7  2016/06/20 15:24:42  bigatti
// -- renamed MinPolyXX --> MinPolyQuotXX
// -- added MinPolyDef(f, x)
// -- added default to MinPolyDef in both cases
//
// Revision 1.6  2016/06/14 08:26:01  bigatti
// -- added MINPOLY_DEBUG flag for enabling debugging info
//
// Revision 1.5  2016/06/10 15:50:42  bigatti
// -- new version of FrobeniusMat
//
// Revision 1.4  2016/06/06 11:45:42  bigatti
// -- added initialization for timer
//
// Revision 1.3  2016/04/14 11:35:32  bigatti
// -- added FrobeniusMat (by Elisa Palezzato)
//
// Revision 1.2  2016/03/18 15:13:47  abbott
// Summary: Commented out two unused variables (for timing)
//
// Revision 1.1  2016/03/18 11:54:00  bigatti
// -- first import
//
