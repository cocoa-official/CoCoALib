//   Copyright (c)  2005,2008,2016  John Abbott and Anna M. Bigatti

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

#include "CoCoA/MatrixOps.H"

#include "CoCoA/BigRatOps.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MatrixFp.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-prime.H"
#include "CoCoA/NumTheory-CRT.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/ring-AutomaticConversion.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/ToString.H"
#include "CoCoA/assert.H"
#include "CoCoA/bool3.H"
#include "CoCoA/combinatorics.H"
#include "CoCoA/convert.H"
#include "CoCoA/error.H"
#include "CoCoA/geobucket.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/time.H"
#include "CoCoA/utils.H" // for len
#include "CoCoA/VectorOps.H"
#include "CoCoA/verbose.H"

// BELOW IS JUST AN EXPERIMENT!  2023-11-30
//extern "C"
//{
#include "TmpFactorDir/linalg/Zdet.h" // enough?
//}

#include <algorithm>
using std::min;
#include <cmath>
using std::ldexp;
#include <iostream>
using std::ostream;
#include <limits>
using std::numeric_limits;
//#include <vector>
using std::vector;

namespace CoCoA
{

  // namespace // anonymous
  // {

  //   int signature(const vector<int>& perm)
  //   {
  //     int ans = 1;
  //     int n = len(perm);
  //     for (int i=0; i < n; ++i)
  //       for (int j=i+1; j < n; ++j)
  //         if (perm[i] > perm[j]) ans = -ans;
  //     return ans;
  //   }

  // } // end of namespace anonymous

  RingElem det(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    RingElem d(RingOf(M));
    M->myDet(d);
    return d;
  }


  // More sophisticated impl in MatrixOps-IsZeroDet.C
  // bool IsZeroDet(const ConstMatrixView& M)
  // {
  //   // Currently just a very naive impl
  //    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
  //   return IsZero(det(M));
  // }


  RingElem DetOfSmallMat(ConstMatrixView M)
  {
    const long n = NumRows(M); // we know matrix is square
    if (NumCols(M) != n)  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (n > 5)  CoCoA_THROW_ERROR1(ERR::ArgTooBig);
    switch (n)
    {
    case 0:
      return one(RingOf(M));
    case 1:
      return M(0,0);
    case 2:
      return det2x2(M);
    case 3:
      return det3x3(M);
    case 4:      
      return det4x4(M);
    case 5:
      return det5x5(M);
    default:
      CoCoA_THROW_ERROR1(ERR::ShouldNeverGetHere);
      /*just to keep compiler quiet*/ throw "NEVER GET HERE";
    }
  }


  RingElem det2x2(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (NumRows(M) != 2)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    CoCoA_ASSERT(IsCommutative(RingOf(M)));
    return M(0,0)*M(1,1) - M(0,1)*M(1,0);
  }
  

  RingElem det3x3(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (NumRows(M) != 3)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    CoCoA_ASSERT(IsCommutative(RingOf(M)));
    return M(2,2) *(M(0,0)*M(1,1) - M(0,1)*M(1,0)) 
         - M(2,1) *(M(0,0)*M(1,2) - M(0,2)*M(1,0))
         + M(2,0) *(M(0,1)*M(1,2) - M(0,2)*M(1,1));
  }
  

  RingElem det4x4(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (NumRows(M) != 4)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    CoCoA_ASSERT(IsCommutative(RingOf(M)));

    // 2x2 dets for rows 2,3
    const RingElem det_col01 = M(2,0)*M(3,1) - M(2,1)*M(3,0);
    const RingElem det_col02 = M(2,0)*M(3,2) - M(2,2)*M(3,0);
    const RingElem det_col03 = M(2,0)*M(3,3) - M(2,3)*M(3,0);
    const RingElem det_col12 = M(2,1)*M(3,2) - M(2,2)*M(3,1);
    const RingElem det_col13 = M(2,1)*M(3,3) - M(2,3)*M(3,1);
    const RingElem det_col23 = M(2,2)*M(3,3) - M(2,3)*M(3,2);

    // 3x3 dets for rows 1,2,3
    const RingElem det_col012 = M(1,0)*det_col12 - M(1,1)*det_col02 + M(1,2)*det_col01;
    const RingElem det_col013 = M(1,0)*det_col13 - M(1,1)*det_col03 + M(1,3)*det_col01;
    const RingElem det_col023 = M(1,0)*det_col23 - M(1,2)*det_col03 + M(1,3)*det_col02;
    const RingElem det_col123 = M(1,1)*det_col23 - M(1,2)*det_col13 + M(1,3)*det_col12;

    // 4x4 det
    const RingElem det = M(0,0)*det_col123 - M(0,1)*det_col023 + M(0,2)*det_col013 - M(0,3)*det_col012;
    return det;
  }


  RingElem det5x5(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (NumRows(M) != 5)  CoCoA_THROW_ERROR1(ERR::BadMatrixSize);
    CoCoA_ASSERT(IsCommutative(RingOf(M)));

    // 2x2 dets for rows 3,4
    const RingElem det_col01 = M(3,0)*M(4,1) - M(3,1)*M(4,0);
    const RingElem det_col02 = M(3,0)*M(4,2) - M(3,2)*M(4,0);
    const RingElem det_col03 = M(3,0)*M(4,3) - M(3,3)*M(4,0);
    const RingElem det_col04 = M(3,0)*M(4,4) - M(3,4)*M(4,0);
    const RingElem det_col12 = M(3,1)*M(4,2) - M(3,2)*M(4,1);
    const RingElem det_col13 = M(3,1)*M(4,3) - M(3,3)*M(4,1);
    const RingElem det_col14 = M(3,1)*M(4,4) - M(3,4)*M(4,1);
    const RingElem det_col23 = M(3,2)*M(4,3) - M(3,3)*M(4,2);
    const RingElem det_col24 = M(3,2)*M(4,4) - M(3,4)*M(4,2);
    const RingElem det_col34 = M(3,3)*M(4,4) - M(3,4)*M(4,3);

    // 3x3 dets for rows 2,3,4
    const RingElem det_col012 = M(2,0)*det_col12 - M(2,1)*det_col02 + M(2,2)*det_col01;
    const RingElem det_col013 = M(2,0)*det_col13 - M(2,1)*det_col03 + M(2,3)*det_col01;
    const RingElem det_col014 = M(2,0)*det_col14 - M(2,1)*det_col04 + M(2,4)*det_col01;
    const RingElem det_col023 = M(2,0)*det_col23 - M(2,2)*det_col03 + M(2,3)*det_col02;
    const RingElem det_col024 = M(2,0)*det_col24 - M(2,2)*det_col04 + M(2,4)*det_col02;
    const RingElem det_col034 = M(2,0)*det_col34 - M(2,3)*det_col04 + M(2,4)*det_col03;
    const RingElem det_col123 = M(2,1)*det_col23 - M(2,2)*det_col13 + M(2,3)*det_col12;
    const RingElem det_col124 = M(2,1)*det_col24 - M(2,2)*det_col14 + M(2,4)*det_col12;
    const RingElem det_col134 = M(2,1)*det_col34 - M(2,3)*det_col14 + M(2,4)*det_col13;
    const RingElem det_col234 = M(2,2)*det_col34 - M(2,3)*det_col24 + M(2,4)*det_col23;

    // 4x4 dets for rows 1,2,3,4
    const RingElem det_col0123 = M(1,0)*det_col123 - M(1,1)*det_col023 + M(1,2)*det_col013 - M(1,3)*det_col012;
    const RingElem det_col0124 = M(1,0)*det_col124 - M(1,1)*det_col024 + M(1,2)*det_col014 - M(1,4)*det_col012;
    const RingElem det_col0134 = M(1,0)*det_col134 - M(1,1)*det_col034 + M(1,3)*det_col014 - M(1,4)*det_col013;
    const RingElem det_col0234 = M(1,0)*det_col234 - M(1,2)*det_col034 + M(1,3)*det_col024 - M(1,4)*det_col023;
    const RingElem det_col1234 = M(1,1)*det_col234 - M(1,2)*det_col134 + M(1,3)*det_col124 - M(1,4)*det_col123;

    // 5x5 det
    const RingElem det = M(0,0)*det_col1234 - M(0,1)*det_col0234 + M(0,2)*det_col0134 - M(0,3)*det_col0124 + M(0,4)*det_col0123;
    
    return det;
  }



  // Compute det as "alternating" sum of products

  namespace // anonymous
  {

    RingElem DetDirect_generic(ConstMatrixView M)
    {
      CoCoA_ASSERT(NumRows(M) == NumCols(M));
      const long n = NumRows(M);
      const ring R = RingOf(M);

      vector<int> cols(n);
      for (long i=0; i < n; ++i) { cols[i] = i; }

      RingElem det(R);
      do
      {
        CheckForInterrupt("DetDirect_generic");
        RingElem term = signature(cols)*one(R);
        for (int i=0; i < n; ++i)
          term *= M(i, cols[i]);
        det += term;
      } while (std::next_permutation(cols.begin(),cols.end()) );
      return det;
    }

    
    // Supposedly faster version using geobuckets -- BUG: messy & undocumented :-(
    RingElem DetDirect_SparsePoly(ConstMatrixView M)
    {
      CoCoA_ASSERT(NumRows(M) == NumCols(M));
      const long n = NumRows(M);
      const ring R = RingOf(M);

      vector<int> cols(n);
      for (long i=0; i < n; ++i) { cols[i] = i; }

      RingElem det(R);
      geobucket ans(R);
      do
      {
        RingElem term = one(R);
        for (int i=0; i < n; ++i)
          term *= M(i, cols[i]);
//      det += signature(cols)*term;
        if (IsZero(term)) continue;
        if (NumTerms(term) == 1)
          ans.myAddMulLM(term, RingElem(R,signature(cols)), 1);
        else det += signature(cols)*term;
      } while ( std::next_permutation(cols.begin(),cols.end()) );

      AddClear(det,ans);
      return det;
    }

  } // end of namespace anon


  RingElem DetDirect(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (IsSparsePolyRing(RingOf(M))) return DetDirect_SparsePoly(M);
    return DetDirect_generic(M);
  }


  // Known defect: this algm is valid only over IntegralDomains
  RingElem DetByGauss(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    const long N = NumRows(M);
    const ring R(RingOf(M));
    // ??? this function works only within integral domains!
    if (!IsIntegralDomain(R))  CoCoA_THROW_ERROR1(ERR::NotIntegralDomain);
// Is it worth adding code for this special case?  General code works fine.
//     // Handle 0x0 matrix specially: its det is 1.
//     if (N == 0) { d = one(R); return; }
    const ring K((IsField(R) ? R : NewFractionField(R)));
    matrix Gss(CanonicalHom(R,K)( M ));
    RingElem c(K);
    RingElem determinant(one(K));
    for (long col=0; col < N; ++col)
    {
      // Pick a good pivot row
      long PivotRow = -1;
      for (long r=col ; r < N; ++r)
      {
        if (!IsZero(Gss(r,col))) { PivotRow = r; break; }
      }

      if (PivotRow == -1)
      {
        return zero(R);
      }
      if (PivotRow != col)
      {
        Gss->mySwapRows(PivotRow,col);
        determinant *= -1;
      }
      c = Gss(col,col);
      determinant *= c;
      for (long i=col+1; i < N; ++i)
      {
        CheckForInterrupt("DetByGauss(Mat)");
        Gss->myAddRowMul(i, col, -Gss(i,col)/c);
      }
    }
    if (R==K) return determinant;
    else return num(determinant)/den(determinant);
  }


  // Simple (but probably not very fast).
  RingElem DetByMinors(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    CheckForInterrupt("DetByMinors");
    const long n = NumRows(M);
    if (n < 6) return DetOfSmallMat(M);
    const ring& R = RingOf(M);
//    if (n == 0) return one(R);
//    if (n == 1) return M(0,0);

    // SLUG: should choose a good row/col along which to expand!
    // Naive impl: expand along first row. -- should try to pick best row/col
    vector<long> rows(n-1);
    vector<long> cols(n-1);
    const long ExpandRow = 0;
    for (long i=0; i < n-1; ++i) { rows[i] = i+1; }
    for (long j=0; j < n-1; ++j) { cols[j] = j+1; }
    RingElem D = zero(R);
    int sign = 1; // depends on which row/col we expand along!!!!!
    for (long j=0; j < n; ++j)
    {
      if (!IsZero(M(ExpandRow,j)))
        D += sign*M(ExpandRow,j)*DetByMinors(submat(M,rows,cols));
      sign = -sign;
      if (j != n-1) cols[j] = j;
    }
    return D;
  }


  
  struct HadamardRowCol_BigInt
  {
    HadamardRowCol_BigInt(const BigInt& R, const BigInt& C): myRowBoundSq(R), myColBoundSq(C) {}
    BigInt myRowBoundSq;
    BigInt myColBoundSq;
  };

  // SQUARE of min of row and col Hadamard bounds
  HadamardRowCol_BigInt HadamardBoundSq(const std::vector<std::vector<BigInt>>& M)
  {
    const int n = len(M);
    // ASSUMES vec-vec is "square"
    BigInt RowBoundSq(1);
    for (int i=0; i < n; ++i)
    {
      CheckForInterrupt("HadamardBoundSq");
      BigInt L2;
      for (int j=0; j < n; ++j)
        L2 += power(M[i][j],2);
      RowBoundSq *= L2;
    }

    BigInt ColBoundSq(1);
    for (int j=0; j < n; ++j)
    {
      CheckForInterrupt("HadamardBoundSq");
      BigInt L2;
      for (int i=0; i < n; ++i)
        L2 += power(M[i][j],2);
      ColBoundSq *= L2;
    }

    return HadamardRowCol_BigInt(RowBoundSq, ColBoundSq);
  }


  HadamardRowCol HadamardBoundSq(const ConstMatrixView& M)
  {
//    CoCoA_ASSERT(IsZero(characteristic(RingOf(M))));  ???
    const int n = NumRows(M);
    if (NumCols(M) != n)  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    const ring& R = RingOf(M);
    if (characteristic(R) == 0)
    {
      try
      {
        vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
        for (int i=0; i < n; ++i)
          for (int j=0; j < n; ++j)
            Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));
        const HadamardRowCol_BigInt bounds = HadamardBoundSq(Mcopy);
        return HadamardRowCol(RingElem(R, bounds.myRowBoundSq), RingElem(R, bounds.myColBoundSq));
      }
      catch (...) {} // discard ERR if conversion to BigInt failed
    }
    
    RingElem RowBoundSq = one(R);
    for (int i=0; i < n; ++i)
    {
      CheckForInterrupt("HadamardBoundSq");
      RingElem L2 = zero(R);;
      for (int j=0; j < n; ++j)
        L2 += power(M(i,j),2);
      RowBoundSq *= L2;
    }

    RingElem ColBoundSq = one(R);
    for (int j=0; j < n; ++j)
    {
      CheckForInterrupt("HadamardBoundSq");
      RingElem L2 = zero(R);
      for (int i=0; i < n; ++i)
        L2 += power(M(i,j),2);
      ColBoundSq *= L2;
    }

    return HadamardRowCol(RowBoundSq, ColBoundSq);
  }




  long HadamardRowScale(const ConstMatrixView& M); //fwd decl
  long HadamardColScale(const ConstMatrixView& M); //fwd decl
  BigInt DetBound_GPS(const ConstMatrixView& M);


  RingElem DetByCRT(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("DetByCRT");
    const ring& ZZ = RingOf(M);
    CoCoA_ASSERT(IsZZ(ZZ));
    const int n = NumRows(M);
    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));

    const HadamardRowCol_BigInt H = HadamardBoundSq(Mcopy); // square of hadamard bounds (both row and col bounds)
    const BigInt CRTBound = 2*FloorSqrt(min(H.myRowBoundSq, H.myColBoundSq)); // factor of 2 to allow for +/- sign
    VERBOSE(80) << "Hrow(sq) = " << FloatStr(H.myRowBoundSq) << "  Hcol(sq) = " << FloatStr(H.myColBoundSq) << "   GPS = " << FloatStr(DetBound_GPS(M)) << std::endl;

 const long RowScale = HadamardRowScale(M);
 const long ColScale = HadamardColScale(M);
 VERBOSE(80) << "RowScale = " << RowScale << "   ColScale = " << ColScale << std::endl;
 VERBOSE(80) << "Scaled Hrow = " << FloatStr(FloorSqrt(H.myRowBoundSq/power(2,RowScale))) <<
   "   Scaled Hcol = " << FloatStr(FloorSqrt(H.myColBoundSq/power(2,ColScale))) << std::endl;

 VERBOSE(80) << "Actually using CRTBound = " << FloatStr(CRTBound) << std::endl;
    int NumPrimes = 0;

    bool UsingHeuristic = false;
    PrimeSeqForCRT PSeq;  // this produces the primes we shall use for CRTing
    CRTMill CRT;
    while (true)
    {
      const bool finished = (CombinedModulus(CRT) > CRTBound);
      const long log2modulus = FloorLog2(CombinedModulus(CRT));
      if (finished || (UsingHeuristic && log2modulus > 50 && (IsZero(CombinedResidue(CRT)) || FloorLog2(CombinedResidue(CRT))+50 < log2modulus)))
      {
        VERBOSE(80) << "Returning; num iters: " << NumPrimes << std::endl;
        return RingElem(ZZ, CombinedResidue(CRT));
      }

      CheckForInterrupt("DetByCRT");
//      p = NextPrime(p);
      const SmallPrime p = NextPrime(PSeq);
      ++NumPrimes;
      VERBOSE(85) << "Iter " << NumPrimes << "  Using prime p=" << p << std::endl;
      SmallFpImpl ModP(p);
      MatrixFp Mp(ModP, n,n);
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          Mp(i,j) = ModP.myReduce(Mcopy[i][j]);

      const long dp = det(Mp);
      CRT.myAddInfo(dp,p, CRTMill::CoprimeModulus);
    }
  }


    RingElem Det_OLD(const ConstMatrixView& M)
    {
    VerboseLog VERBOSE("Det_OLD");
    const ring& ZZ = RingOf(M);
    CoCoA_ASSERT(IsZZ(ZZ));
    const int n = NumRows(M);
    Zmat Mcopy = Zmat_ctor(n,n);
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        mpz_set(Mcopy->entry[i][j], mpzref(ConvertTo<BigInt>(M(i,j))));  // hopefully safe
    BigInt D;
    Zdet(mpzref(D), Mcopy);
    Zmat_dtor(Mcopy);
    return RingElem(ZZ,D);
    }


  struct IPentry
  {
    IPentry(double cos2, long i1, long i2): cosine2(cos2), row1(i1), row2(i2) {}
    double cosine2;
    long row1;
    long row2;
  };

  inline bool ByCosine2(const IPentry& A, const IPentry& B)
  {
    return A.cosine2 > B.cosine2;
  }

  // Returns an integer K such that abs val of true det is less than
  // row Hadamard bound times 2^(-K).
  // Complexity is CUBIC in matrix dimension!!   ASSUMES ring is ZZ
  long HadamardRowScale(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("HadamardRowScale");
//    const double t0 = CpuTime();
    CoCoA_ASSERT(IsZZ(RingOf(M)));
    const int n = NumRows(M);

    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));
//    const double t1 = CpuTime();
//    VERBOSE(81) << "Copy time: " << t1-t0 << std::endl;

    vector<BigInt> L2(n);
    for (int i=0; i < n; ++i)
    {
      BigInt L2row;
      for (int j=0; j < n; ++j)
        L2row += power(Mcopy[i][j],2);
      L2[i] = L2row;
    }
//    const double t2 = CpuTime();
//    VERBOSE(81) << "L2 time: " << t2-t1 << std::endl;

    vector< vector<double> > Mdbl(n, vector<double>(n));
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
      {
        long mij_exp;
        long len_exp;
        double mij_mant = mpz_get_d_2exp(&mij_exp, mpzref(Mcopy[i][j]));
        double len_mant = mpz_get_d_2exp(&len_exp, mpzref(L2[i]));
        if (IsOdd(len_exp)) { ++len_exp; len_mant /= 2.0; }
        len_mant = sqrt(len_mant); len_exp /= 2; // exact division!
        if (len_exp - mij_exp > 50) { Mdbl[i][j] = 0.0; continue; }
        const double shift = ldexp(1.0, mij_exp - len_exp);
        Mdbl[i][j] = (mij_mant/len_mant)*shift;
//        VERBOSE(89) << "Entry (" << i << "," << j << ") = " << Mdbl[i][j] << std::endl;
      }
//    const double t3 = CpuTime();
//    VERBOSE(81) << "Mdbl time: " << t3-t2 << std::endl;

    vector<IPentry> tbl; tbl.reserve((n*n-n)/2);
    for (int i1=0; i1 < n; ++i1)
    {
      double BiggestInnerProd = 0.0;
      for (int i2=i1+1; i2 < n; ++i2)
      {
        double InnerProd = 0.0;
        for (int j=0; j < n; ++j)
          InnerProd += Mdbl[i1][j]*Mdbl[i2][j];
        tbl.push_back(IPentry(pow(InnerProd,2),i1,i2));
        if (std::abs(InnerProd) > BiggestInnerProd) BiggestInnerProd = std::abs(InnerProd);
      }
      if (i1 == 0)
      {
//        VERBOSE(81) << "BiggestInnerProd = " << BiggestInnerProd << std::endl;
        
        if (BiggestInnerProd < 0.2) { /*VERBOSE(80) << "Not worth it" << std::endl;*/ return 0;} // heuristic - not worth computing redn factor
      }
    }
//    const double t4 = CpuTime();
//    VERBOSE(81) << "IP tbl time: " << t4-t3 << std::endl;
    std::sort(tbl.begin(), tbl.end(), ByCosine2);
//    const double t5 = CpuTime();
//    VERBOSE(81) << "sort time: " << t5-t4 << std::endl;
    // If all pairs are practically orthog, just return 0 = log(1).
    if (tbl[0].cosine2 < 1.0/n) { /*VERBOSE(80) << "Practically orthog" << std::endl;*/ return 0; }
    
    vector<int> TreeIndex(n); for (int i=0; i < n; ++i) TreeIndex[i] = i;
    int NumTrees = n;
    double RednFactor = 1.0;
    constexpr double eps = 1.0/(1048576.0*1048576.0);
    const double Shift256 = pow(2.0, 256);
    const double LowerLimit = 1.0/Shift256;
    long exp2 = 0;
    int NumBigJumps = 0;
    for (int i=0; i < len(tbl); ++i)
    {
//      if (i < 10) VERBOSE(81) << "cos2 = " << tbl[i].cosine2 << "   (i1,i2)=("<<tbl[i].row1<<","<<tbl[i].row2<<")"<<std:: endl;
      if (tbl[i].cosine2 < 0.01) break;
      const int index1 =TreeIndex[tbl[i].row1];
      const int index2 =TreeIndex[tbl[i].row2];
      if (index1 == index2) continue;
      if (index1 < index2) TreeIndex[tbl[i].row2] = index1;
      else TreeIndex[tbl[i].row1] = index2;
      if (tbl[i].cosine2 < 1.0-eps)
      RednFactor *= (1-tbl[i].cosine2);
      else { RednFactor *= eps; ++NumBigJumps; }
      if (RednFactor < LowerLimit) { RednFactor *= Shift256; exp2 += 256; }
      --NumTrees;
      if (NumTrees == 1) break;
    }
    int exp2rest = 0;
    if (RednFactor != 1) frexp(RednFactor, &exp2rest);
    VERBOSE(81) << "NumBigJumps = " << NumBigJumps << std::endl;
//    VERBOSE(81) << "RednFactor = " << RednFactor << "  exp2rest = " << exp2rest << std::endl;
//    VERBOSE(81) << "Returning " << exp2-exp2rest << std::endl;
//    VERBOSE(81) << "TIME: " << CpuTime()-t0 << std::endl;
    return exp2-exp2rest;  // subtract 1??
  }


  // Lazy man's impl -- at least it's "obviously correct" (right?)
  long HadamardColScale(const ConstMatrixView& M)
  {
    return HadamardRowScale(transpose(M));
  }

  
  // ASSUMES ring is ZZ (or QQ?)
  // SHOULD PREPROCESS MATRIX:
  // (1) make sure that all row sums and all col sums are >= 0
  //     (by changing signs of a row/col --> does not affect abs(det))
  // (2) rescale rows/cols so that most entries are about the same size
  //     (changes det by a known factor).  JAA is not sure how to do this well.
  BigInt DetBound_GPS(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("DetBound_GPS");
    CoCoA_ASSERT(IsZZ(RingOf(M)));
    const int n = NumRows(M);

    vector< vector<BigInt> > Mcopy(n, vector<BigInt>(n, BigInt(0)));
    vector<BigInt> RowSum(n);
    vector<BigInt> ColSum(n);
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
      {
        Mcopy[i][j] = ConvertTo<BigInt>(M(i,j));
        RowSum[i] += Mcopy[i][j];
        ColSum[j] += Mcopy[i][j];
      }

    BigInt SumAll;
    for (int i=0; i < n; ++i) SumAll += RowSum[i];
    VERBOSE(150) << "(INPUT MAT) Entry sum = " << FloatStr(SumAll) << std:: endl;

    // Make sure all RowSums and ColSums are >= 0
    while (true)
    {
      int MinRow = 0; // index of min value
      for (int i=1; i < n; ++i)
        if (RowSum[i] < RowSum[MinRow]) MinRow = i;
      int MinCol = 0; // index of min value
      for (int j=1; j < n; ++j)
        if (ColSum[j] < ColSum[MinCol]) MinCol = j;
      if (RowSum[MinRow] >= 0 && ColSum[MinCol] >= 0) break;
      if (RowSum[MinRow] <= ColSum[MinCol])
      {
        // Negate MinRow
        RowSum[MinRow] = -RowSum[MinRow];
        for (int j=0; j < n; ++j)
        {
          Mcopy[MinRow][j] = -Mcopy[MinRow][j];
          ColSum[j] += 2*Mcopy[MinRow][j];
        }
      }
      else
      {
        // negate MinCol
        ColSum[MinCol] = -ColSum[MinCol];
        for (int i=0; i < n; ++i)
        {
          Mcopy[i][MinCol] = -Mcopy[i][MinCol];
          RowSum[i] += 2*Mcopy[i][MinCol];
        }
      }
    }

    VERBOSE(150) << "After first stage" << std::endl;
    SumAll = 0;
    for (int i=0; i < n; ++i) SumAll += RowSum[i];
    VERBOSE(150) << "Entry sum = " << FloatStr(SumAll) << std:: endl;

    // Now check for row-col pairs
    bool changed = true;
    while (changed)
    {
      changed = false;
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
        {
          if (RowSum[i]+ColSum[j] >= 2*Mcopy[i][j]) continue;
          changed = true;
//          VERBOSE(22) << "Pivot " << i << "," << j << std::endl;
          // negate col j
          ColSum[j] = -ColSum[j];
          for (int ii=0; ii < n; ++ii)
          {
            Mcopy[ii][j] = -Mcopy[ii][j];
            RowSum[ii] += 2*Mcopy[ii][j];
          }
          // negate row i
          RowSum[i] = -RowSum[i];
          for (int jj=0; jj < n; ++jj)
          {
            Mcopy[i][jj] = -Mcopy[i][jj];
            ColSum[jj] += 2*Mcopy[i][jj];
          }
        }
    }
    
    VERBOSE(150) << "After second stage" << std::endl;
    SumAll = 0;
    for (int i=0; i < n; ++i) SumAll += RowSum[i];
    VERBOSE(150) << "Entry sum = " << FloatStr(SumAll) << std:: endl;

    BigInt SumEntries;
    BigInt SumSquares;
    for (int i=0; i < n; ++i)
      for (int j=0; j < n; ++j)
      {
        SumEntries += Mcopy[i][j];
        SumSquares += power(Mcopy[i][j],2);
      }

    const BigInt delta = power(SumEntries,2) - n*SumSquares;
    VERBOSE(150) << "alpha = " << SumEntries << "  over " << n << std::endl;
    VERBOSE(150) << "beta = " << SumSquares << "  over " << n << std::endl;
    VERBOSE(150) << "delta = " << delta << "  over " << n*n*(n-1) << std::endl;
    if (delta <= 0)
    {
      // Always worse than Hadanard's bound
      if (IsEven(n)) return power(SumSquares,n/2);
      return FloorSqrt(power(SumSquares,n));
    }
    VERBOSE(150) << "Good case" << std::endl;
    if (IsOdd(n))
      return (abs(SumEntries)*power((n*(n-1))*SumSquares-delta, (n-1)/2))/(power(n,n)*power((n-1),(n-1)/2)); // int division
    return FloorSqrt(power(SumEntries,2)*power((n*(n-1))*SumSquares-delta,n-1)/power(n-1,n-1))/power(n,n); // integer division (twice)
  }



//////////////////////////////////////////////////////////////////
// Bareiss

  namespace /*anonymous*/
  {
    RingElem DetByBareiss_ZZ(const ConstMatrixView& M)
    {
      CoCoA_ASSERT(IsZZ(RingOf(M)));
      const int n = NumRows(M);
      CoCoA_ASSERT(NumCols(M) == n);
      const ring R = RingOf(M);
      if (n == 0) return one(R);
      if (n == 1) return M(0,0);
      vector< vector< BigInt > > M2(n, vector<BigInt>(n));
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          IsInteger(M2[i][j], M(i,j)); // ignore return value (must be true)
      BigInt d(1);
      int sign = 1;
      for (int k=0; k < n-1; ++k)
      {
        // Find a non-zero pivot in k-th column
        int row = -1;
        for (int i=k; i < n; ++i)
          if (!IsZero(M2[i][k])) { row = i; break; }
        if (row == -1) return zero(R);
        if (row != k) { swap(M2[k], M2[row]); sign = -sign; }
        BigInt tmp; // temporary workspace used in inner loop below
        for (int i=k+1; i < n; ++i)
          for (int j=k+1; j < n; ++j)
          {
//  This block effectively does the following, but is usefully faster (about 4x)
//            M2[i][j] = (M2[i][j]*M2[k][k]-M2[i][k]*M2[k][j])/d;
            mpz_mul(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(M2[k][k]));
            mpz_mul(mpzref(tmp), mpzref(M2[i][k]), mpzref(M2[k][j]));
            mpz_sub(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(tmp));
            mpz_divexact(mpzref(M2[i][j]), mpzref(M2[i][j]), mpzref(d));
          }
        d = M2[k][k];
      }
      return RingElem(R,sign*M2[n-1][n-1]);
    }
  } // end of anonymous namespace


  RingElem DetByBareiss(const ConstMatrixView& M)
  {
    const ring R = RingOf(M);
    if (IsZZ(R)) return DetByBareiss_ZZ(M); // a bit faster
    CoCoA_ASSERT(IsIntegralDomain(R));
    const int n = NumRows(M);
    CoCoA_ASSERT(NumCols(M) == n);
    if (n == 0) return one(R);
    if (n == 1) return M(0,0);
    matrix M2 = NewDenseMat(M);
    RingElem d = one(R);
    int sign = 1;
    for (int k=0; k < n-1; ++k)
    {
      // Find a non-zero pivot in k-th column
      int row = -1;
      for (int i=k; i < n; ++i)
        if (!IsZero(M2(i,k))) { row = i; break; }
      if (row == -1) return zero(R);
      if (row != k) { M2->mySwapRows(k, row); sign = -sign; }
      // Now use pivot row to reduce all lower rows
      for (int i=k+1; i < n; ++i)
        for (int j=k+1; j < n; ++j)
        {
          CheckForInterrupt("DetByBareiss");
          SetEntry(M2,i,j,(M2(i,j)*M2(k,k)-M2(i,k)*M2(k,j))/d );
        }
      d = M2(k,k);
    }
    if (sign == 1)
      return M2(n-1,n-1);
    return -M2(n-1,n-1);
  }


  

  namespace // anonymous
  {

    BigInt FindRowAndColScales(vector<BigInt>& Rscale, vector<BigInt>& Cscale, const ConstMatrixView& M)
    {
      VerboseLog VERBOSE("FindRowAndColScales");

      CoCoA_ASSERT(IsQQ(RingOf(M)));
      CoCoA_ASSERT(NumCols(M) == NumRows(M));
      const int n = NumRows(M);
      vector< vector<BigInt> > denom(n, vector<BigInt>(n));
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
          denom[i][j] = ConvertTo<BigInt>(den(M(i,j)));

//    VERBOSE(25) << "denom: " << denom << std::endl;
      vector<BigInt> RowScale(n, BigInt(1));
      vector<BigInt> ColScale(n, BigInt(1));

      // Spit notionally into k=sqrt(n) blocks;
      // compute gcd of each block, then take largest lcm of pairs of blocks
      const int sqrtn = 1+FloorSqrt(n);
      vector<BigInt> BlockGCD(sqrtn);
      // Compute preliminary row factors
      for (int i=0; i < n; ++i) // row indexer
      {
        int end=n;
        for (int k=0; k < sqrtn; ++k)
        {
          BigInt g;
          const int BlockSize = end/(sqrtn-k); // integer division!
          for (int j=end-BlockSize; j < end; ++j)
            g = gcd(g, denom[i][j]);
          BlockGCD[k] = g;
          end -= BlockSize;
        }
//      VERBOSE(25) << "BlockGCD: " << BlockGCD << std::endl;
        BigInt scale = BlockGCD[0];
        for (int k=1; k < sqrtn; ++k)
          scale = lcm(scale, BlockGCD[k]);
        RowScale[i] = scale;
        if (scale != 1)
          for (int j=0; j < n; ++j)
            denom[i][j] /= gcd(scale, denom[i][j]);
      }
      if (VerbosityLevel() >= 90)
      {
        BigInt Stage1Fac(1);
        for (int i=0; i < n; ++i)
          Stage1Fac *= RowScale[i];
        VERBOSE(90) << "Stage1 RowFac = " << FloatStr(Stage1Fac) << std::endl;
      }
      if (VerbosityLevel() >= 91)
      {
        vector<long> LogRowScale(n);
        for (int i=0; i < n; ++i)
          LogRowScale[i] = FloorLog2(RowScale[i]);
        VERBOSE(91) << "Stage 1 LogRowScale: " << LogRowScale << std::endl;
      }

      // Compute preliminary col factors
      for (int j=0; j < n; ++j) // col indexer
      {
        int end=n;
        for (int k=0; k < sqrtn; ++k)
        {
          BigInt g;
          const int BlockSize = end/(sqrtn-k); // integer division!
          for (int i=end-BlockSize; i < end; ++i)
            g = gcd(g, denom[i][j]);
          BlockGCD[k] = g;
          end -= BlockSize;
        }
//      VERBOSE(25) << "BlockGCD: " << BlockGCD << std::endl;
        BigInt scale = BlockGCD[0];
        for (int k=1; k < sqrtn; ++k)
          scale = lcm(scale, BlockGCD[k]);
        ColScale[j] = scale;
        if (scale != 1)
          for (int i=0; i < n; ++i)
            denom[i][j] /= gcd(scale, denom[i][j]);
      }
      if (VerbosityLevel() >= 90)
      {
        BigInt Stage1Fac(1);
        for (int j=0; j < n; ++j)
          Stage1Fac *= ColScale[j];
        VERBOSE(90) << "Stage1 ColFac = " << FloatStr(Stage1Fac) << std::endl;
      }
      if (VerbosityLevel() >= 91)
      {
        vector<long> LogColScale(n);
        for (int j=0; j < n; ++j)
          LogColScale[j] = FloorLog2(ColScale[j]);
        VERBOSE(91) << "Stage 1 LogColScale: " << LogColScale << std::endl;
      }

//    VERBOSE(25) << "denom: " << denom << std::endl;
//    VERBOSE(25) << "RowScale: " << RowScale << std::endl;

#if 0
      // First idea: gcd(lcm(first_half), lcm(second_half))
      // Works OK, but not great
      const int n2 = n/2; // integer division!
      // Fill RowScale -- phase 1
      for (int i=0; i < n; ++i)
      {
        BigInt block1(1);
        for (int j=0; j < n2; ++j)
          block1 = lcm(block1, denom[i][j]);
        BigInt block2(1);
        for (int j=n2; j < n; ++j)
          block2 = lcm(block2, denom[i][j]);
        const BigInt scale = gcd(block1, block2);
        RowScale[i] = scale;
        if (scale != 1)
          for (int j=0; j < n; ++j)
            denom[i][j] /= gcd(scale, denom[i][j]);
      }
    
      if (VerbosityLevel() >= 90)
      {
        vector<long> LogRowScale(n);
        for (int i=0; i < n; ++i)
          LogRowScale[i] = FloorLog2(RowScale[i]);
        VERBOSE(90) << "Stage 1 LogRowScale: " << LogRowScale << std::endl;
      }

      // Fill ColScale -- phase 1
      for (int j=0; j < n; ++j)
      {
        BigInt block1(1);
        for (int i=0; i < n2; ++i)
          block1 = lcm(block1, denom[i][j]);
        BigInt block2(1);
        for (int i=n2; i < n; ++i)
          block2 = lcm(block2, denom[i][j]);
        const BigInt scale = gcd(block1, block2);
        ColScale[j] = scale;
        if (scale != 1)
          for (int i=0; i < n; ++i)
            denom[i][j] /=  gcd(scale, denom[i][j]);
      }

      if (VerbosityLevel() >= 90)
      {
        vector<long> LogColScale(n);
        for (int i=0; i < n; ++i)
          LogColScale[i] = FloorLog2(ColScale[i]);
        VERBOSE(90) << "Stage 1 LogColScale: " << LogColScale << std::endl;
      }
#endif

// #if 0
//       THESE TWO LOOPS ARE PROBABLY NOT A GOOD IDEA!!!
//         // First reduce by rows, then by cols
//         // Uh oh: what is the complexity of this loop???
//         for (int i=0; i < n; ++i)
//       {
//         for (int j1=0; j1 < n; ++j1)
//         {
//           if (IsOne(denom[i][j1])) continue;
//           for (int j2=j1+1; j2 < n; ++j2)
//           {
//             const BigInt g = gcd(denom[i][j1], denom[i][j2]);
//             if (IsOne(g)) continue;
//             RowScale[i] *= g;
//             for (int j3=0; j3 < n; ++j3)
//               denom[i][j3] = denom[i][j3]/gcd(g, denom[i][j3]);
//           }
//         }
//       }
//       VERBOSE(30) << "RowScale: " << RowScale << std::endl;

//       // Now columns
//       for (int j=0; j < n; ++j)
//       {
//         for (int i1=0; i1 < n; ++i1)
//         {
//           if (IsOne(denom[i1][j])) continue;
//           for (int i2=i1+1; i2 < n; ++i2)
//           {
//             const BigInt g = gcd(denom[i1][j], denom[i2][j]);
//             if (IsOne(g)) continue;
//             ColScale[j] *= g;
//             for (int i3=0; i3 < n; ++i3)
//               denom[i3][j] = denom[i3][j]/gcd(g, denom[i3][j]);
//           }
//         }
//       }
//       VERBOSE(30) << "ColScale: " << ColScale << std::endl;
//       END OF NOT A GOOD IDEA
// #endif
    
      // Compute row-by-row LCMs
      vector<BigInt> LCMrow(n);
      BigInt ByRowsFac(1);
      for (int i=0; i < n; ++i)
      {
        BigInt RowFac(1);
        for (int j=0; j < n; ++j)
          RowFac = lcm(RowFac, denom[i][j]);
        LCMrow[i] = RowFac;
        ByRowsFac *= RowFac;
      }

      // Compute col-by-col LCMs
      vector<BigInt> LCMcol(n);
      BigInt ByColsFac(1);
      for (int j=0; j < n; ++j)
      {
        BigInt ColFac(1);
        for (int i=0; i < n; ++i)
          ColFac = lcm(ColFac, denom[i][j]);
        LCMcol[j] = ColFac;
        ByColsFac *= ColFac;
      }
      VERBOSE(85) << "ByRowsFac = " << FloatStr(ByRowsFac) << "   ByColsFac = " << FloatStr(ByColsFac) << std::endl;

      // Now pick whichever is smaller from row-by-row or col-by-col
      if (ByRowsFac <= ByColsFac)
      {
        for (int i=0; i < n; ++i)
          RowScale[i] *= LCMrow[i];
      }
      else
      {
        for (int j=0; j < n; ++j)
          ColScale[j] *= LCMcol[j];
      }


      BigInt OverallFactor(1);
      for (int i=0; i <n; ++i) OverallFactor *= RowScale[i];
      for (int j=0; j <n; ++j) OverallFactor *= ColScale[j];

      if (VerbosityLevel() >= 90)
      {
        vector<long> LogRowScale(n);
        for (int i=0; i < n; ++i)
          LogRowScale[i] = FloorLog2(RowScale[i]);
        VERBOSE(90) << "Final LogRowScale: " << LogRowScale << std::endl;
        vector<long> LogColScale(n);
        for (int j=0; j < n; ++j)
          LogColScale[j] = FloorLog2(ColScale[j]);
        VERBOSE(90) << "Final LogColScale: " << LogColScale << std::endl;
      }
      VERBOSE(80) << "Overall factor " << FloatStr(OverallFactor) << std::endl;

      swap(Rscale, RowScale);  // really assignment
      swap(Cscale, ColScale);  // really assignment
      return OverallFactor;
    }


    matrix ScaleRowsAndCols(const ConstMatrixView& M, const vector<BigInt>& RowScale, const vector<BigInt>& ColScale)
    {
      const int n = NumRows(M);
      matrix ans = NewDenseMat(RingZZ(), n,n);
      for (int i=0; i < n; ++i)
        for (int j=0; j < n; ++j)
        {
          const RingElem entry = (RowScale[i]*ColScale[j])*M(i,j);
          CoCoA_ASSERT(IsOne(den(entry)));
          SetEntry(ans,i,j, num(entry));
        }

      return ans;
    }

  } // end of namespace anonymous
  
  RingElem DetOverQQ(const ConstMatrixView& M)
  {
    VerboseLog VERBOSE("DetOverQQ");
    CoCoA_ASSERT(IsQQ(RingOf(M)));
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const int n = NumRows(M);
    if (n == 0) return one(RingOf(M));
    if (n == 1) return M(0,0);
    VERBOSE(80) << "UNTRANSPOSED" << std::endl;
    vector<BigInt> Rscale;
    vector<BigInt> Cscale;
    const BigInt OverallFactor = FindRowAndColScales(Rscale, Cscale, M);
    VERBOSE(80) << std::endl;
    VERBOSE(80) << "TRANSPOSED" << std::endl;
    vector<BigInt> trRscale;
    vector<BigInt> trCscale;
    const BigInt trOverallFactor = FindRowAndColScales(trRscale, trCscale, transpose(M));

    if (OverallFactor <= trOverallFactor)
    {
      return RingElem(RingOf(M), BigRat(ConvertTo<BigInt>(det(ScaleRowsAndCols(M, Rscale, Cscale))), OverallFactor));
    }
    else
      return RingElem(RingOf(M), BigRat(ConvertTo<BigInt>(det(ScaleRowsAndCols(M, trCscale, trRscale))), trOverallFactor));
      
    // BigInt DetZZ = ConvertTo<BigInt>(det(MoverZZ));
    // BigInt CombinedDenom(1);
    // for (int i=0; i < n; ++i)
    //   CombinedDenom *= Rscale[i]*Cscale[i];
    // return RingElem(RingOf(M), BigRat(DetZZ, CombinedDenom));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-det.C,v 1.13 2024/08/28 15:23:58 bigatti Exp $
// $Log: MatrixOps-det.C,v $
// Revision 1.13  2024/08/28 15:23:58  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2
//
// Revision 1.12  2024/08/21 20:27:05  abbott
// Summary: Made HadamardBoundSq interruptible
//
// Revision 1.11  2024/07/15 17:03:56  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.10  2024/03/08 19:59:59  abbott
// Summary: Used constexpr (redmine 1511)
//
// Revision 1.9  2023/12/01 20:21:47  abbott
// Summary: Renamed fields in record produces by HadamardBoundSq
//
// Revision 1.8  2022/12/05 21:44:29  abbott
// Summary: Commented out old impl of IsZeroDet
//
// Revision 1.7  2022/05/14 16:00:00  abbott
// Summary: Corrected comment
//
// Revision 1.6  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/10/04 08:56:01  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.4  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.3  2020/10/06 19:15:30  abbott
// Summary: Added CheckForInterrupt in DetDirect_generic
//
// Revision 1.2  2020/09/28 11:15:58  abbott
// Summary: Added DetByMinors; improved DetDirect so it works in all rings
//
// Revision 1.1  2020/09/22 18:14:29  abbott
// Summary: Split off from old (and lengthy) MatrixOps.C; all fn related to dets
//
//
