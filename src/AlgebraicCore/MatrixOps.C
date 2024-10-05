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

#include "CoCoA/CanonicalHom.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MatrixFp.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/matrix.H"
#include "CoCoA/utils.H" // for len


//#include <vector>
using std::vector;

namespace CoCoA
{


  namespace // anonymous
  {
  // WARNING!! Pivot selection strategy is simple rather than clever!
    long RankAndGauss(matrix& M, const int ToDoCols)   // add ErrCtx here
  {
    const ring R = RingOf(M);
    if (!IsField(R)) CoCoA_THROW_ERROR1(ERR::ReqField); // add ErrCtx
    if (ToDoCols > NumCols(M)) CoCoA_THROW_ERROR2(ERR::OutOfRange, "ToDoCols");  // CoCoA_ASSERT?
    const long Mrows = NumRows(M);

    long rank = 0;
    for (long j=0; j < ToDoCols; ++j)
    {
      // Look for a pivot in col j.
      long PivotRow=rank;
      while (PivotRow < Mrows && M(PivotRow, j) == 0)
        ++PivotRow;
      if (PivotRow == Mrows) continue; // col was zero, try next col.
      if (PivotRow != rank)  SwapRows(M, rank, PivotRow);
      M->myRowMul(rank, 1/M(rank, j));  // make pivot entry = 1
      for (long i=0; i < Mrows; ++i)
      {
        CheckForInterrupt("RankAndGauss");
        if (i == rank) continue;
        if (M(i, j) == 0) continue;
        M->myAddRowMul(i, rank, -M(i,j));
      }
      ++rank;
    }
    return rank;
  }
  } // end of namespace anonymous


  std::vector<RingElem> GetRow(ConstMatrixView M, long i)
  {
    if (i < 0 || i >= NumRows(M))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    const long ncols = NumCols(M);
    vector<RingElem> v; v.reserve(ncols);
    for (long j=0; j < ncols; ++j)
      v.push_back(M(i,j));
    return v;
  }

  std::vector<RingElem> GetCol(ConstMatrixView M, long j)
  {
    if (j < 0 || j >= NumCols(M))  CoCoA_THROW_ERROR1(ERR::BadIndex);
    const long nrows = NumRows(M);
    vector<RingElem> v; v.reserve(nrows);
    for (long i=0; i < nrows; ++i)
      v.push_back(M(i,j));
    return v;
  }

  std::vector< std::vector<RingElem> > GetRows(ConstMatrixView M)
  {
    const long nrows = NumRows(M);
    vector< vector<RingElem> > v; v.reserve(nrows);
    for (long i=0; i < nrows; ++i)
      v.push_back(GetRow(M,i));
    return v;
  }

  std::vector< std::vector<RingElem> > GetCols(ConstMatrixView M)
  {
    const long ncols = NumCols(M);
    vector< vector<RingElem> > v; v.reserve(ncols);
    for (long j=0; j < ncols; ++j)
      v.push_back(GetCol(M,j));
    return v;
  }



  // norms moved to MatrixOps-norms.C
  // det moved to MatrixOps-det.C
  // rk,rank moved to MatrixOps-rank.C  
  // adj moved to MatrixOps-adj.C


  matrix inverse(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    if (IsFiniteField(RingOf(M)))
      return InverseFp(M);
    return InverseByGauss(M);
  }


  // Should restriction to full rank be in the name?
  matrix PseudoInverse(ConstMatrixView M)
  {
    // BUG??? Would it make sense to generalize to non fields???
    const ring R = RingOf(M);
    if (!IsField(R))  CoCoA_THROW_ERROR1(ERR::ReqField);

    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    const long rank = rk(M);
    if (rank < Nrows && rank < Ncols) 
      CoCoA_THROW_ERROR2(ERR::NYI, "non full rank matrix");

    // Easy case: a square full rank matrix
    if (Nrows == Ncols)
      return inverse(M);

    if (Nrows < Ncols)
      return transpose(M)*inverse(M*transpose(M));
    else
      return inverse(transpose(M)*M)*transpose(M);
  }


  matrix LinSolve(ConstMatrixView M, ConstMatrixView rhs)
  {
    const ring R = RingOf(M);
    if (RingOf(rhs) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumRows(M) != NumRows(rhs))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    if (IsField(R)) return LinSolveByGauss(M, rhs);
    if (IsTrue3(IsPID3(R))) return LinSolveByHNF(M, rhs);
    if (IsPolyRing(R) && IsField(CoeffRing(R))) return LinSolveByModuleRepr(M, rhs);

    CoCoA_THROW_ERROR2(ERR::NYI, "LinSolve over non-field, non-gcddomain, non-polynomial-ring");
    return LinSolve(M,rhs); // never reached -- just to keep compiler quiet
  }


  matrix LinKer(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (IsField(R)) return LinKerByGauss(M);
    //    if (IsTrue3(IsPID3(R))) return LinSolveByHNF(M, rhs);
    //    if (IsPolyRing(R) && IsField(CoeffRing(R))) return LinSolveByModuleRepr(M, rhs);

    CoCoA_THROW_ERROR2(ERR::NYI, "LinKer over non-field");
    return LinKer(M); // never reached -- just to keep compiler quiet
  }



  matrix LinSolveByGauss(ConstMatrixView M, ConstMatrixView rhs)
  {
    const ring R = RingOf(M);
    if (RingOf(rhs) != R)  CoCoA_THROW_ERROR1(ERR::MixedRings);
    if (NumRows(M) != NumRows(rhs))  CoCoA_THROW_ERROR1(ERR::IncompatDims);
    if (!IsField(R))  CoCoA_THROW_ERROR1(ERR::ReqField);
    const long Mrows = NumRows(M);
    const long Mcols = NumCols(M);
    const long RHScols = NumCols(rhs);
    matrix tmp = NewDenseMat(ConcatHor(M, rhs));

    // Do row redn by Gauss
    const long rank = RankAndGauss(tmp, Mcols);

    // Now tmp has been row reduced, so get answer out.
    matrix ans = NewDenseMat(R, Mcols, RHScols); // initially full of zeroes
    long col=0;
    for (long i=0; i < rank; ++i)
    {
      while (tmp(i,col) == 0) { ++col; }
      for (long j=0; j < RHScols; ++j)
        SetEntry(ans, col, j, tmp(i, j+Mcols));
    }
    for (long i=rank; i < Mrows; ++i)
    {
      for (long j=0; j < RHScols; ++j)
        {
          if (tmp(i, j+Mcols) != 0)
            return NewDenseMat(R,0,0); // to indicate that no soln exists
        }
    }
    return ans;
  }


  matrix LinKerByGauss(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (!IsField(R))  CoCoA_THROW_ERROR1(ERR::ReqField);
    matrix tmp = NewDenseMat(M);
  
    const long Mrows = NumRows(M);
    const long Mcols = NumCols(M);

    // Do row redn by Gauss
    const long rank = RankAndGauss(tmp, Mcols);

    // Now tmp has been row reduced, so get answer out.
    matrix ans = NewDenseMat(R, Mcols, Mcols-rank); // initially full of zeroes
    long row=0;
    long anscol=0;
    vector<long> PivotCols; // the i-th pivot is in col PivotCols[i]
    ConstMatrixView Z = ZeroMat(R, std::max(Mcols-Mrows, (long)0), Mcols);
    ConstMatrixView SqTmp = ConcatVer(tmp, Z); // make it square
    for (long j=0; j < Mcols; ++j) // we consider only Mcols x Mcols
      if (SqTmp(row,j) != 0) // j-th col with pivot
      {
        PivotCols.push_back(j);
        ++row;
      }
      else // j-th col without pivot
      {
        for (long i=0; i < len(PivotCols); ++i)  // "copy" j-th column
          SetEntry(ans, PivotCols[i], anscol, -SqTmp(i, j));
        SetEntry(ans, j, anscol, one(R));
        ++anscol;
      }
    return ans;
  }


  matrix LinSolveByHNF(ConstMatrixView M, ConstMatrixView rhs)
  {
    // HNF works only for PIDs: i.e. ZZ or k[x]
    // NB Could work in k[x,y,z] if M is univariate!
    CoCoA_THROW_ERROR1(ERR::NYI);
    return NewDenseMat(RingOf(M),0,0); // never reached -- just to keep compiler quiet
  }

  
  matrix LinSolveByModuleRepr(ConstMatrixView M, ConstMatrixView rhs)
  {
    // Works in k[x,y,z] where k is a field.  Probably slow.
    CoCoA_THROW_ERROR1(ERR::NYI);
    return NewDenseMat(RingOf(M),0,0); // never reached -- just to keep compiler quiet
  }



  matrix InverseByGauss(ConstMatrixView M)
  {
    // this code works only if the base ring is an integral domain
    CoCoA_ASSERT(IsIntegralDomain(RingOf(M)));
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    const ring R(RingOf(M));
    if (!IsIntegralDomain(R))
      CoCoA_THROW_ERROR1(ERR::NotIntegralDomain);
    const ring    K(IsField(R) ?  R : NewFractionField(R));
    const RingHom R2K(R==K ? IdentityHom(K) : EmbeddingHom(K));
    const long N = NumRows(M);
    matrix Gss(CanonicalHom(R,K)( M ));
    matrix inv = NewDenseMat(IdentityMat(K, N));
    RingElem c(K);
    for (long j=0; j < N; ++j)
    {
      CheckForInterrupt("InverseByGauss");
      if (IsZero(Gss(j,j)))
      {
        long i=j+1;
        for ( ; i < N; ++i)
          if (!IsZero(Gss(i,j))) break;
        if (i == N)  CoCoA_THROW_ERROR1(ERR::NotInvMatrix);
        Gss->mySwapRows(i,j);
        inv->mySwapRows(i,j);
      }
      c = 1/Gss(j,j);
      Gss->myRowMul(j, c);
      inv->myRowMul(j, c);
      for (long i=0; i < N; ++i)
        if (i != j)
        {
          c = -Gss(i,j);
          Gss->myAddRowMul(i, j, c); // AddRowMul(Gss, i, j, c);
          inv->myAddRowMul(i, j, c);
        }
    }
    if (R==K) return inv;
    matrix inv_R = NewDenseMat(R,N,N);
    for (long i=0; i < N; ++i)
      for (long j=0; j < N; ++j)
        SetEntry(inv_R, i, j, num(inv(i,j))/den(inv(i,j)));
    return inv_R;
  }



//////////////////////////////////////////////////////////////////



  bool IsZero(const ConstMatrixView& M)
  { return M == ZeroMat(RingOf(M), NumRows(M), NumCols(M)); }


  bool IsZeroRow(const ConstMatrixView& M, long i)
  {
    M->myCheckRowIndex(i, CoCoA_ERROR_CONTEXT);
    return M->myIsZeroRow(i);
  }


  bool IsZeroCol(const ConstMatrixView& M, long j)
  {
    M->myCheckColIndex(j, CoCoA_ERROR_CONTEXT);
    return M->myIsZeroCol(j);
  }


  bool IsSymmetric(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    return M->IamSymmetric();
  }


  bool IsAntiSymmetric(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    return M->IamAntiSymmetric();
  }


  bool IsDiagonal(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    return M->IamDiagonal();
  }

  
  bool IsUpperTriangular(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    return M->IamUpperTriangular();
  }

  
  bool IsLowerTriangular(const ConstMatrixView& M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
    return M->IamLowerTriangular();
  }

  
  bool IsMat0x0(const ConstMatrixView& M)
  {
    return NumRows(M) == 0 && NumCols(M) == 0;
  }


  bool HasNegEntry(const ConstMatrixView& M)
  {
    return M->IhaveNegEntry();
  }
  

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps.C,v 1.28 2024/08/28 15:25:56 bigatti Exp $
// $Log: MatrixOps.C,v $
// Revision 1.28  2024/08/28 15:25:56  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2
//
// Revision 1.27  2024/08/02 15:59:13  bigatti
// Summary: using CoCoA_THROW_ERROR1/2;   and some comments
//
// Revision 1.26  2024/07/15 17:03:57  bigatti
// Summary: Error codes:  added CoCoA_THROW_ERROR1(ErrID);
// Removed BadRowIndex, BadColIndex;  Renamed  ERR::Empty --> ERR::ReqNonEmpty
//
// Revision 1.25  2024/07/11 15:42:42  bigatti
// *** empty log message ***
//
// Revision 1.24  2024/07/02 15:31:12  bigatti
// Summary: Changed error codes: LogZero into ReqNonZero
//
// and Not... into ReqNonNegative, ReqNonNegativeGrading, ReqPositive, ReqPositiveGrading
//
// Revision 1.23  2024/03/16 10:15:12  abbott
// Summary: Use ErrorContext when calling arg-check fns (redmine 1793)
//
// Revision 1.22  2023/10/25 19:22:49  abbott
// Summary: inverse now calls InverseFp if mat is over Fp; added CheckForInterrupt to InverseByGauss
//
// Revision 1.21  2022/08/05 20:44:48  abbott
// Summary: Took NumRows/NumCols out of for-loop exit conds
//
// Revision 1.20  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.19  2021/10/04 08:56:29  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.18  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.17  2020/09/28 11:13:37  abbott
// Summary: Further splitting of MatrixOps (redmine 1196)
//
// Revision 1.16  2020/09/22 18:13:53  abbott
// Summary: Renamed FrobeniusNorm2 to FrobeniusNormSq; for det stuff into new file MatrixOps-det.C
//
// Revision 1.15  2020/07/28 08:03:00  abbott
// Summary: MAde determinant interruptible
//
// Revision 1.14  2020/06/22 15:43:37  abbott
// Summary: Use new CoCoA_STATIC_ERROR_MESG macro
//
// Revision 1.13  2020/06/20 19:12:16  abbott
// Summary: AutomaticConversionHom now requires 3rd arg (FnName of caller)
//
// Revision 1.12  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.11  2020/02/27 10:53:47  abbott
// Summary: Added new fn IsZeroDet
//
// Revision 1.10  2020/02/13 16:34:10  abbott
// Summary: Added PowerOfDiagMat (in anon namespace)
//
// Revision 1.9  2020/02/12 09:15:10  abbott
// Summary: Renamed RowReducedForm to rref
//
// Revision 1.8  2020/01/26 14:42:28  abbott
// Summary: Revised includes after splitting NumTheory
//
// Revision 1.7  2019/09/27 14:43:39  abbott
// Summary: Cleaned DetOfSmallMat
//
// Revision 1.6  2019/09/16 17:26:57  abbott
// Summary: Added GetRow, GetRows, GetCol, GetCols, DetOfSmallMat, RowReducedForm
//
// Revision 1.5  2019/03/18 11:16:27  abbott
// Summary: Added include (after changing NumTheory); commented out (unused) calls to CpuTime
//
// Revision 1.4  2018/06/26 10:01:10  abbott
// Summary: Increased VerbosityLevel needed to print blurb
//
// Revision 1.3  2018/06/12 13:55:33  abbott
// Summary: Some cleaning
//
// Revision 1.2  2018/05/18 12:15:04  bigatti
// -- renamed IntOperations --> BigIntOps
//
// Revision 1.1  2018/05/17 15:25:53  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.24  2018/04/18 14:27:15  abbott
// Summary: Minor cleaning
//
// Revision 1.23  2018/03/02 13:44:11  abbott
// Summary: Major revision to conversion mat(QQ) --> mat(ZZ) for determinant
//
// Revision 1.22  2018/02/27 17:30:22  abbott
// Summary: Renamed NumTheory_prime to NumTheory-prime; changed includes
//
// Revision 1.21  2018/02/27 10:58:59  abbott
// Summary: Added include NumTheory_prime; added first impl for converting det over QQ to det over ZZ
//
// Revision 1.20  2018/02/15 17:03:53  abbott
// Summary: Removed some verbosity
//
// Revision 1.19  2018/02/12 14:52:11  abbott
// Summary: Revised det2x2, det3x3, DetByGauss.  Added det4x4, det5x5, and DetByCRT
//
// Revision 1.18  2017/04/07 14:20:43  bigatti
// -- 2 small changes to unused lines to keep compiler quiet
//
// Revision 1.17  2016/10/27 14:05:39  abbott
// Summary: Added code for det of 0x0 and 1x1 matrices; see issue 956
//
// Revision 1.16  2016/06/10 12:01:12  bigatti
// -- just a change of sign to the result of LinKer
//
// Revision 1.15  2015/12/11 15:51:23  bigatti
// -- added IsLowerTriangular
//
// Revision 1.14  2015/12/11 13:03:20  bigatti
// -- added IsUpperTriangular, IhaveNegEntry
// -- removed useless checks
//
// Revision 1.13  2015/06/26 14:59:02  abbott
// Summary: Corrected assertion inside InverseByGauss
// Author: JAA
//
// Revision 1.12  2015/05/20 13:27:18  abbott
// Summary: adj now calls AdjByDetOfMinors if det is zero
// Author: JAA
//
// Revision 1.11  2015/05/19 07:26:36  abbott
// Summary: Minor improvement to DetDirect
// Author: JAA
//
// Revision 1.10  2015/05/15 15:59:04  bigatti
// -- fixed swapped arguments in DetDirect
//
// Revision 1.9  2015/05/11 15:49:03  bigatti
// -- now InverseByGauss works also if matrix is over IntegralDomain (and invertible)
//
// Revision 1.8  2015/04/27 10:08:48  bigatti
// Summary: changed myAddMul --> myAddMulLM
//
// Revision 1.7  2015/04/13 16:16:12  abbott
// Summary: Changed "rank" --> "rk"; adjoint" --> "adj"; added "AdjDirect", "DetDirect"
// Author: JAA
//
// Revision 1.6  2014/08/26 12:55:58  abbott
// Summary: Cleaned up DetByGauss; added DetByBareiss
// Author: JAA
//
// Revision 1.5  2014/07/30 14:06:24  abbott
// Summary: Changed BaseRing into RingOf
// Author: JAA
//
// Revision 1.4  2014/07/08 08:35:16  abbott
// Summary: Removed AsFractionField
// Author: JAA
//
// Revision 1.3  2014/07/07 12:23:22  abbott
// Summary: Removed AsPolyRing
// Author: JAA
//
// Revision 1.2  2014/04/17 13:38:47  bigatti
// -- MatrixViews --> MatrixView
//
// Revision 1.1  2014/04/11 15:42:37  abbott
// Summary: Renamed from MatrixArith
// Author: JAA
//
// Revision 1.46  2014/04/08 15:40:50  abbott
// Summary: Replaced test for IsZero by IsZeroDivisor
// Author: JAA
//
// Revision 1.45  2014/01/16 16:10:56  abbott
// Removed a blank line.
//
// Revision 1.44  2012/11/23 17:28:26  abbott
// Modified LinSolve so that it returns a 0x0 matrix when no soln exists
// (previously it threw an exception).
//
// Revision 1.43  2012/10/16 09:45:44  abbott
// Replaced RefRingElem by RingElem&.
//
// Revision 1.42  2012/10/03 15:25:05  abbott
// Replaced swap by assignment in DetByGauss; new impl of swap did not work
// in that instance (since one value was a temporary).
//
// Revision 1.41  2012/07/10 12:59:34  bigatti
// -- added two lines to keep compiler quiet
//
// Revision 1.40  2012/07/10 09:48:41  bigatti
// -- fixes of some naive errors
//
// Revision 1.39  2012/07/10 09:23:30  bigatti
// -- separated gauss code from LinSolveByGauss
// -- added LinKerByGauss
//
// Revision 1.38  2012/06/19 15:43:27  abbott
// Added division of a matrix by a scalar.
//
// Revision 1.37  2012/06/11 08:20:33  abbott
// Added multiplication on the right by a scalar.
//
// Revision 1.36  2012/06/10 22:57:31  abbott
// Added negation for matrices -- same as doing (-1)*M.
//
// Revision 1.35  2012/05/30 16:04:55  bigatti
// -- applied "3" convention on bool3 functions and member fields
//
// Revision 1.34  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.33  2012/05/28 09:18:21  abbott
// Created IntOperations which gathers together all operations on
// integers (both big and small).  Many consequential changes.
//
// Revision 1.32  2012/05/04 15:39:29  abbott
// Corrected a comment.
//
// Revision 1.31  2012/04/27 14:49:33  abbott
// Added LinSolve family (incl. LinSolveByGauss, LinSolveByHNF, LinSolveByModuleRepr).
//
// Revision 1.30  2012/04/16 09:21:19  abbott
// Added missing include directive.
//
// Revision 1.29  2012/04/13 16:24:35  abbott
// Added solve and SolveByGauss.
//
// Revision 1.28  2012/04/11 14:03:24  abbott
// Very minor change: slight improvement to readability.
//
// Revision 1.27  2012/03/16 14:42:46  bigatti
// -- fixed AdjointByDetOfMinors
// -- fixed adjoint (field + det(M)=0)
//
// Revision 1.26  2012/02/10 10:26:39  bigatti
// -- changed RingZ.H, RingQ.H --> RingZZ.H, RingQQ.H
//
// Revision 1.25  2011/11/09 14:09:53  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.24  2011/08/24 10:25:53  bigatti
// -- renamed QQ --> BigRat
// -- sorted #include
//
// Revision 1.23  2011/08/14 15:52:17  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.22  2011/05/13 16:47:20  abbott
// Added power fn for matrices: partial impl, cannot yet handle negative powers.
//
// Revision 1.21  2011/05/03 13:48:02  abbott
// Now using CanonicalHom inside DetByGauss.
// Cleaner and avoids a mysterious compiler warning.
//
// Revision 1.20  2011/03/22 16:44:19  bigatti
// -- fixed check in det
//
// Revision 1.19  2011/03/16 15:41:06  bigatti
// -- minor cleaning
//
// Revision 1.18  2011/03/10 11:25:46  bigatti
// -- now using long instead of size_t, and len(v) instead of v.size()
//
// Revision 1.17  2011/03/03 13:50:22  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.16  2011/03/01 14:13:24  bigatti
// -- added f*M
//
// Revision 1.15  2011/02/28 14:08:49  bigatti
// -- added det3x3
// -- using apply mapping matrix (in DetByGauss)
//
// Revision 1.14  2011/02/10 15:27:06  bigatti
// -- commented #include vector  (included in MatrixArith.H)
//
// Revision 1.13  2011/02/09 16:48:27  bigatti
// -- added + and - for matrices
//
// Revision 1.12  2009/12/23 18:55:16  abbott
// Removed some useless comments.
//
// Revision 1.11  2009/12/03 17:40:36  abbott
// Added include directives for ZZ.H (as a consequence of removing
// the directive from ring.H).
//
// Revision 1.10  2009/06/25 16:59:42  abbott
// Minor improvement to some error messages (better coherence & comprehensibility).
//
// Revision 1.9  2008/07/09 16:09:11  abbott
// Removed pointless bogus function declarations.
//
// Revision 1.8  2008/04/22 14:42:03  abbott
// Two very minor changes.
//
// Revision 1.7  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
// Revision 1.6  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.5  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.4  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.3  2007/10/30 17:14:08  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/10/30 15:54:15  bigatti
// -- fixed index too big in RankByGauss(ConstMatrix M)
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.10  2007/03/08 18:22:29  cocoa
// Just whitespace cleaning.
//
// Revision 1.9  2007/03/05 21:06:07  cocoa
// New names for homomorphism pseudo-ctors: removed the "New" prefix.
//
// Revision 1.8  2007/03/02 10:47:53  cocoa
// First stage of RingZ modifications -- tests do not compile currently, Anna will fix this.
//
// Revision 1.7  2007/02/10 18:44:03  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.6  2006/12/21 13:48:32  cocoa
// Made all increment/decrement calls prefix (except where the must be postfix).
//
// Revision 1.5  2006/11/27 16:18:32  cocoa
// -- moved classes declarations from .H to .C (DenseMatrix, DiagMatrix,
//    FieldIdeal, SpecialMatrix)
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/08/17 09:39:07  cocoa
// -- added: elimination ordering matrix for non-homogeneous input
//
// Revision 1.2  2006/07/17 16:58:05  cocoa
// -- added: NewMatrixElim(size_t NumIndets, std::vector<size_t> IndetsToElim)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.8  2006/05/02 14:39:20  cocoa
// -- Changed "not" into "!" becuase of M$Windoze (by M.Abshoff)
//
// Revision 1.7  2006/04/27 13:45:30  cocoa
// Changed name of NewIdentityRingHom to NewIdentityHom.
// Changed name of member functions which print out their own object
// into myOutputSelf (to distinguish from "transitive" myOutput fns).
//
// Revision 1.6  2006/04/10 13:20:43  cocoa
// -- fixed buglets for Elimination orderings
//
// Revision 1.5  2006/04/05 16:45:29  cocoa
// -- added comment and CoCoA_ASSERT in NewPositiveMatrix
// -- added IsPositiveGrading
//
// Revision 1.4  2006/04/05 14:49:20  cocoa
// -- fixed: NewPositiveMatrix (tested and used in OrdvArith.C)
//
// Revision 1.3  2006/01/19 15:48:49  cocoa
// -- fixed RankByGauss by Stefan Kaspar
//
// Revision 1.2  2005/12/31 12:22:17  cocoa
// Several minor tweaks to silence the Microsoft compiler:
//  - added some missing #includes and using directives
//  - moved some function defns into the right namespace
//  - etc.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:31  cocoa
// Imported files
//
// Revision 1.5  2005/04/20 15:40:48  cocoa
// Major change: modified the standard way errors are to be signalled
// (now via a macro which records filename and line number).  Updated
// documentation in error.txt accordingly.
//
// Improved the documentation in matrix.txt (still more work to be done).
//
// Revision 1.4  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/30 17:15:14  cocoa
// Cleaned the SpecialMatrix code; a simple test compiles and
// runs fine.
//
// Revision 1.1  2005/03/11 18:38:32  cocoa
// -- first import
//
