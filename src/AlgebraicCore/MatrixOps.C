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
