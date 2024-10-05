//   Copyright (c)  2020  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/error.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/BigRat.H"
#include "CoCoA/BigRatOps.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/matrix.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/random.H"
#include "CoCoA/utils.H"


//#include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {
    // Returns -1 to mean no non-zero pivot possible, o/w returns a row index
    int PickPivotRow(const vector<vector<BigInt>>& M, int col, const vector<bool>& RowUsed)
    {
      const int nrows = len(M);
  
      int MaxRow = -1;
      int count = 0;
      for (int i=0; i < nrows; i++)
      {
        if (RowUsed[i] || IsZero(M[i][col])) continue;
        ++count;
        if (MaxRow == -1 || CmpAbs(M[MaxRow][col], M[i][col]) < 0)
          MaxRow = i;
      }
      if (count < 2) return MaxRow; // this is -1 if all entries are 0, o/w only 1 non-zero pivot

      // Pick randomly an unused row with non-zero entry, avoiding MaxRow
      count = RandomLong(1, count-1);
      for (int PivotRow=0; PivotRow < nrows; ++PivotRow)
      {
        if (RowUsed[PivotRow] || IsZero(M[PivotRow][col]) || PivotRow == MaxRow)
          continue;
        if (--count == 0) return PivotRow;
      }
      return -99; // NEVER EXECUTED; just to keep compiler quiet
    }

  
    // WARNING!  Modifies its arg!!
    // Mindlessly translated from ancient C code in TmpFactorDir/linalg/Zkernel.c
    // Randomized algorithm!
    vector<vector<BigInt>> ZZKernel(vector<vector<BigInt>>& M)
    {
      const int nrows = len(M);
      const int ncols = len(M[0]);
      vector<vector<BigInt>> U(nrows, vector<BigInt>(nrows));
      vector<bool> RowUsed(nrows); // initially all false
      int NumRowsRemaining = nrows;
      for (int i=0; i < nrows; ++i) U[i][i] = 1;
      BigInt q; // used in inner loop (see call to quorem)
      for (int col=0; col < ncols; ++col)
      {
        if (NumRowsRemaining == 0) break;
        int PivotRow = PickPivotRow(M, col, RowUsed);
        if (PivotRow < 0) continue; // happens only if col is zero
        bool finished = false;
        while (!finished)
        {
          if (M[PivotRow][col] < 0)
          {
            for (int j=col; j < ncols; ++j) negate(M[PivotRow][j]);
            for (int j=0;   j < nrows; ++j) negate(U[PivotRow][j]);
          }
          finished = true;
          for (int i=0; i < nrows; ++i)
          {
            if (RowUsed[i] || i == PivotRow) continue;
            quorem(q, M[i][col], M[i][col], M[PivotRow][col]);
            if (IsZero(q)) continue;
            finished = false;
            for (int j=0; j < nrows; ++j)
              U[i][j] -= q*U[PivotRow][j];
            for (int j=col+1; j < ncols; ++j)
              M[i][j] -= q*M[PivotRow][j];
          }
          PivotRow = PickPivotRow(M, col, RowUsed);
        }
        RowUsed[PivotRow] = true;
        --NumRowsRemaining;
      }
      // At this point dimension of kernel is rows_remaining
      // Move the Z-basis of the kernel into the first rows of U

      int j = 0;
      for (int i=0; i < NumRowsRemaining; ++i, ++j)
      {
        while (RowUsed[j]) ++j;
        swap(U[i],U[j]);
      }
      U.resize(NumRowsRemaining);
      return U;
    }

  } // end of namespace anonymous


  // Allow matrix with integer or rational entries
  // Result is full-rank matrix K (over ZZ) such that M*K = 0
  matrix LinKerZZ(ConstMatrixView M)
  {
    const ring& R = RingOf(M);
    if (!IsZero(characteristic(R))) CoCoA_THROW_ERROR("char not 0", "LinKerZZ");
    const int nrows = NumRows(M);
    const int ncols = NumCols(M);
    // Convert transpose(M) into vec<vec<BigInt>>
    vector<vector<BigInt>> M_ZZ(ncols, vector<BigInt>(nrows));
    vector<BigInt> D(ncols);
    vector<BigRat> col_j(nrows);
    for (int j=0; j < ncols; ++j)
    {
      for (int i=0; i < nrows; ++i)
        col_j[i] = ConvertTo<BigRat>(M(i,j)); // may throw
      D[j] = CommonDenom(col_j);
      for (int i=0; i < nrows; ++i)
        M_ZZ[j][i] = num(col_j[i])*(D[j]/den(col_j[i]));
    }
    auto K = ZZKernel(M_ZZ);
    if (len(K) > 0)
    {
      // Have to rescale kernel entries by common denoms
      // WARNING K is TRANSPOSED!!!
      const int krows = len(K);
      const int kcols = ncols;
      for (int i=0; i < krows; ++i)
        for (int j=0; j < kcols; ++j)
          K[i][j] *= D[j];  // assume that *= is clever if rhs==1
    }
    return NewDenseMatTranspose(RingZZ(), K);
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-LinKerZZ.C,v 1.5 2022/02/18 14:11:54 abbott Exp $
// $Log: MatrixOps-LinKerZZ.C,v $
// Revision 1.5  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2022/02/07 17:30:51  abbott
// Summary: Fixed bug (redmine 1658)
//
// Revision 1.3  2021/04/26 13:57:14  abbott
// Summary: Fixed incorrect vector size; changed name of variable; redmine 1590
//
// Revision 1.2  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.1  2020/02/28 08:56:18  abbott
// Summary: Added LinKerZZ
//
//
