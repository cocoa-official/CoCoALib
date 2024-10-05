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

#include "CoCoA/MachineInt.H"
#include "CoCoA/error.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/matrix.H"
#include "CoCoA/ring.H"


namespace CoCoA
{

// OBSOLETE -- new rref is "better"

  // matrix RowReducedForm(ConstMatrixView M_orig)
  // {
  //   if (!IsField(RingOf(M_orig))) CoCoA_THROW_ERROR(ERR::NYI, "row reduction over non field");
  //   matrix M = NewDenseMat(M_orig);
  //   if (NumRows(M) < 2) return M;
  //   const long ncols = NumCols(M);
  //   const long nrows = NumRows(M);
  //   long OutputRow = 0;
  //   for (long j=0; j < ncols; ++j)
  //   {
  //     // Pick pivot row; if poss, pick +1 or -1 as pivot
  //     long PivotRow = -1; // -1 means no non-zero pivot found
  //     for (long i=OutputRow; i < nrows; ++i)
  //     {
  //       if (IsZero(M(i,j))) continue;
  //       if (IsOne(M(i,j)) || IsMinusOne(M(i,j))) { PivotRow = i; break; }
  //       if (PivotRow == -1) PivotRow = i;
  //     }
  //     if (PivotRow == -1) continue; // no pivot in this column
  //     if (PivotRow != OutputRow)
  //       SwapRows(M, OutputRow, PivotRow); // now pivot row has index OutputRow
  //     PivotRow = OutputRow;
  //     // Make rest of col zero
  //     for (long i=OutputRow+1; i < nrows; ++i)
  //       AddRowMul(M,i, PivotRow, -M(i,j)/M(PivotRow,j));
  //     ++OutputRow;
  //   }
  //   return M;
  // }

  //-------------------------------------------------------


  matrix rref(ConstMatrixView M)
  {
    return RREFByGauss(M);
  }
  
  // Simple rather than smart/fast
  matrix RREFByGauss(ConstMatrixView M_orig)
  {
    if (!IsField(RingOf(M_orig))) CoCoA_THROW_ERROR(ERR::ReqField, "RREFByGauss");
    matrix M = NewDenseMat(M_orig);
    const int nrows = NumRows(M);
    const int ncols = NumCols(M);
    int CurrCol = 0;
    for (int r=0; r < nrows; ++r)
    {
      int PivotRow;
      do
      {
        PivotRow = r;
        while (PivotRow < nrows && IsZero(M(PivotRow,CurrCol))) { ++PivotRow; }
        if (PivotRow < nrows) break;
        ++CurrCol;
      }
      while (CurrCol < ncols);
      if (CurrCol == ncols) return M;

      SwapRows(M, r, PivotRow);  // does nothing if r == PivotRow
      const RingElem pivot = M(r,CurrCol);
      for (int j=CurrCol; j < ncols; ++j)
        SetEntry(M, r,j, M(r,j)/pivot);
      for (int i=0; i < nrows; ++i)
      {
        if (i == r) continue;
        // Next lines are efficient form of AddRowMul(M, i, r, -M(i,CurrCol));
        // (efficient because skips cols before CurrCol)
        if (IsZero(M(i,CurrCol))) continue;
        const RingElem scale = M(i,CurrCol);
        for (int j=CurrCol; j < ncols; ++j)
          SetEntry(M, i,j, M(i,j)-scale*M(r,j));
      }
    }
    return M;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-rref.C,v 1.8 2024/07/11 15:42:42 bigatti Exp $
// $Log: MatrixOps-rref.C,v $
// Revision 1.8  2024/07/11 15:42:42  bigatti
// *** empty log message ***
//
// Revision 1.7  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2020/06/17 15:49:23  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.5  2020/02/12 10:38:56  abbott
// Summary: Added missing include
//
// Revision 1.4  2020/02/12 09:15:10  abbott
// Summary: Renamed RowReducedForm to rref
//
// Revision 1.3  2020/01/27 19:55:16  abbott
// Summary: Some comments; minor tidying
//
// Revision 1.2  2020/01/26 14:33:52  abbott
// Summary: Added include MachineInt
//
// Revision 1.1  2020/01/09 18:35:41  abbott
// Summary: Added rref (Row reduced echelon form)
//
//
