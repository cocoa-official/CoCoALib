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

#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MachineInt.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/error.H"
#include "CoCoA/interrupt.H"


//#include <vector>
using std::vector;


namespace CoCoA
{

  long rk(const ConstMatrixView& M)
  {
    if (!IsIntegralDomain(RingOf(M)))
      CoCoA_THROW_ERROR(ERR::NotIntegralDomain, "rk(Mat)");
    return M->myRank();
  }


  // // WARNING!! Pivot selection strategy is simple rather than clever!
  // long RankAndGauss(matrix& M, const int ToDoCols)
  // {
  //   const ring R = RingOf(M);
  //   if (!IsField(R)) CoCoA_THROW_ERROR(ERR::NotField, "gauss");  
  //   if (ToDoCols > NumCols(M)) CoCoA_THROW_ERROR(ERR::BadColIndex, "gauss");  
  //   const long Mrows = NumRows(M);

  //   long rank = 0;
  //   for (long j=0; j < ToDoCols; ++j)
  //   {
  //     // Look for a pivot in col j.
  //     long PivotRow=rank;
  //     while (PivotRow < Mrows && M(PivotRow, j) == 0)
  //       ++PivotRow;
  //     if (PivotRow == Mrows) continue; // col was zero, try next col.
  //     if (PivotRow != rank)  SwapRows(M, rank, PivotRow);
  //     M->myRowMul(rank, 1/M(rank, j));  // make pivot entry = 1
  //     for (long i=0; i < Mrows; ++i)
  //     {
  //       CheckForInterrupt("RankAndGauss");
  //       if (i == rank) continue;
  //       if (M(i, j) == 0) continue;
  //       M->myAddRowMul(i, rank, -M(i,j));
  //     }
  //     ++rank;
  //   }
  //   return rank;
  // }



  long RankByGauss(std::vector<long>& IndepRows, ConstMatrixView M)
  {
    // ??? this function works only within integral domains!
    const ring R(RingOf(M));
    if (!IsIntegralDomain(R))
      CoCoA_THROW_ERROR(ERR::NotIntegralDomain, "RankByGauss(v,Mat)  over non-integral domain");
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    if (Nrows == 0 || Ncols==0) return 0;
    // Henceforth M has at least 1 row & 1 col
    const ring K((IsField(R) ? R : NewFractionField(R)));
    const RingHom R2K(R==K ? IdentityHom(K) : EmbeddingHom(K));
    matrix Gss(NewDenseMat(K, Nrows, Ncols));

    // below copy M via RingHom into Gss
    {
      // This will eventually become a separate function
      for (long row=0; row < Nrows; ++row)
        for (long col=0; col < Ncols; ++col)
          SetEntry(Gss, row, col, R2K(M(row,col)));
    }
    IndepRows.clear();
    RingElem pivot(K);
    long row=0;
    for (long col=0; col < Ncols; ++col)
    {
      if (IsZero(Gss(row,col)))
      {
        long i=row+1;
        for ( ; i < Nrows; ++i)
          if (!IsZero(Gss(i,col)))
          {
            Gss->mySwapRows(i,row);
            break;
          }
        if (i==Nrows) continue;
      }
      IndepRows.push_back(row);
      pivot = Gss(row,col);
      for (long i=row+1; i < Nrows; ++i)
      {
        CheckForInterrupt("RankByGauss");
        Gss->myAddRowMul(i, row, -Gss(i,col)/pivot);
      }
      ++row;
      if (row == Nrows) return row;
    }
    return row;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-rank.C,v 1.4 2023/01/01 11:38:59 abbott Exp $
// $Log: MatrixOps-rank.C,v $
// Revision 1.4  2023/01/01 11:38:59  abbott
// Summary: Fixed bug when matrix has 0 cols (redmine ??)
//
// Revision 1.3  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.1  2020/09/28 11:13:37  abbott
// Summary: Further splitting of MatrixOps (redmine 1196)
//
//
