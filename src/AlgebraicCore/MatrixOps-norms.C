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

#include "CoCoA/MachineInt.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"


//#include <vector>
using std::vector;

namespace CoCoA
{

  // The square of the Frobenius norm of a matrix.
  RingElem FrobeniusNormSq(ConstMatrixView M)
  {
    RingElem FrNorm2 = zero(RingOf(M));
    const long nrows = NumRows(M);
    const long ncols = NumCols(M);
    for (long i=0; i < nrows; ++i)
      for (long j=0; j < ncols; ++j)
        FrNorm2 += M(i,j)*M(i,j);
    return FrNorm2;
  }


  // Compute the induced infty-norm of a matrix
  RingElem OperatorNormInfinity(ConstMatrixView M)
  {
    const ring R = RingOf(M);
    if (!IsOrderedDomain(R))
      CoCoA_THROW_ERROR(ERR::ReqOrdDom, "OperatorNormInfinity(Mat)");
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    RingElem MaxRowSize = zero(R);

    for (long i=0; i < Nrows; ++i)
    {
      RingElem RowSize(zero(R));
      for (long j=0; j < Ncols; ++j) {	RowSize += abs(M(i,j)); }
      if (RowSize > MaxRowSize)
        MaxRowSize = RowSize;
    }
    return MaxRowSize;
  }


  RingElem OperatorNorm1(ConstMatrixView M)
  {
    return OperatorNormInfinity(transpose(M));
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-norms.C,v 1.5 2024/07/02 20:42:07 abbott Exp $
// $Log: MatrixOps-norms.C,v $
// Revision 1.5  2024/07/02 20:42:07  abbott
// Summary: Renaming of errors (redmine 308)
//
// Revision 1.4  2022/08/05 20:44:30  abbott
// Summary: Took NumRows/NumCols out of for-loop exit conds
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
