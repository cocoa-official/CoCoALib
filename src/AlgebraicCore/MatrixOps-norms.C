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
      CoCoA_THROW_ERROR1(ERR::ReqOrdDom);
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    RingElem MaxRowSize = zero(R);

    for (long i=0; i < Nrows; ++i)
    {
      RingElem RowSize(zero(R));
      for (long j=0; j < Ncols; ++j)  { RowSize += abs(M(i,j)); }
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
