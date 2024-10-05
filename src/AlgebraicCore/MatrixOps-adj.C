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

#include "CoCoA/BigIntOps.H"
#include "CoCoA/DenseMatrix.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/MatrixView.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/combinatorics.H"
#include "CoCoA/error.H"
#include "CoCoA/matrix.H"


//#include <vector>
using std::vector;
#include<algorithm>
using std::next_permutation;

namespace CoCoA
{


  matrix adj(ConstMatrixView M)
  {
    if (NumRows(M) != NumCols(M))  CoCoA_THROW_ERROR1(ERR::ReqSquareMatrix);
///JAA???return AdjDirect(M);
///JAA???return AdjByDetOfMinors(M);
    if (!IsIntegralDomain(RingOf(M)) || (IsPolyRing(RingOf(M)) && NumIndets(RingOf(M)) > 1))
      return AdjByDetOfMinors(M);
    if (!IsInvertible(det(M)))
      return AdjByDetOfMinors(M);
    return AdjByInverse(M);
  }


  matrix AdjByInverse(ConstMatrixView M)
  {
    // This code works only if the matrix is invertible!
    // It is morally equivalent to swap(lhs, inverse(M)*det(M));
    CoCoA_ASSERT(IsIntegralDomain(RingOf(M))); // needed for this method, not mathematically necessary BUG BUG!!!
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    if (IsField(RingOf(M)))
    {
      const RingElem d = det(M);
      if (IsZero(d)) CoCoA_THROW_ERROR(ERR::DivByZero, "AdjByInverse");
      return d*inverse(M);
    }
    FractionField K(NewFractionField(RingOf(M)));
    RingHom R2K(EmbeddingHom(K));
    const matrix ans_K = adj(R2K(M));
    const long Nrows = NumRows(M);
    const long Ncols = NumCols(M);
    matrix ans = NewDenseMat(RingOf(M), Nrows, Ncols);
    for (long i=0; i < Nrows; ++i)
      for (long j=0; j < Ncols; ++j)
        SetEntry(ans, i,j, num(ans_K(i,j)));
    return ans;
  }


  // Simple (but probably not very fast).
  matrix AdjByDetOfMinors(ConstMatrixView M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    matrix adj = NewDenseMat(RingOf(M), n,n);
    vector<long> rows(n-1);
    for (long i=0; i < n-1; ++i) { rows[i] = i+1; }

    vector<long> cols(n-1);
    for (long i=0; i < n; ++i)
    {
      if (i > 0) rows[i-1] = i-1;
      for (long j=0; j < n-1; ++j) { cols[j] = j+1; }
      for (long j=0; j < n; ++j)
      {
        if (j > 0) cols[j-1] = j-1;
        SetEntry(adj, j,i, SmallPower(-1,i+j)*det(submat(M,rows,cols)));
      }
    }
    return adj;
  }


  matrix AdjDirect(ConstMatrixView M)
  {
    CoCoA_ASSERT(NumRows(M) == NumCols(M));
    const long n = NumRows(M);
    const ring R = RingOf(M);
    matrix adj = NewDenseMat(R, n,n);

    for (int col = 0; col < n; ++col)
    {
      for (int row = 0; row < n; ++row)
      {
        vector<int> rows(n-1);
        for (long i=0; i < n-1; ++i) { rows[i] = i + (i >= row); }

        RingElem entry(R);

        do
        {
          RingElem term = one(R);
          if ((row+col)%2 == 1) term = -term;
          for (int i=0; i < n-1; ++i)
            term *= M(rows[i], i+(i>=col));
          entry += signature(rows)*term;
        } while ( std::next_permutation(rows.begin(),rows.end()) );

        SetEntry(adj, col,row, entry); // transposed!!
      }
    }
    return adj;
  }


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-adj.C,v 1.6 2024/08/28 15:23:58 bigatti Exp $
// $Log: MatrixOps-adj.C,v $
// Revision 1.6  2024/08/28 15:23:58  bigatti
// Summary: now using CoCoA_THROW_ERROR1/2
//
// Revision 1.5  2022/02/18 14:11:54  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.4  2021/10/04 08:54:44  abbott
// Summary: Replaced calls to apply by direct applic of ringhom (redmine 1467, 1598)
//
// Revision 1.3  2021/01/07 15:16:51  abbott
// Summary: Corrected copyright
//
// Revision 1.2  2020/10/06 19:15:02  abbott
// Summary: Replaced power(-1,...) by SmallPower(-1,...)
//
// Revision 1.1  2020/09/28 11:13:36  abbott
// Summary: Further splitting of MatrixOps (redmine 1196)
//
//
