//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Original author: 2015  Mario Albert

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

#include "CoCoA/TmpResolutionMinimization.H"
#include "CoCoA/matrix.H"

using std::make_pair;

namespace CoCoA
{

  std::pair<long, long> ResolutionMinimization::myFindPivot(ConstMatrixView m) const
  {
    const long nrows = NumRows(m);
    const long ncols = NumCols(m);
    //iterating through complete matrix
    for (long j = 0; j < ncols; ++j)
    {
      for (long i = 0; i < nrows; ++i)
      {
        const RingElemAlias mij = m->myEntry(i,j);
//        if (!myRing->myIsZero(raw(mij)) && myRing->myIsConstant(raw(mij)))  // SLIGHTLY FASTER
        // if we found invertible element return his position
        if (myRing->myIsInvertible(raw(mij)))
//        if (IsInvertible(m(i, j)))
        {
          return make_pair(i, j);
        }
      }
    }
    // nothing found return (-1, -1)
    return make_pair(-1, -1);
  }


  void ResolutionMinimization::myManipulateMatrix(matrix& m, long r, const std::vector<RingElem>& PivotColumn) const
  {
    const long nrows = NumRows(m);
    CoCoA_ASSERT(0 <= r && r < nrows);
    CoCoA_ASSERT(len(PivotColumn) == nrows);
    CoCoA_ASSERT(IsInvertible(PivotColumn[r]));
    // m behaves as a pointer -> no need to give something back
    const RingElem ScaleFactor = (-1)/PivotColumn[r];
    for (long i = 0; i < nrows; ++i)
    {
      if (i != r)
      {
        m->myAddRowMul(i, r, ScaleFactor*PivotColumn[i]);
      }
    }
  }


  void ResolutionMinimization::myMinimization()
  {
    std::vector<matrix>::iterator begin(myResolution.begin());
    ++begin;
    std::vector<matrix>::iterator resEnd(myResolution.end());
    // iterating through the resolution, starting with the second one
    for (std::vector<matrix>::iterator it = begin; it != myResolution.end(); ++it)
    {
      // until we find a PivotElemen reduce resolution
      std::pair<long, long> PivotElement(myFindPivot(*it));
      while (PivotElement.first != -1)
      {
        // delete the pivot column & row and manipulate the rest of the matrix
        std::vector<RingElem> PivotColumn(myGetColumn(*it, PivotElement.second));
        DeleteCol(*it, PivotElement.second);
        myManipulateMatrix(*it, PivotElement.first, PivotColumn);
        DeleteRow(*it, PivotElement.first);
        // delete column in previous matrix (index of col is index of pivot row)
        CoCoA_ASSERT(it != myResolution.begin());
        std::vector<matrix>::iterator TmpIter(it);
        --TmpIter;
        DeleteCol(*TmpIter, PivotElement.first);
        // delete row in subsequent matrix (index of row is index of pivot column)
        if (it != (--myResolution.end()))
        {
          std::vector<matrix>::iterator TmpIter(it);
          ++TmpIter;
          DeleteRow(*TmpIter, PivotElement.second);
        }
        PivotElement = myFindPivot(*it);
      }
      // if matrix has zero dimension we can skip everything behind
      if (NumRows(*it) == 0 || NumCols(*it) == 0) {
        resEnd = it;
        break;
      }
    }
    myResolution.erase(resEnd, myResolution.end());
  }


  std::vector<RingElem> ResolutionMinimization::myGetColumn(ConstMatrixView m, long col) const
  {
    CoCoA_ASSERT(0 <= col && col < NumCols(m));
    std::vector<RingElem> res;
    const long nrows = NumRows(m);
    res.reserve(nrows);
    for (long i = 0; i < nrows; ++i)
    {
      res.push_back(m(i, col));
    }
    return res;
  }

} // end of namespace CoCoA
