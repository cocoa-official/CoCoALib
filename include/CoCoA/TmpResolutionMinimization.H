#ifndef CoCoA_TmpResolutionMinimization_H
#define CoCoA_TmpResolutionMinimization_H

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

#include "CoCoA/SparsePolyRing.H"

namespace CoCoA
{

  class matrix; // fwd decl -- defined in matrix.H

  /**
   * This class takes a non-minimal free resolution and converts it to a minimal free resolution.
   * The resolution must be represented by a vector of matrices.
   */
  class ResolutionMinimization
  {
  public:
    /**
     * Constructor: It takes the base ring and a vector of matrices
     * which represents the free resolution and initialises the class.
     * The first map of the resolution must be in the first component of
     * the vector of matrices.
     *
     * @param ring
     * @param matrices
     */
    ResolutionMinimization(SparsePolyRing ring, const std::vector<matrix>& matrices):
        myRing(ring),
        myResolution(matrices)
    {}

    /**
     * Reduces the given free resolution to a minimal one.
     */
    void myMinimization();

    /**
     * Returns complete resolution as vector<matrix>
     *
     * @return a vector of matrices which represents a free resolution
     */
    inline const std::vector<matrix>& myGetResolution() const
    {
      return myResolution;
    }

    /**
     * Returns map from F_i -> F_(i-1) as matrix (starting with zero).
     * No check that index is valid!!!
     *
     * @param i
     *
     * @return the map F_i as matrix
     */
    inline const matrix& myGetResolution(long i) const
    {
      return myResolution[i];
    }

  private: // impl details
    /**
     * Finds a pivot element for the reduction process.
     * A pivot element is an invertible entry in the matrix m.
     * If there is no pivot element it returns (-1, -1).
     *
     * @param m
     *
     * @return a pair of longs which specifies the position of the pivot element.
     */
    std::pair<long, long> myFindPivot(ConstMatrixView m) const;

    /*
     * Manipulating the matrix m except row r:
     * Subtract -(PivotColumn[i]/PivotColumn[r]) from i-th row
     * Attention: m is a reference!
     *
     * @param m the matrix which should be manipulated
     * @param r the index of the PivotRow
     * @param PivotColumn the PivotColumn
     */
    void myManipulateMatrix(matrix& m, long r, const std::vector<RingElem>& PivotColumn) const;

    /**
     * Returns the col-th column of m.
     * @param m
     * @param col
     *
     * @return a vector of RingElems which represents the col-th column of m
     */
    std::vector<RingElem> myGetColumn(ConstMatrixView m, long col) const;

    SparsePolyRing myRing;

    std::vector<matrix> myResolution;
  };
}
#endif
