#ifndef CoCoA_DenseMatrix_H
#define CoCoA_DenseMatrix_H

//   Copyright (c)  2005-2009,2011  John Abbott,  Anna M. Bigatti

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


#include <vector>
//using std::vector;


namespace CoCoA
{

  class BigInt;           // fwd decl -- defined in BigInt.H
  class BigRat;           // fwd decl -- defined in BigRat.H
  class ConstMatrixView;  // fwd decl -- defined in matrix.H
  class RingElem;         // fwd decl -- defined in ring.H
  class matrix;           // fwd decl -- defined in matrix.H
  class ring;             // fwd decl -- defined in ring.H

  matrix NewDenseMat(const ring& R, long NumRows, long NumCols);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<long> >& M);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigInt> >& M);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<BigRat> >& M);
  matrix NewDenseMat(const ring& R, const std::vector< std::vector<RingElem> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<long> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<BigInt> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<BigRat> >& M);
  matrix NewDenseMatTranspose(const ring& R, const std::vector< std::vector<RingElem> >& M);
  matrix NewDenseMat(const ConstMatrixView& M);

} // end of namespace CoCoA

#endif
