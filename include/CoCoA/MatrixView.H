#ifndef CoCoA_MatrixView_H
#define CoCoA_MatrixView_H

//   Copyright (c)  2005,2008,2014  John Abbott,  Anna M. Bigatti

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
// using std::vector;


namespace CoCoA
{

  class BigInt;           // fwd decl -- defined in BigInt.H
  class BigRat;           // fwd decl -- defined in BigRat.H
  class ConstMatrixView;  // fwd decl -- defined in matrix.H
  class RingElemAlias;    // fwd decl -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;
  class MachineInt;       // fwd decl -- defined in MachineInt.H
  class MatrixView;       // fwd decl -- defined in matrix.H
  class RingElem;         // fwd decl -- defined in ring.H
  class ring;             // fwd decl -- defined in ring.H


  MatrixView transpose(MatrixView M);
  ConstMatrixView transpose(ConstMatrixView M);

  MatrixView submat(MatrixView M, const std::vector<long>& rows, const std::vector<long>& cols);
  ConstMatrixView submat(ConstMatrixView M, const std::vector<long>& rows, const std::vector<long>& cols);

  MatrixView FirstCols(MatrixView M, long ncols);
  ConstMatrixView FirstCols(ConstMatrixView M, long ncols);
  MatrixView FirstRows(MatrixView M, long nrows);
  ConstMatrixView FirstRows(ConstMatrixView M, long nrows);
  MatrixView ColMat(MatrixView M, long col);
  ConstMatrixView ColMat(ConstMatrixView M, long col);
  MatrixView RowMat(MatrixView M, long row);
  ConstMatrixView RowMat(ConstMatrixView M, long row);

  MatrixView ColMat(std::vector<RingElem>& v);
  ConstMatrixView ColMat(const std::vector<RingElem>& v);
  MatrixView RowMat(std::vector<RingElem>& v);
  ConstMatrixView RowMat(const std::vector<RingElem>& v);
  MatrixView DiagMat(std::vector<RingElem>& DiagEntries);
  ConstMatrixView DiagMat(const std::vector<RingElem>& DiagEntries);

  MatrixView MatByRows(const MachineInt& nrows, const MachineInt& ncols, std::vector<RingElem>& v);
  ConstMatrixView MatByRows(const MachineInt& nrows, const MachineInt& ncols, const std::vector<RingElem>& v);
  MatrixView MatByCols(const MachineInt& nrows, const MachineInt& ncols, std::vector<RingElem>& v);
  ConstMatrixView MatByCols(const MachineInt& nrows, const MachineInt& ncols, const std::vector<RingElem>& v);

//   MatrixView MatByRows(std::vector< std::vector<RingElem> >& vv);
//   ConstMatrixView MatByRows(const std::vector< std::vector<RingElem> >& vv);
//   MatrixView MatByCols(std::vector< std::vector<RingElem> >& vv);
//   ConstMatrixView MatByCols(const std::vector< std::vector<RingElem> >& vv);

  /***************************************************************************/
  // BlockMat2x2 -- block matrices  ( A  B )
  //                                ( C  D )
  enum class ZeroBlockIndicator { zeroes };
  constexpr ZeroBlockIndicator zeroes = ZeroBlockIndicator::zeroes;
  MatrixView BlockMat2x2(MatrixView A, MatrixView B, MatrixView C, MatrixView D);
  MatrixView BlockMat2x2(MatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, MatrixView D);
  MatrixView BlockMat2x2(ZeroBlockIndicator, MatrixView B, MatrixView C,  ZeroBlockIndicator);
  ConstMatrixView BlockMat2x2(ConstMatrixView A, ConstMatrixView B, ConstMatrixView C, ConstMatrixView D);
  ConstMatrixView BlockMat2x2(ConstMatrixView A, ZeroBlockIndicator, ZeroBlockIndicator, ConstMatrixView D);
  ConstMatrixView BlockMat2x2(ZeroBlockIndicator, ConstMatrixView B, ConstMatrixView C,  ZeroBlockIndicator);

  MatrixView ConcatVer(MatrixView M1, MatrixView M2);
  ConstMatrixView ConcatVer(ConstMatrixView M1, ConstMatrixView M2);
  MatrixView ConcatHor(MatrixView M1, MatrixView M2);
  ConstMatrixView ConcatHor(ConstMatrixView M1, ConstMatrixView M2);
  MatrixView ConcatDiag(MatrixView M1, MatrixView M2);
  ConstMatrixView ConcatDiag(ConstMatrixView M1, ConstMatrixView M2);
  MatrixView ConcatAntiDiag(MatrixView M1, MatrixView M2);
  ConstMatrixView ConcatAntiDiag(ConstMatrixView M1, ConstMatrixView M2);

} // end of namespace CoCoA

#endif
