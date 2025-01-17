#ifndef CoCoA_MatrixOps_H
#define CoCoA_MatrixOps_H

//   Copyright (c)  2005,2008,2012  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/BigInt.H"
#include "CoCoA/matrix.H"

#include <vector>
// using std::vector

namespace CoCoA
{

  class BigInt;          // fwd decl -- defined in BigInt.H
  class BigRat;          // fwd decl -- defined in BigRat.H
  class ConstMatrixView; // fwd decl -- defined in matrix.H
  class RingElemAlias   ;// fwd decl -- defined in ring.H
  typedef const RingElemAlias& ConstRefRingElem;
  class MachineInt;      // fwd decl -- defined in MachineInt.H
  class MatrixView;      // fwd decl -- defined in matrix.H
  class RingElem;        // fwd decl -- defined in ring.H
  class matrix;          // fwd decl -- defined in matrix.H
  class ring;            // fwd decl -- defined in ring.H

  std::vector<RingElem> GetRow(ConstMatrixView M, long i);
  std::vector<RingElem> GetCol(ConstMatrixView M, long j);
  std::vector< std::vector<RingElem> > GetRows(ConstMatrixView M);
  std::vector< std::vector<RingElem> > GetCols(ConstMatrixView M);

  std::vector<RingElem> FlattenByRows(ConstMatrixView M);
  std::vector<RingElem> FlattenByCols(ConstMatrixView M);


  matrix operator*(ConstMatrixView Mleft, ConstMatrixView Mright);
  matrix operator+(ConstMatrixView Mleft, ConstMatrixView Mright);
  matrix operator-(ConstMatrixView Mleft, ConstMatrixView Mright);
/////  void mul(matrix& lhs, ConstMatrixView M1, ConstMatrixView M2);
  matrix power(ConstMatrixView M, long n);
  matrix power(ConstMatrixView M, const BigInt& N);

  matrix operator*(ConstRefRingElem x, ConstMatrixView M);
  matrix operator*(const BigRat& x, ConstMatrixView M);
  matrix operator*(const BigInt& x, ConstMatrixView M);
  matrix operator*(const MachineInt& x, ConstMatrixView M);
  matrix operator*(ConstMatrixView M, ConstRefRingElem x);
  matrix operator*(ConstMatrixView M, const BigRat& x);
  matrix operator*(ConstMatrixView M, const BigInt& x);
  matrix operator*(ConstMatrixView M, const MachineInt& x);
  matrix operator-(const ConstMatrixView& M);

  matrix operator/(ConstMatrixView M, ConstRefRingElem x);
  matrix operator/(ConstMatrixView M, const BigRat& x);
  matrix operator/(ConstMatrixView M, const BigInt& x);
  matrix operator/(ConstMatrixView M, const MachineInt& x);

  inline matrix& operator+=(matrix& Mlhs, ConstMatrixView Mrhs) { Mlhs = Mlhs+Mrhs; return Mlhs; } // assignment does not copy because matrix is ref-counted
  inline matrix& operator-=(matrix& Mlhs, ConstMatrixView Mrhs) { Mlhs = Mlhs-Mrhs; return Mlhs; } //
  inline matrix& operator*=(matrix& Mlhs, ConstMatrixView Mrhs) { Mlhs = Mlhs*Mrhs; return Mlhs; } //
  inline matrix& operator*=(matrix& Mlhs, ConstRefRingElem r) { Mlhs = Mlhs*r; return Mlhs; }      //
  inline matrix& operator*=(matrix& Mlhs, const BigRat& q) { Mlhs = Mlhs*q; return Mlhs; }         //
  inline matrix& operator*=(matrix& Mlhs, const BigInt& N) { Mlhs = Mlhs*N; return Mlhs; }         //
  inline matrix& operator*=(matrix& Mlhs, const MachineInt& n) { Mlhs = Mlhs*n; return Mlhs; }     //
  inline matrix& operator/=(matrix& Mlhs, ConstRefRingElem r) { Mlhs = Mlhs/r; return Mlhs; }      //
  inline matrix& operator/=(matrix& Mlhs, const BigRat& q) { Mlhs = Mlhs/q; return Mlhs; }         //
  inline matrix& operator/=(matrix& Mlhs, const BigInt& N) { Mlhs = Mlhs/N; return Mlhs; }         //
  inline matrix& operator/=(matrix& Mlhs, const MachineInt& n) { Mlhs = Mlhs/n; return Mlhs; }     //

  RingElem FrobeniusNormSq(ConstMatrixView A);
  RingElem OperatorNormInfinity(ConstMatrixView M);
  RingElem OperatorNorm1(ConstMatrixView M);

  RingElem det(const ConstMatrixView& M);
  RingElem Det_OLD(const ConstMatrixView& M);  // !!!BUG!!! JUST EXPERIMENTAL
  bool IsZeroDet(const ConstMatrixView& M);
  long rk(const ConstMatrixView& M);
  matrix inverse(ConstMatrixView M);
  matrix adj(ConstMatrixView M);
  matrix PseudoInverse(ConstMatrixView M);
  matrix LinSolve(ConstMatrixView M, ConstMatrixView rhs);
  matrix LinSolveByGauss(ConstMatrixView M, ConstMatrixView rhs);
  matrix LinSolveByHNF(ConstMatrixView M, ConstMatrixView rhs);
  matrix LinSolveByModuleRepr(ConstMatrixView M, ConstMatrixView rhs);
  matrix LinKer(ConstMatrixView M);   // result is matrix over QQ such that M*K =0
  matrix LinKerZZ(ConstMatrixView M); // result is matrix over ZZ such that M*K=0
  matrix LinKerByGauss(ConstMatrixView M);
  matrix rref(ConstMatrixView M);
  matrix RREFByGauss(ConstMatrixView M);

  RingElem DetOfSmallMat(ConstMatrixView M);
  RingElem det2x2(ConstMatrixView M);
  RingElem det3x3(ConstMatrixView M);
  RingElem det4x4(ConstMatrixView M);
  RingElem det5x5(ConstMatrixView M);
  RingElem DetDirect(ConstMatrixView M);
  RingElem DetByGauss(ConstMatrixView M);
  RingElem DetByBareiss(const ConstMatrixView& M);
  RingElem DetByMinors(const ConstMatrixView& M);
  RingElem DetByCRT(const ConstMatrixView& M);
  RingElem DetOverQQ(const ConstMatrixView& M);

  long RankByGauss(std::vector<long>& IndepRows, ConstMatrixView M);
  matrix InverseByGauss(ConstMatrixView M);
  matrix AdjDirect(ConstMatrixView M);
  matrix AdjByDetOfMinors(ConstMatrixView M);
  matrix AdjByInverse(ConstMatrixView M); // base ring must be integral domain

  matrix rref(ConstMatrixView M);

  matrix KroneckerProd(ConstMatrixView M1, ConstMatrixView M2);

  struct HadamardRowCol
  {
    HadamardRowCol(ConstRefRingElem R, ConstRefRingElem C): myRowBoundSq(R), myColBoundSq(C) {}
    RingElem myRowBoundSq;
    RingElem myColBoundSq;
  };
  HadamardRowCol HadamardBoundSq(const ConstMatrixView& M);


  bool IsZero(const ConstMatrixView& M);
  bool IsZeroRow(const ConstMatrixView& M, long i);
  bool IsZeroCol(const ConstMatrixView& M, long j);
  bool IsSymmetric(const ConstMatrixView& M);
  bool IsAntiSymmetric(const ConstMatrixView& M);
  bool IsDiagonal(const ConstMatrixView& M);
  bool IsLowerTriangular(const ConstMatrixView& M);
  bool IsUpperTriangular(const ConstMatrixView& M);
  bool IsMat0x0(const ConstMatrixView& M);
  bool HasNegEntry(const ConstMatrixView& M);

//   void GramSchmidtRows(MatrixView& M);
//   void GramSchmidtRows(MatrixView& M, long row);

} // end of namespace CoCoA

#endif
