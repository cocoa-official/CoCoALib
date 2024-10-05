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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/MatrixOps.H,v 1.17 2023/12/01 20:21:47 abbott Exp $
// $Log: MatrixOps.H,v $
// Revision 1.17  2023/12/01 20:21:47  abbott
// Summary: Renamed fields in record produces by HadamardBoundSq
//
// Revision 1.16  2022/06/08 15:13:45  abbott
// Summary: Added comment (for assign-op operators)
//
// Revision 1.15  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.14  2021/06/21 11:19:56  abbott
// Summary: Removed unsightly space
//
// Revision 1.13  2021/01/07 14:16:30  abbott
// Summary: Corrected copyright
//
// Revision 1.12  2020/10/02 18:40:37  abbott
// Summary: Added new fns FlattenByRows & FlattenByCols
//
// Revision 1.11  2020/09/28 11:15:58  abbott
// Summary: Added DetByMinors; improved DetDirect so it works in all rings
//
// Revision 1.10  2020/09/22 18:13:53  abbott
// Summary: Renamed FrobeniusNorm2 to FrobeniusNormSq; for det stuff into new file MatrixOps-det.C
//
// Revision 1.9  2020/06/26 15:03:01  abbott
// Summary: Added defns for += & similar ops (all inline)
//
// Revision 1.8  2020/06/22 15:42:27  abbott
// Summary: Commented out mul(,,)
//
// Revision 1.7  2020/05/26 12:06:18  abbott
// Summary: Renamed TensorMat to KroneckerProd; doc & tests updated
//
// Revision 1.6  2020/02/28 08:55:53  abbott
// Summary: Added LinKerZZ
//
// Revision 1.5  2020/02/27 10:53:47  abbott
// Summary: Added new fn IsZeroDet
//
// Revision 1.4  2020/02/12 09:14:46  abbott
// Summary: Renamed RowReducedForm to rref
//
// Revision 1.3  2020/01/09 18:35:41  abbott
// Summary: Added rref (Row reduced echelon form)
//
// Revision 1.2  2019/09/16 17:27:19  abbott
// Summary: Added GetRow, GetRows, GetCol, GetCols, DetOfSmallMat, RowReducedForm
//
// Revision 1.1  2018/05/17 15:25:53  bigatti
// -- renamed MatrixOperations --> MatrixOps
//
// Revision 1.8  2018/02/27 11:04:08  abbott
// Summary: Added DetOverQQ
//
// Revision 1.7  2018/02/15 14:35:42  abbott
// Summary: Added det4x4 and det5x5, DetByCRT
//
// Revision 1.6  2015/12/11 15:51:23  bigatti
// -- added IsLowerTriangular
//
// Revision 1.5  2015/12/11 13:03:20  bigatti
// -- added IsUpperTriangular, IhaveNegEntry
// -- removed useless checks
//
// Revision 1.4  2015/04/13 15:39:16  abbott
// Summary: Changed "rank" --> "rk"; adjoint" --> "adj"; added "DetDirect"
// Author: JAA
//
// Revision 1.3  2014/08/26 12:56:14  abbott
// Summary: Added DetByBareiss
// Author: JAA
//
// Revision 1.2  2014/04/17 14:09:46  bigatti
// -- moved Is... functions from matrix.H to here
//
// Revision 1.1  2014/04/11 15:42:37  abbott
// Summary: Renamed from MatrixArith
// Author: JAA
//
// Revision 1.23  2012/10/24 11:42:18  abbott
// Revised fwd decl of ConstRefRingElem.
//
// Revision 1.22  2012/10/16 09:45:45  abbott
// Replaced RefRingElem by RingElem&.
//
// Revision 1.21  2012/07/10 12:00:47  bigatti
// -- added LinKer, LinKerByGauss
//
// Revision 1.20  2012/06/19 15:43:27  abbott
// Added division of a matrix by a scalar.
//
// Revision 1.19  2012/06/11 08:20:33  abbott
// Added multiplication on the right by a scalar.
//
// Revision 1.18  2012/06/10 22:57:31  abbott
// Added negation for matrices -- same as doing (-1)*M.
//
// Revision 1.17  2012/04/27 14:49:34  abbott
// Added LinSolve family (incl. LinSolveByGauss, LinSolveByHNF, LinSolveByModuleRepr).
//
// Revision 1.16  2012/04/13 16:24:35  abbott
// Added solve and SolveByGauss.
//
// Revision 1.15  2011/11/09 13:47:56  bigatti
// -- renamed MachineInteger --> MachineInt
//
// Revision 1.14  2011/08/24 10:20:28  bigatti
// -- renamed QQ --> BigRat
//
// Revision 1.13  2011/08/14 15:52:18  abbott
// Changed ZZ into BigInt (phase 1: just the library sources).
//
// Revision 1.12  2011/05/13 16:47:20  abbott
// Added power fn for matrices: partial impl, cannot yet handle negative powers.
//
// Revision 1.11  2011/03/03 13:50:22  abbott
// Replaced several occurrences of std::size_t by long; there's still more
// work to do though!
//
// Revision 1.10  2011/03/01 14:13:24  bigatti
// -- added f*M
//
// Revision 1.9  2011/02/28 14:08:49  bigatti
// -- added det3x3
// -- using apply mapping matrix (in DetByGauss)
//
// Revision 1.8  2011/02/18 12:54:59  bigatti
// -- commented out GrammSchmidtRows
//
// Revision 1.7  2011/02/09 16:48:27  bigatti
// -- added + and - for matrices
//
// Revision 1.6  2008/04/21 11:23:11  abbott
// Separated functions dealing with matrices and PPOrderings into a new file.
// Added matrix norms, and completed adjoint.
//
// Revision 1.5  2008/04/18 15:35:57  abbott
// (long overdue) Major revision to matrices
//
// Revision 1.4  2008/04/16 17:24:17  abbott
// Further cleaning of the new matrix code.  Updated documentation too.
//
// Revision 1.3  2008/04/08 15:26:42  abbott
// Major revision to matrix implementation: added matrix views.
// Lots of changes.
//
// Revision 1.2  2007/10/30 17:14:12  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.8  2007/03/08 18:42:06  cocoa
// Cleaned up whitespace.
//
// Revision 1.7  2007/03/07 13:07:10  cocoa
// -- minor: added comment for NewMatrixMinimize
//
// Revision 1.6  2006/12/06 17:14:12  cocoa
// -- removed #include "config.H"
//
// Revision 1.5  2006/11/24 17:27:46  cocoa
// -- reorganized includes of header files
//
// Revision 1.4  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.3  2006/08/17 09:39:08  cocoa
// -- added: elimination ordering matrix for non-homogeneous input
//
// Revision 1.2  2006/07/17 16:58:05  cocoa
// -- added: NewMatrixElim(size_t NumIndets, std::vector<size_t> IndetsToElim)
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.4  2006/04/10 13:20:43  cocoa
// -- fixed buglets for Elimination orderings
//
// Revision 1.3  2006/04/05 17:18:33  cocoa
// -- added IsPositiveGrading
//
// Revision 1.2  2006/03/27 12:21:26  cocoa
// Minor silly changes to reduce number of complaints from some compiler or other.
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.4  2005/04/19 15:39:55  cocoa
// Matrices now use reference counts.
//
// Revision 1.3  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.2  2005/03/30 17:15:14  cocoa
// Cleaned the SpecialMatrix code; a simple test compiles and
// runs fine.
//
// Revision 1.1  2005/03/11 18:40:54  cocoa
// -- first import
//

#endif
