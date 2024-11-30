#ifndef CoCoA_MatrixForOrdering_H
#define CoCoA_MatrixForOrdering_H

//   Copyright (c)  2008,2015  John Abbott,  Anna M. Bigatti

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

  class ConstMatrixView; // fwd decl -- defined in matrix.H
  class ConstMatrix;     // fwd decl -- defined in matrix.H
  class MachineInt;      // fwd decl -- defined in MachineInt.H
  class matrix;          // fwd decl -- defined in matrix.H

  bool IsTermOrdering(const ConstMatrixView& M);
  bool IsPositiveGrading(const ConstMatrixView& M);  ///< is M a positive grading?
  //  bool IsNonNegGrading(const ConstMatrixView& M);    // removed 2019-09

  ConstMatrix LexMat(const MachineInt& n);
  ConstMatrix XelMat(const MachineInt& n);
  ConstMatrix RevLexMat(const MachineInt& n);
  ConstMatrix StdDegRevLexMat(const MachineInt& n);
  ConstMatrix StdDegLexMat(const MachineInt& n);

  matrix MakeTermOrdMat(ConstMatrixView M);
  matrix MakeTermOrdMat(ConstMatrixView M, const MachineInt& GrDim);

  matrix ElimMat(const std::vector<long>& IndetsToElim, const MachineInt& NumIndets);
  matrix ElimMat(const std::vector<long>& IndetsToElim, const ConstMatrixView& GradingM);
  matrix ElimHomogMat(const std::vector<long>& IndetsToElim, const ConstMatrixView& GradingM);

} // end of namespace CoCoA

#endif
