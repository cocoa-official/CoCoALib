#ifndef CoCoA_MatrixSpecial_H
#define CoCoA_MatrixSpecial_H

//   Copyright (c)  2011,2016 Anna Bigatti

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

#include "CoCoA/ring.H"
// need only fwd decl of ring. RingElem, ConstRefRingElem (but last is a typedef)

#include <vector>
// using std::vector

namespace CoCoA
{

  class MachineInt;      // forward declaration -- defined in MachineInt.H
  class matrix;          // forward declaration -- defined in matrix.H
////  class ring;            // forward declaration -- defined in ring.H
////  class RingElem;        // forward declaration -- defined in ring.H
////  class ConstRefRingElem; // forward declaration -- defined in ring.H

  matrix JacobianMat(const std::vector<RingElem>& polys);
  matrix JacobianMat(const std::vector<RingElem>& polys, const std::vector<RingElem>& inds);

  matrix LawrenceMat(ConstMatrixView M);

  matrix SylvesterMat(ConstRefRingElem f, ConstRefRingElem g, ConstRefRingElem x);

  matrix RandomUnimodularMat(const ring& R, const MachineInt& N, const MachineInt& Niters=MachineInt(0));

  matrix HilbertMat(const MachineInt& N);

} // end of namespace CoCoA

#endif // CoCoA_MatrixSpecial_H
