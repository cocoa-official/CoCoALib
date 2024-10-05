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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/MatrixSpecial.H,v 1.9 2022/02/18 14:11:46 abbott Exp $
// $Log: MatrixSpecial.H,v $
// Revision 1.9  2022/02/18 14:11:46  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2020/05/26 12:06:18  abbott
// Summary: Renamed TensorMat to KroneckerProd; doc & tests updated
//
// Revision 1.7  2020/02/27 17:39:45  bigatti
// -- added LawrenceMat
//
// Revision 1.6  2019/10/11 19:54:28  abbott
// Summary: Renamed jacobian to JacobianMat
//
// Revision 1.5  2019/10/11 12:54:47  abbott
// Summary: Added SylvesterMat
//
// Revision 1.4  2017/11/08 14:04:39  abbott
// Summary: Added new fn HilbertMat
//
// Revision 1.3  2016/10/27 14:03:04  abbott
// Summary: Added RandomUnimodularMat
//
// Revision 1.2  2011/03/21 07:58:30  bigatti
// -- added TensorMat
// -- changed size into len
//
// Revision 1.1  2011/02/10 15:30:14  bigatti
// -- first import: jacobian
//


#endif // CoCoA_MatrixSpecial_H
