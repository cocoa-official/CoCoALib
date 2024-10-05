#ifndef CoCoA_RandomSparseNonSing01Mat_H
#define CoCoA_RandomSparseNonSing01Mat_H

//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti

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


namespace CoCoA
{

  class matrix;  // fwd decl -- defined in matrix.H
  class MachineInt;  // fwd decl -- define in MachineInt.H
  class ring;  // fwd decl -- define in ring.H

  matrix RandomSparseNonSing01Mat(const ring& R, const MachineInt& N);

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RandomSparseNonSing01Mat.H,v 1.4 2022/02/18 14:11:48 abbott Exp $
// $Log: RandomSparseNonSing01Mat.H,v $
// Revision 1.4  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2020/01/26 14:17:25  abbott
// Summary: Removed useless include
//
// Revision 1.2  2018/01/16 11:42:14  abbott
// Summary: Changed include guard
//
// Revision 1.1  2017/11/14 14:59:35  abbott
// Summary: New fns: RandomSmallPrime, RandomSparseNonSing01Matrix
//
//

#endif
