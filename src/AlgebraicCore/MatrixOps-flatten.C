//   Copyright (c)  2020  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/MatrixOps.H"
#include "CoCoA/matrix.H"
#include "CoCoA/MachineInt.H"

//#include <vector>
using std::vector;

namespace CoCoA
{

  std::vector<RingElem> FlattenByRows(ConstMatrixView M)
  {
    const long r = NumRows(M);
    const long c = NumCols(M);
    const long NumEntries = r*c; // BUG??? could overflow?
    if (NumEntries == 0) return vector<RingElem>();
    vector<RingElem> v;  v.reserve(NumEntries);
    for (long i=0; i < r; ++i)
      for (long j=0; j < c; ++j)
        v.push_back(M(i,j));
    return v;
  }

  std::vector<RingElem> FlattenByCols(ConstMatrixView M)
  {
    const long r = NumRows(M);
    const long c = NumCols(M);
    const long NumEntries = r*c; // BUG??? could overflow?
    if (NumEntries == 0) return vector<RingElem>();
    vector<RingElem> v;  v.reserve(NumEntries);
    for (long j=0; j < c; ++j)
      for (long i=0; i < r; ++i)
        v.push_back(M(i,j));
    return v;
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/MatrixOps-flatten.C,v 1.2 2022/02/18 14:11:55 abbott Exp $
// $Log: MatrixOps-flatten.C,v $
// Revision 1.2  2022/02/18 14:11:55  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.1  2020/10/02 18:40:37  abbott
// Summary: Added new fns FlattenByRows & FlattenByCols
//
//
