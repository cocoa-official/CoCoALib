#ifndef CoCoA_toric_H
#define CoCoA_toric_H

//   Copyright (c)  2011  Anna Bigatti

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

  class ideal;        // forward declaration -- defined in ideal.H

  ideal SequentialToric_C(const ideal& I, const std::vector<long>& indices);
  ideal SequentialToric_C(const SparsePolyRing& P, ConstMatrixView M);
  void EndToric_C();
  

} // end of namespace CoCoA

#endif
