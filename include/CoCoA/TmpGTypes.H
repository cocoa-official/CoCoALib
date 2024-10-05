#ifndef CoCoA_TmpGTypes_H
#define CoCoA_TmpGTypes_H

//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti
//   Original author: 2007  Massimo Caboara

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
// using vector;

namespace CoCoA
{

  class RingElem;  // fwd declaration -- defined in ring.H
  class ModuleElem;// fwd declaration -- defined in module.H

  // don't use these types: not for public use
  typedef std::vector<RingElem> PolyList;
  typedef std::vector<ModuleElem> VectorList;
}

#endif
