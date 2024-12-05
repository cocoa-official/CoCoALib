#ifndef CoCoA_PPMonoidEvOv_H
#define CoCoA_PPMonoidEvOv_H

//   Copyright (c)  2005,2007,2010  John Abbott,  Anna M. Bigatti

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

  class symbol;         // forward declaration -- defined in symbol.H
  class PPMonoid;       // forward declaration -- defined in PPMonoid.H
  class PPOrdering;     // forward declaration -- defined in PPOrdering.H
  class PPOrderingCtor; // forward declaration -- defined in PPOrdering.H

  // Pseudo-ctor
  PPMonoid NewPPMonoidEvOv(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  PPMonoid NewPPMonoidEvOv(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);

} // end of namespace CoCoA

#endif
