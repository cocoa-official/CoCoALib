#ifndef CoCoA_PPMonoidOv_H
#define CoCoA_PPMonoidOv_H

//   Copyright (c)  2005-2017  John Abbott and Anna M. Bigatti
//   Author:  2005,2007,2010  John Abbott

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


#include "CoCoA/OrdvArith.H"
#include <vector>
// using vector;

namespace CoCoA
{

  class symbol;         // forward declaration -- defined in symbol.H
  class PPMonoid;       // forward declaration -- defined in PPMonoid.H
  class PPOrdering;     // forward declaration -- defined in PPOrdering.H
  class PPOrderingCtor; // forward declaration -- defined in PPOrdering.H

  // Pseudo-ctor
  PPMonoid NewPPMonoidOv(const std::vector<symbol>& IndetNames, const PPOrdering& ord);
  PPMonoid NewPPMonoidOv(const std::vector<symbol>& IndetNames, const PPOrderingCtor& ord);

  bool IsPPMonoidOv(const PPMonoid& PPM);
  const OrdvArith::reference& OrdvA(const PPMonoid& PPM);

} // end of namespace CoCoA

#endif
