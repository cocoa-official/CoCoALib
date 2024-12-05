#ifndef CoCoA_UPoly_H
#define CoCoA_UPoly_H

//   Copyright (c)  2014  John Abbott, and Anna M. Bigatti

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

#include <vector>
//using std::ostream;

namespace CoCoA
{

  void HornerRecursiveIterQQ2(BigInt& AnsN, const std::vector<BigInt>& c, const BigInt& n, const BigInt& d);

} // end of namespace CoCoA

#endif
