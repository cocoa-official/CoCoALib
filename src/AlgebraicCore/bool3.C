//   Copyright (c)  2006,2012  John Abbott and Anna M. Bigatti

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


#include "CoCoA/bool3.H"
#include "CoCoA/OpenMath.H"

#include <iostream>
using std::ostream;

namespace CoCoA
{

  std::ostream& operator<<(std::ostream& out, bool3 flag)
  {
    if (!out) return out;  // short-cut for bad ostreams
    if (IsFalse3(flag)) return out << "false3";
    if (IsTrue3(flag)) return out << "true3";
    return out << "uncertain3";
  }


  OpenMathOutput& operator<<(OpenMathOutput& out, bool3 flag)
  {
    if (IsFalse3(flag)) return out << OpenMathSymbol("logic3", "false");
    if (IsTrue3(flag)) return out << OpenMathSymbol("logic3", "true");
    return out << OpenMathSymbol("logic3", "uncertain");
  }

} // end of namespace CoCoA
