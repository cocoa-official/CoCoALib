#ifndef CoCoA_ToString_H
#define CoCoA_ToString_H

//   Copyright (c)  2011,2014  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/MachineInt.H"

#include <string>
//using std::string;
#include <sstream>
//using std::ostringstream;

namespace CoCoA
{

  class BigInt;      // fwd decl -- defined in BigInt.H
  class BigRat;      // fwd decl -- defined in BigRat.H

  template <typename T>
  std::string ToString(const T& x)
  {
    std::ostringstream buffer;
    buffer << x;
    return buffer.str();
  }

  std::string ScientificStr(const BigInt& N, const MachineInt& SigFig = 5);
  std::string ScientificStr(const BigRat& q, const MachineInt& SigFig = 5);

  std::string FloatStr(const BigInt& N, const MachineInt& SigFig = 5);
  std::string FloatStr(const BigRat& q, const MachineInt& SigFig = 5);

  std::string DecimalStr(const BigInt& N, const MachineInt& DecimalPlaces=3);
  std::string DecimalStr(const BigRat& q, const MachineInt& DecimalPlaces=3);

} // end of namespace CoCoA

#endif
