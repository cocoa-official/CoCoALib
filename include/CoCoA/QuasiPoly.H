#ifndef CoCoA_QuasiPoly_H
#define CoCoA_QuasiPoly_H

//   Copyright (c)  2014  John Abbott,  Anna M. Bigatti

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


#include <iosfwd>
//using std::ostream;
#include <vector>
//using std::vector;

namespace CoCoA
{

  class RingElem;   // fwd decl -- defined in ring.H
  class BigInt;     // fwd decl -- defined in BigInt.H
  class MachineInt; // fwd decl -- defined in MachineInt.H

  class QuasiPoly
  {
  public:
    explicit QuasiPoly(const std::vector<RingElem>& v);
    RingElem operator()(const MachineInt& n) const;
    RingElem operator()(const BigInt& N) const;
  private: // data members
    std::vector<RingElem> myConstituents;
    friend const std::vector<RingElem>& constituents(const QuasiPoly& p);
  };

  const std::vector<RingElem>& constituents(const QuasiPoly& p);
  std::ostream& operator<<(std::ostream& out, const QuasiPoly& p);

} // end of namespace CoCoA

#endif
