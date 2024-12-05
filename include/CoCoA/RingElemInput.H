#ifndef CoCoA_RingElemInput_H
#define CoCoA_RingElemInput_H

//   Copyright (c)  2014,2019  John Abbott & Anna M. Bigatti

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

#include <istream>
// using std::istream;
#include <string>
// using std::string;
#include <vector>
// using std::vector;

namespace CoCoA
{
  class ring;          // fwd decl -- defined in ring.H
  class RingElem;      // fwd decl -- defined in ring.H

  RingElem ReadExpr(const ring& P, std::istream& in);
  RingElem ReadExpr(const ring& P, const std::string& s);

  RingElem ReadExprSemicolon(const ring& P, std::istream& in);
  RingElem ReadExprSemicolon(const ring& P, const std::string& s);

  std::vector<RingElem> RingElems(const ring& P, std::istream& in);
  std::vector<RingElem> RingElems(const ring& P, const std::string& s);
  std::vector<RingElem> RingElemVec(const ring& P, std::istream& in);
  std::vector<RingElem> RingElemVec(const ring& P, const std::string& s);
  

} // end of namespace CoCoA

#endif
