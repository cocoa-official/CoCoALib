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

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/RingElemInput.H,v 1.8 2022/02/18 14:11:48 abbott Exp $
// $Log: RingElemInput.H,v $
// Revision 1.8  2022/02/18 14:11:48  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.7  2020/10/30 19:32:44  abbott
// Summary: Added RingElemVec (redmine 1509)
//
// Revision 1.6  2020/01/10 15:39:33  abbott
// Summary: Corrected copyright notice
//
// Revision 1.5  2019/10/18 14:11:19  bigatti
// -- Renamed ReadExprVector --> RingElems
//
// Revision 1.4  2019/10/09 16:37:36  bigatti
// -- added ReadExprVector
//
// Revision 1.3  2014/07/15 16:54:35  abbott
// Summary: Corrected name of include guard
// Author: JAA
//
// Revision 1.2  2014/03/21 15:57:19  bigatti
// -- Ring is now first argument of ReadExpr
//
// Revision 1.1  2014/01/30 15:15:01  bigatti
// -- first import
//

#endif
