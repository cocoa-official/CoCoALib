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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/bool3.C,v 1.6 2022/02/18 14:12:01 abbott Exp $
// $Log: bool3.C,v $
// Revision 1.6  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.5  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2016/11/11 14:15:34  abbott
// Summary: Added short-cut to operator<< when ostream is in bad state
//
// Revision 1.3  2012/05/29 07:48:53  abbott
// Implemented simplification of bool3:
//  changed names of the constants,
//  changes names of the testing fns.
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2006/11/02 13:25:43  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.1.1.1  2006/05/30 11:39:37  cocoa
// Imported files
//
// Revision 1.2  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.1  2006/03/12 21:28:33  cocoa
// Major check in after many changes
//
