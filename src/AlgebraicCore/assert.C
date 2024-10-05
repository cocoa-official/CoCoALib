//   Copyright (c)  2005  John Abbott and Anna M. Bigatti

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


#include "CoCoA/assert.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <sstream>
using std::ostringstream;


namespace CoCoA
{

  void AssertionFailed(const char* const cond, const char* const OrigFile, unsigned long OrigLine)
  {
    ostringstream message;
    message << endl
            << "===========================================================================" << endl
            << "===========================================================================" << endl
            << "==== CoCoA Assertion failed: [[" << cond << "]]" << endl
            << "==== File: " << OrigFile << endl
            << "==== Line: " << OrigLine << endl
            << "===========================================================================" << endl
            << "===========================================================================" << endl
            << endl;
    cerr << message.str();

    // Throw a CoCoA error in an unusual way... (so that the "right" file and line no. info appear).
    ThrowException(ErrorInfo(ERR::AssertFail, "AssertionFailed", OrigFile, OrigLine));
  }

} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/assert.C,v 1.7 2022/02/18 14:12:01 abbott Exp $
// $Log: assert.C,v $
// Revision 1.7  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/02/10 19:42:06  abbott
// Summary: Renamed fn params
//
// Revision 1.5  2021/01/07 15:23:38  abbott
// Summary: Corrected copyright
//
// Revision 1.4  2020/06/19 19:39:21  abbott
// Summary: Now all throws go through new template fn ThrowException; seems to work
//
// Revision 1.3  2010/12/17 16:11:06  abbott
// Output related to assertions is now on standard C++ streams (rather
// than GlobalErrput, etc).
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.1  2006/10/06 14:04:57  cocoa
// The new assert header and implementation files.
// A new test.
//
//
