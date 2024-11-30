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
