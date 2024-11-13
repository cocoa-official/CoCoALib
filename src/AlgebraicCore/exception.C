//   Copyright (c)  2015,2020  John Abbott and Anna M. Bigatti

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


#include "CoCoA/exception.H"
#include "CoCoA/utils.H"

#include <iostream>
#include <sstream>
using std::ostringstream;
//#include <string>
using std::string;

namespace CoCoA
{

  // Must define this because of the virtual fns
  exception::~exception()
  {}
  

  void exception::myOutputSelf(std::ostream& out) const
  {
    if (!out) return;  // short-cut for bad ostreams
    out << "CoCoA::exception(\"" << myMessage << '"';
    if (!myContext.empty()) out << ", context=\"" << myContext << '"';
    out << ")";
  }

  std::ostream& operator<<(std::ostream& out, const exception& exc)
  {
    exc.myOutputSelf(out);
    return out;
  }

  void PrintInFrame(std::ostream& out, const exception& exc)
  {
    if (!out) return;  // short-cut for bad ostreams
    ostringstream buffer;
    buffer << ">>>> " << exc << " <<<<";
    const long n = len(buffer.str());
    const string HorizontalLine(n, '-');
    using std::endl;
    out << endl
        << HorizontalLine << endl
        << buffer.str() << endl
        << HorizontalLine << endl;
  }


} // end of namespace CoCoA
