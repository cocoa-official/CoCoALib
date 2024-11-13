//   Copyright (c)  2017  John Abbott and Anna M. Bigatti

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


#include "CoCoA/LogStream.H"

#include <iostream>

namespace CoCoA
{

  std::ostream* LogStreamForThisBlock::ourActiveLogStreamPtr = &std::cout;

  
  std::ostream& operator<<(std::ostream& out, const LogStreamForThisBlock& /*dummy*/)
  {
    if (!out) return out;  // short-cut for bad ostreams (makes sense)
    return out << "LogStreamForThisBlock";
  }

} // end of namespace CoCoA
