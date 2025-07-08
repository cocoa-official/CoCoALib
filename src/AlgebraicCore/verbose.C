//   Copyright (c)  2016  John Abbott and Anna M. Bigatti

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


#include "CoCoA/verbose.H"
#include "CoCoA/error.H"
#include "CoCoA/LogStream.H"

#include <algorithm>
using std::min;
#include <fstream>
using std::ofstream;
#include <iostream>
//#include <string>
using std::string;

namespace // anonymous
{
  
  // Global "bit-bucket" ostream
  std::ofstream DevNull;

} // end of namespace anonymous

namespace CoCoA
{

  // static data members (effectively global variables)
  long VerboseLog::ourNestingDepth = 0;
  long VerboseLog::ourVerbosityLevel = 0; // default 0 ==> print nothing
  

  VerboseLog::VerboseLog(const char* const FnName):
      myFnName(FnName)
  {
    ++ourNestingDepth;
  }

  VerboseLog::~VerboseLog()
  {
    --ourNestingDepth;
  }

  std::ostream& VerboseLog::operator()(long level)
  {
    if (level < 1)  CoCoA_THROW_ERROR1(ERR::ReqPositive);
    if (level > ourVerbosityLevel)  return DevNull;
    LogStream() << string(std::min(100L,ourNestingDepth)-1,' '); // indent by ourNestingDepth-1
// min indent 3 spaces (2spaces for CoCoA5 verbosity)
    LogStream() << "   [D" << ourNestingDepth << ",L" << level << "," << myFnName << "]  ";
    return LogStream();
  }


  //  long SetVerbosityLevel(long NewLevel)
  void SetVerbosityLevel(long NewLevel)
  {
    if (NewLevel < 0)  CoCoA_THROW_ERROR1(ERR::ReqNonNegative);
    //    const long OldValue = VerboseLog::ourVerbosityLevel;
    VerboseLog::ourVerbosityLevel = NewLevel;
    //    return OldValue;
  }

  long VerbosityLevel() noexcept
  {
    return VerboseLog::ourVerbosityLevel;
  }

} // end of namespace CoCoA
