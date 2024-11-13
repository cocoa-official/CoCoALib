#ifndef CoCoA_verbose_H
#define CoCoA_verbose_H

//   Copyright (c)  2016  John Abbott, and Anna M. Bigatti

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

#include "CoCoA/config.H"

#include <iosfwd>
//using std::ostream;
#include <string>
//using std::string;

namespace CoCoA
{

  class VerboseLog
  {
  public:
    VerboseLog(const char* const FnName);
    ~VerboseLog();
    std::ostream& operator()(long level);
  public: // disable copy ctor & assignment
    VerboseLog(const VerboseLog&) = delete;
    VerboseLog& operator=(const VerboseLog&) = delete;
  private:
    const std::string myFnName;
    static long ourNestingDepth;
    static long ourVerbosityLevel;
    friend void SetVerbosityLevel(long);
    friend long VerbosityLevel() noexcept;
  };

  void SetVerbosityLevel(long NewLevel); // only log messages with level <= NewLevel will be printed
  long VerbosityLevel() noexcept; // return current verbosity level

} // end of namespace CoCoA

#endif
