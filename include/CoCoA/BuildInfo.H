#ifndef CoCoA_BuildInfo_H
#define CoCoA_BuildInfo_H

//   Copyright (c)  2005,2007  John Abbott,  Anna M. Bigatti

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


#include<iosfwd>
//using std::ostream;
#include <string>

namespace CoCoA
{
  // String constants containing the build information.
  namespace BuildInfo
  {
    // These are probably NOT noexcept because the 1st call creates a std::string object.
    const std::string& version();       ///< Which version of CoCoALib is this?
    const std::string& GitRevision();   ///< Which Git revision of CoCoALib is this?
    const std::string& compiler();
    const std::string& CompilationFlags();
    const std::string& CompilationPlatform();
    const std::string& CompilationPreprocessorDefines();
    std::string CoCoA_32BIT_LONG_Flag();

    void PrintAll(std::ostream& out);       ///< Print all the build information on stream out.
  }

} // end of namespace CoCoA

#endif
