//   Copyright (c) 2012  Anna M. Bigatti,  John Abbott
//
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

#include "Banner.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/ExternalLibs.H"
#include "CoCoA/utils.H"

#include <sstream>
#include <vector>
using std::vector;

using std::string;

namespace CoCoA
{

  namespace // anonymous
  {

    string ExtLibNames()
    {
      const vector<ExternalLibInfo>& LibInfo = ExternalLibs();
      string ExtLibs;

      const int NumLibs = len(LibInfo);

      for (int i=0; i < NumLibs; ++i)
      {
        if (!ExtLibs.empty()) ExtLibs += ", ";
        ExtLibs += LibInfo[i].myName;
//        ExtLibs += LibInfo[i].myName + "(v." + LibInfo[i].myVersion + ")";
      }
      return ExtLibs;
    }

  } // end of anonymous namespace


  std::string CoCoA5Banner()
  {
    const string ExtLibs = ExtLibNames();

    string banner(
      "   ______      ______      ___         ______\n"
      "  / ____/___  / ____/___  /   |       / ____/\n"
      " / /   / __ `/ /   / __ `/ /| |______/___ `  \n"
      "/ /___/ /_/ / /___/ /_/ / ___ /_____/___/ /  \n"
      "`____/`____/`____/`____/_/  |_|    /_____/   \n\n"
      "With CoCoALib");

    if ( !ExtLibs.empty() ) banner += " and external libraries " + ExtLibs;
    banner += string("\n"
                     "indent(VersionInfo(), 2); -- for information about this version");
    // Initial ring is defined in Interpreter.C (near line 1229 setTopLevelVar("P",...))
    banner += string("\n"
                     "Initial CurrentRing is  P ::= QQ[x,y,z];");
    return banner;
  }


  std::string CoCoA5BannerNonFixedWidthFonts()
  {
    const string ExtLibs = ExtLibNames();
    string banner("\n"
                  "/******      C o C o A - 5      ******/\n");
    if ( !ExtLibs.empty() ) banner += "\nWith external libraries: " + ExtLibs;
    banner += string("\n"
                     "indent(VersionInfo()); -- for information about this version");
    return banner;
  }

} // end of namespace CoCoA
