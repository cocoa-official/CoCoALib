//   Copyright (c) 2014  John Abbott,  Anna M. Bigatti
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

#include "VersionInfo.H"
//-----------------------------------------------------------------
// Version number defined in   src/CoCoA-5/release-files/version
//-----------------------------------------------------------------

// CPP trick to convert a macro value into a string
#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

namespace CoCoA
{

  const std::string& CoCoAVersion()
  {
    static const std::string version(STRINGIFY(COCOA5_VER_MAJ) "."
                                     STRINGIFY(COCOA5_VER_MIN) "."
                                     STRINGIFY(COCOA5_VER_MINMIN));
    return version;
  }


} // end of namespace CoCoA
