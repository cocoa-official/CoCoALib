//  Copyright (c)  2012  John Abbott,  Anna M Bigatti

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


#include <iostream>
#include <string>
#include "OnlineHelp.H"

// This global is needed by something
std::string packageDir = "packages";

int main()
{
  try
  {
    CoCoA::OnlineHelp::PrintAllExamplesWithoutOutput(std::cout);
  }
  catch (...)//(const CoCoA::ErrorInfo& err)
  {
//    ANNOUNCE(std::cerr, err);
    std::cerr << "---------------------------------------------------------------------------------" << std::endl;
    std::cerr << "Failed to open CoCoA manual XML file; it ought to be in CoCoAManual/CoCoAHelp.xml" << std::endl;
    std::cerr << "---------------------------------------------------------------------------------" << std::endl;
    exit(1);
  }
}
