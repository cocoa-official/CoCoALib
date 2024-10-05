#ifndef CoCoA_ExternalLibs_H
#define CoCoA_ExternalLibs_H

//   Copyright (c)  2017  John Abbott, Anna M. Bigatti
//   Author: 2017  John Abbott

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
#include <vector>
//using std::vector;

namespace CoCoA
{

  struct ExternalLibInfo
  {
  public:
    ExternalLibInfo(const std::string& LibName, const std::string& LibFile, const std::string& LibVersion, const std::string& LibWebsite);
  public: // data members
    std::string myName;
    std::string myFile;
    std::string myVersion;
    std::string myWebsite;
  };

  std::ostream& operator<<(std::ostream& out, const ExternalLibInfo& info);

  const std::vector<ExternalLibInfo>& ExternalLibs();


} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/ExternalLibs.H,v 1.3 2022/02/18 14:11:45 abbott Exp $
// $Log: ExternalLibs.H,v $
// Revision 1.3  2022/02/18 14:11:45  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.2  2017/04/27 16:20:11  bigatti
// -- separated name and lib file
//
// Revision 1.1  2017/04/26 14:52:37  abbott
// Summary: Fns to find out about external libs in CoCoALib
//
//

#endif
