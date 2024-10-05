//   Copyright (c)  2018  Anna Bigatti, John Abbott
//   Main author: John Abbott

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

#include "globals.H"


namespace CoCoA
{

  bool GlobalFlag_SuppressPrompt = false; // can be set by command line flag --no-prompt


  std::ofstream GlobalStatusLogStream;

  bool SystemCommandPermit::ourGlobalFlag_AllowSysCmd = false;

  void SystemCommandPermit::EnableCommand()
  {
    ourGlobalFlag_AllowSysCmd = true;
  }

  void SystemCommandPermit::DisableCommand()
  {
    ourGlobalFlag_AllowSysCmd = false;
  }

  bool SystemCommandPermit::IsEnabled()
  {
    return ourGlobalFlag_AllowSysCmd;
  }

} // namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/CoCoA-5/globals.C,v 1.5 2022/02/22 20:39:28 abbott Exp $
// $Log: globals.C,v $
// Revision 1.5  2022/02/22 20:39:28  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.4  2021/02/22 19:40:28  abbott
// Summary: Added flag to handle --no-prompt commandline option (redmine 500)
//
// Revision 1.3  2020/03/06 18:41:37  abbott
// Summary: Changed name GlobalStatusOStream to GlobalStatusLogStream; changed name of CLI flag
//
// Revision 1.2  2020/02/21 14:03:38  abbott
// Summary: Implemented status messages (redmine 1399) via a file (supplied on command line)
//
// Revision 1.1  2019/03/04 13:16:27  abbott
// Summary: Added new files globals.H/globals.C
//
//
