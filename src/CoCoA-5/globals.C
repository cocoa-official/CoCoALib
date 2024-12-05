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
