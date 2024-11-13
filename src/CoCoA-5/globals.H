#ifndef CoCoA5_globals_H
#define CoCoA5_globals_H

//   Copyright (c)  2018  John Abbott, Anna M. Bigatti

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


#include <fstream>
// using std::ofstream;

namespace CoCoA
{

  extern bool GlobalFlag_SuppressPrompt; // corr to commandline arg --no-prompt

  extern std::ofstream GlobalStatusLogStream; // for status log stream to GUI

  class SystemCommandPermit
  {
  public:
    static void EnableCommand();
    static void DisableCommand();
    static bool IsEnabled();
  private: // data members
    static bool ourGlobalFlag_AllowSysCmd; // false = disabled, true = enabled
  };

} // namespace CoCoA

// Close include guard:
#endif
