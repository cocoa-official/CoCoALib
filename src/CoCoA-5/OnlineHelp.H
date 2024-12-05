#ifndef CoCoA5_OnlineHelp_H
#define CoCoA5_OnlineHelp_H

//   Copyright (c)  2010  John Abbott,  Anna M. Bigatti

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

#include <string>
//using std::string;
#include <vector>
//using std::vector;


namespace CoCoA
{
namespace OnlineHelp
{

  void PrintMan(std::ostream &out, std::string keyword);
  void ReloadMan(std::ostream &out);
  void ReloadMan(std::ostream &out, const std::vector<std::string>& FileNames);
  void PrintAllExamples(std::ostream &out);
  void PrintAllExamplesWithoutOutput(std::ostream &out);
  void PrintWordlist(std::ostream &out);
  std::vector<std::string> WordList();
  const std::string& CoCoAManFileName();

} // namespace OnlineHelp
} // namespace CoCoA

// Close include guard:
#endif
