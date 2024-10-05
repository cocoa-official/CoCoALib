//   Copyright (c) 2020  John Abbott,  Anna M. Bigatti
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


////////////////////////////////////////////////////////////
// Checks that COCOA5_VER_MAJ, MIN & MINMIN are alphanumeric
// Exit with 0 if so, otherwise exit with 1
// Makefile compiles and runs this prog!!!
////////////////////////////////////////////////////////////

#include <string>

// CPP trick to convert a macro value into a string
#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

// true is str is non-empty and alphanumeric
bool IsAlphanumeric(const std::string& str)
{
  if (str.empty()) return false;
  for (char ch: str)
    if (!isalnum(ch)) return false;
  return true;
}

int main()
{
  using std::string;
  if (!IsAlphanumeric(string(STRINGIFY(COCOA5_VER_MAJ)))) return 1;
  if (!IsAlphanumeric(string(STRINGIFY(COCOA5_VER_MIN)))) return 1;
  if (!IsAlphanumeric(string(STRINGIFY(COCOA5_VER_MINMIN)))) return 1;
  return 0;
}
