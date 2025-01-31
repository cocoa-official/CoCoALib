#ifndef CoCoA_GlobalIO_H
#define CoCoA_GlobalIO_H

//   Copyright (c)  2014  John Abbott,  Anna Bigatti

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
// using std::istream;
// using std::ostream;

namespace CoCoA
{

  std::istream& GlobalInput();
  std::ostream& GlobalOutput();
  std::ostream& GlobalErrput();
  std::ostream& GlobalLogput();

  std::istream& SetGlobalInput(std::istream& in);
  std::ostream& SetGlobalOutput(std::ostream& out);
  std::ostream& SetGlobalErrput(std::ostream& err);
  std::ostream& SetGlobalLogput(std::ostream& log);


  //  void InputFailCheck(const std::istream& in, const char* const FnName);
  void InputFailCheck(const std::istream& in, const std::string& FnName);

}

#endif
