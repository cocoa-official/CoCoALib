//   Copyright (c)  2006  John Abbott and Anna M. Bigatti

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


#include "CoCoA/utils.H"

#include <algorithm>
using std::copy;
//#include <string>
using std::string;

namespace CoCoA
{

  // Currently there is nothing I can usefully put here.

  std::string fold(const std::string& str, long MaxLineLen)
  {
    if (MaxLineLen < 1)  CoCoA_THROW_ERROR1(ERR::ReqPolyRing);
    string ans;
    const long EndOfStr = len(str);
    long BlockStart = 0;
    while (BlockStart < EndOfStr)
    {
      if (BlockStart != 0)
        ans += '\n';
      ans.append(str, BlockStart, MaxLineLen);
      BlockStart += MaxLineLen; // BUG??? could overflow???
    }
    return ans;
  }


  bool IsDecimal(const std::istream& in) noexcept
  {
    return (in.flags() & std::ios_base::dec);
// return !((in.flags() & std::ios_base::oct) || (in.flags() && std::ios_base::hex));
  }

  bool IsDecimal(const std::ostream& out) noexcept
  {
    return (out.flags() & std::ios_base::dec);
// return !((out.flags() & std::ios_base::oct) || (out.flags() && std::ios_base::hex));
  }


} // end of namespace CoCoA
