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
    if (MaxLineLen < 1) CoCoA_THROW_ERROR(ERR::BadArg, "fold: MaxLineLen must be at least 1");
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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/utils.C,v 1.12 2022/02/18 14:12:03 abbott Exp $
// $Log: utils.C,v $
// Revision 1.12  2022/02/18 14:12:03  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.11  2021/02/10 19:40:01  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.10  2021/01/07 15:23:39  abbott
// Summary: Corrected copyright
//
// Revision 1.9  2020/12/05 13:05:10  abbott
// Summary: Added IsDecimal for istream and ostream
//
// Revision 1.8  2020/06/17 15:49:30  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.7  2020/02/03 17:05:27  abbott
// Summary: Added new fn fold
//
// Revision 1.6  2010/10/22 12:29:04  abbott
// Removed some entirely useless code.
//
// Revision 1.5  2008/04/21 12:53:56  abbott
// Added a comment.
//
// Revision 1.4  2008/02/11 15:14:55  abbott
// Added useless dummy functions to keep ranlib quiet on MacOS.
//
// Revision 1.3  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.2  2007/05/21 12:53:08  abbott
// No real change in the end.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.1  2006/08/07 21:23:25  cocoa
// Removed almost all publicly visible references to SmallExponent_t;
// changed to long in all PPMonoid functions and SparsePolyRing functions.
// DivMask remains to sorted out.
//
