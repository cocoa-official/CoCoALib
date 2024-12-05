//   Copyright (c)  2018  John Abbott,  Anna M. Bigatti

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


#include "CoCoA/VerificationLevel.H"
#include "CoCoA/error.H"

#include <ostream>
//using std::ostream;

namespace CoCoA
{

  VerificationLevel::VerificationLevel(long vl): // NOT noexcept (arg check)
      myLevel(vl)
  {
    if (vl < 0 || vl > 1000)
      CoCoA_THROW_ERROR(ERR::OutOfRange, "VerificationLevel ctor");
  }


  VerificationLevel guaranteed() noexcept
  {
    VerificationLevel ans(0);
    ans.myLevel = -1; // not normally allowed, but we have private access
    return ans;
  }
  

  std::ostream& operator<<(std::ostream& out, const VerificationLevel& vl)
  {
    if (!out) return out;  // short-cut for bad ostreams

    if (level(vl) < 0) return out << "guaranteed";
    return out << "VerificationLevel(" << level(vl) << ")";
  }

} // end of namespace CoCoA
