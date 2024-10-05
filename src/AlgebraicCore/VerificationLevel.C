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


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/VerificationLevel.C,v 1.7 2022/02/18 14:12:01 abbott Exp $
// $Log: VerificationLevel.C,v $
// Revision 1.7  2022/02/18 14:12:01  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.6  2021/02/10 19:40:00  abbott
// Summary: Added noexcept (sometimes instead of throw()) -- redmine 1572
//
// Revision 1.5  2020/06/17 15:49:29  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.4  2020/02/11 16:12:20  abbott
// Summary: Added some checks for bad ostream (even to mem fns myOutput); see redmine 969
//
// Revision 1.3  2018/03/15 10:47:01  abbott
// Summary: Added new fns, IsGuaranteed & level (for VerificationLevel)
//
// Revision 1.2  2018/03/14 14:30:38  abbott
// Summary: Use new error OutOfRange
//
// Revision 1.1  2018/03/13 17:34:07  abbott
// Summary: Added new files VerificationLevel
//
// Revision 1.6  2017/07/21 13:20:18  abbott
// Summary: Copyright template is now Abbott+Bigatti
//
// Revision 1.5  2010/12/17 16:10:11  abbott
// Changed copyright year.
//
// Revision 1.4  2009/06/22 15:15:56  abbott
// Changed copyright year
//
// Revision 1.3  2008/03/12 14:41:28  abbott
// Updated copyright year
//
// Revision 1.2  2007/10/30 17:14:07  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.3  2007/03/03 15:24:39  cocoa
// Changed 2006 to 2007.
//
// Revision 1.2  2006/10/06 14:04:14  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
