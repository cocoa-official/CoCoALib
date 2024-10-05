#ifndef CoCoA_VerificationLevel_H
#define CoCoA_VerificationLevel_H

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

#include "CoCoA/config.H"

#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class VerificationLevel
  {
  public:
    explicit VerificationLevel(long vl); // NOT noexcept (arg check)
  private: // data members
    int myLevel;
    friend VerificationLevel guaranteed() noexcept;
    friend long level(VerificationLevel vl) noexcept;
  };

  VerificationLevel guaranteed() noexcept;
  std::ostream& operator<<(std::ostream& out, const VerificationLevel& VL);

  inline long level(VerificationLevel vl) noexcept
  { return vl.myLevel; }
  
  bool inline IsGuaranteed(VerificationLevel vl) noexcept
  { return (level(vl) < 0); }

} // end of namespace CoCoA



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/VerificationLevel.H,v 1.4 2022/02/18 14:11:51 abbott Exp $
// $Log: VerificationLevel.H,v $
// Revision 1.4  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.3  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.2  2018/03/15 10:47:01  abbott
// Summary: Added new fns, IsGuaranteed & level (for VerificationLevel)
//
// Revision 1.1  2018/03/13 17:34:07  abbott
// Summary: Added new files VerificationLevel
//
//

#endif
