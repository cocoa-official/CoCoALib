#ifndef CoCoA_bool3_H
#define CoCoA_bool3_H

//   Copyright (c)  2003-2010,2012,2021  John Abbott,  Anna M. Bigatti

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


#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class OpenMathOutput; // forward declaration -- defined in OpenMath.H

  /*-----------------------------------------------------------------*/
  /** three way "booleans": false3, uncertain3, and true3.
   **/
  /*-----------------------------------------------------------------*/


  class bool3
  {
  public:
    constexpr bool3() noexcept; // default value is "uncertain3"
    explicit constexpr bool3(bool b) noexcept;
    // Default dtor, copy ctor, and assignment are fine.
    bool3& operator=(bool b) noexcept;

    friend bool IsFalse3(bool3 flag) noexcept;
    friend bool IsUncertain3(bool3 flag) noexcept;
    friend bool IsTrue3(bool3 flag) noexcept;
  private: // impl details
    enum TruthValueSet { false3, uncertain3, true3 }; // not enum class: I want the names to be injected!
    constexpr TruthValueSet FromStdBool(bool b) noexcept { return (b?(true3):(false3)); }
  private: // data member
    TruthValueSet myTruthValue;
  };


  std::ostream& operator<<(std::ostream& out, bool3 flag);
  OpenMathOutput& operator<<(OpenMathOutput& out, bool3 flag);

  //---------------------------------------------------------------------------
  // inline functions

  inline constexpr bool3::bool3() noexcept:
      myTruthValue(uncertain3)
  {}

  inline constexpr bool3::bool3(bool b) noexcept:
      myTruthValue(FromStdBool(b))
  {}


  constexpr bool3 true3(true);
  constexpr bool3 uncertain3;
  constexpr bool3 false3(false);


  inline bool3& bool3::operator=(bool b) noexcept
  { myTruthValue = FromStdBool(b); return *this; }


  inline bool IsFalse3(bool3 flag) noexcept
  { return flag.myTruthValue == bool3::false3; }


  inline bool IsUncertain3(bool3 flag) noexcept
  { return flag.myTruthValue == bool3::uncertain3; }


  inline bool IsTrue3(bool3 flag) noexcept
  { return flag.myTruthValue == bool3::true3; }


} // end of namespace CoCoA



// RCS header/log in the next few lines.
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/bool3.H,v 1.13 2022/02/18 14:11:51 abbott Exp $
// $Log: bool3.H,v $
// Revision 1.13  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.12  2021/01/31 10:09:12  abbott
// Summary: Added noexcept
//
// Revision 1.11  2021/01/28 21:18:53  abbott
// Summary: Major revision: redesign, and uses constexpr
//
// Revision 1.10  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.9  2018/04/05 10:26:12  abbott
// Summary: Corrected comment
//
// Revision 1.8  2018/03/29 09:42:53  bigatti
// -- just initial comment
//
// Revision 1.7  2013/06/03 17:45:07  bigatti
// -- fixed comment
//
// Revision 1.6  2012/05/29 07:45:23  abbott
// Implemented simplification change to bool3:
//  changed names of the constants,
//  changed names of the testing fns.
//
// Revision 1.5  2010/01/21 15:23:36  abbott
// Added two comments.
//
// Revision 1.4  2010/01/21 14:09:09  bigatti
// -- added copy constructor
//
// Revision 1.3  2009/05/14 09:41:25  abbott
// Added explicit ctor for a bool3 from a normal bool.
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.6  2007/03/08 18:42:05  cocoa
// Cleaned up whitespace.
//
// Revision 1.5  2006/12/06 17:17:44  cocoa
// -- removed #include "config.H"
//
// Revision 1.4  2006/11/24 17:22:05  cocoa
// -- removed OpenMathFwd.H
//
// Revision 1.3  2006/11/02 13:25:44  cocoa
// Simplification of header files: the OpenMath classes have been renamed.
// Many minor consequential changes.
//
// Revision 1.2  2006/10/06 14:04:15  cocoa
// Corrected position of #ifndef in header files.
// Separated CoCoA_ASSERT into assert.H from config.H;
// many minor consequential changes (have to #include assert.H).
// A little tidying of #include directives (esp. in Max's code).
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.3  2006/05/12 16:10:58  cocoa
// Added OpenMathFwd.H, and tidied OpenMath.H.
// Many consequential but trivial changes.
//
// Revision 1.2  2006/03/12 21:28:34  cocoa
// Major check in after many changes
//
// Revision 1.1.1.1  2005/10/17 10:46:54  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.2  2005/04/19 14:06:04  cocoa
// Added GPL and GFDL licence stuff.
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.6  2004/11/12 15:49:29  cocoa
// Tidying prior to 0.90 release.
// (a) documentation improved (or marked as poor)
// (b) sundry minor improvements to the code
//
// Revision 1.4  2004/06/29 15:34:29  cocoa
// Tidied the functions dealing with bool3 values.
//
// Revision 1.3  2004/05/27 16:14:02  cocoa
// Minor revision for new coding conventions.
//
// Revision 1.2  2003/10/09 12:16:38  cocoa
// New coding convention for rings.
//
// Revision 1.3  2003/06/23 16:17:57  abbott
// Minor cleaning prior to public release.
// Slightly changed order in a few places.
//
// Revision 1.2  2003/05/14 17:11:03  abbott
// Added missing keyword "inline" -- omitting it was a silly mistake.
//
// Revision 1.1  2003/04/24 15:56:45  abbott
// Initial revision
//

#endif
