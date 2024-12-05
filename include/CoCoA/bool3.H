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

#endif
