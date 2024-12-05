#ifndef CoCoA_exception_H
#define CoCoA_exception_H

//   Copyright (c)  2015,2021  John Abbott,  Anna M. Bigatti

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
// using std::string;
#include <iosfwd>
//using std::ostream;

namespace CoCoA
{

  class exception
  {
  public:
    explicit exception(std::string message): myMessage(message), myContext() {};
    exception(std::string message, std::string context): myMessage(message), myContext(context) {};
///???    const char* what() const noexcept { return myMessage.c_str(); }
    virtual void myOutputSelf(std::ostream& out) const;
    virtual ~exception();
    friend const std::string& message(const exception& e) noexcept;
    friend const std::string& context(const exception& e) noexcept;
  protected: // data members
    std::string myMessage;
    std::string myContext;
  };


  std::ostream& operator<<(std::ostream& out, const exception& exc);
  void PrintInFrame(std::ostream& out, const exception& exc);


  inline void JustBeforeThrowing() noexcept
  { /*Does nothing; just enables gdb intercept*/ }


  template<class EXCEPTION>
  [[noreturn]]
  void ThrowException(const EXCEPTION& exc)
  {
    JustBeforeThrowing();
    // An attempt to have verbose output here failed:
    // e.g. too many reports of InsuffPrec which were then caught & properly handled.
    throw exc;
  }



  // -- inlines:
  inline const std::string& message(const exception& e) noexcept
  {
    return e.myMessage;
  }

  inline const std::string& context(const exception& e) noexcept
  {
    return e.myContext;
  }


} // end of namespace CoCoA

#endif
