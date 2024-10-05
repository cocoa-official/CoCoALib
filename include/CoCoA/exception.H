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



// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/exception.H,v 1.12 2022/10/25 17:20:38 abbott Exp $
// $Log: exception.H,v $
// Revision 1.12  2022/10/25 17:20:38  abbott
// Summary: Added attribute [[noreturn]]
//
// Revision 1.11  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.10  2021/11/03 19:21:41  abbott
// Summary: Changed ctor arg types (was const string&, now string)
//
// Revision 1.9  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.8  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.7  2020/06/19 19:42:59  abbott
// Summary: Cleaned; removed cruft
//
// Revision 1.6  2020/06/19 19:39:21  abbott
// Summary: Now all throws go through new template fn ThrowException; seems to work
//
// Revision 1.5  2020/06/19 14:57:34  abbott
// Summary: Added new fn ThrowException (similar to ThrowError)
//
// Revision 1.4  2020/02/11 17:16:17  abbott
// Summary: All mem fns like myOutput now return void (redmine 969)
//
// Revision 1.3  2017/07/22 12:59:37  abbott
// Summary: Added virt mem fn myOutputSelf; added PrintInFrame.
//
// Revision 1.2  2017/07/15 15:16:33  abbott
// Summary: Removed useless include directive
//
// Revision 1.1  2015/06/26 14:56:25  abbott
// Summary: Created new class "exception"
// Author: JAA
//
//

#endif
