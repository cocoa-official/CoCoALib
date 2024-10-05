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

#include "GlobalIO.H"
#include "CoCoA/error.H"

#include <iostream>
// using std::cin;
// using std::cout;
// using std::clog;
// using std::cerr;

namespace CoCoA
{

  namespace // anonymous namespace
  {
    std::istream* InPtr = &std::cin;
    std::ostream* OutPtr = &std::cout;
    std::ostream* ErrPtr = &std::cerr;
    std::ostream* LogPtr = &std::clog;
  } // end of anonymous namespace


  std::istream& GlobalInput()
  {
    return *InPtr;
  }

  std::ostream& GlobalOutput()
  {
    return *OutPtr;
  }

  std::ostream& GlobalErrput()
  {
    return *ErrPtr;
  }

  std::ostream& GlobalLogput()
  {
    return *LogPtr;
  }


  std::istream& SetGlobalInput(std::istream& in)
  {
    std::istream& prev = GlobalInput();
    InPtr = &in;
    return prev;
  }

  std::ostream& SetGlobalOutput(std::ostream& out)
  {
    std::ostream& prev = GlobalOutput();
    OutPtr = &out;
    return prev;
  }

  std::ostream& SetGlobalErrput(std::ostream& err)
  {
    std::ostream& prev = GlobalErrput();
    ErrPtr = &err;
    return prev;
  }

  std::ostream& SetGlobalLogput(std::ostream& log)
  {
    std::ostream& prev = GlobalLogput();
    LogPtr = &log;
    return prev;
  }


  void InputFailCheck(const std::istream& in, const std::string& FnName)
  {
    if (!in) CoCoA_THROW_ERROR(ERR::InputFail, FnName);
  }

}

// OLD DOCUMENTATION (was in io.txt)
// These files supply four standard global I/O channels for the CoCoA
// library.  These channels mimic the global channels present in the C++
// STL: ``cin``, ``cout``, ``cerr``, and ``clog``.

// The current choices for these four channels can be obtained by calling
// these functions (inside namespace CoCoA):
// ```
//  std::istream& GlobalInput();
//  std::ostream& GlobalOutput();
//  std::ostream& GlobalErrput();
//  std::ostream& GlobalLogput();
// ```

// By default the standard global I/O channels for the CoCoA library are
// the corresponding ones of the standard C++ library.  Alternative
// choices may be specified by calling these functions (see *NOTE*)
// ```
//  std::istream& SetGlobalInput(std::istream& in);
//  std::ostream& SetGlobalOutput(std::ostream& out);
//  std::ostream& SetGlobalErrput(std::ostream& err);
//  std::ostream& SetGlobalLogput(std::ostream& log);
// ```
// In each case the value returned is the previous istream/ostream
// associated with that channel.

// *NOTE* the procedures for changing the settings of the global i/o
// streams maintain a reference to the value supplied.  If you use a
// local variable as argument, make sure that its value is not destroyed
// before you have finished i/o on it.  See the example program ex-io.C.

// As recommended in Meyers's book I have put the globals in an anonymous
// namespace.  I chose to use plain pointers for the global variables
// ``InPtr``, ``OutPtr``, ``ErrPtr`` and ``LogPtr``; references are unsuitable because
// they cannot be reseated in C++, and ``auto_ptr`` is unsuitable because
// we do not want to own the streams.  I do not believe that there can be
// problems with race-conditions when these four global variables are
// initialized since we use only the addresses of ``std::cin``, etc.
// However, there could be race-condition problems with subsequent
// changes to the values.


// It is tempting to make these simple functions into inline functions,
// but inlining is tricky with the global variables in an anonymous namespace.
// And anyway the minor gain in performance will easily be swamped by the
// high costs of actually conducting I/O; so making them inline would
// really be quite pointless.

// It would be nice to avoid the potential pitfalls of dangling references,
// but I do not currently see how to achieve this.



// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/server/GlobalIO.C,v 1.4 2022/02/22 20:40:04 abbott Exp $
// $Log: GlobalIO.C,v $
// Revision 1.4  2022/02/22 20:40:04  abbott
// Summary: Updated copyright message (redmine 1555)
//
// Revision 1.3  2020/06/17 15:49:31  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.2  2019/03/19 12:53:00  bigatti
// -- just some debugging info
//
// Revision 1.1  2014/05/15 12:30:52  abbott
// Summary: New files for impls previously in CoCoA/io.HC but needed only in the server
// Author: JAA
//
//
