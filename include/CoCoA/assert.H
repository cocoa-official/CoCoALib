#ifndef CoCoA_assert_H
#define CoCoA_assert_H

//   Copyright (c)  2005,2016  John Abbott,  Anna M. Bigatti

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

#include "CoCoA/PREPROCESSOR_DEFNS.H"

namespace CoCoA
{

  // This fn prints mesg on std::cerr, and then throws ERR::AssertFail.  Never returns.
  void AssertionFailed(const char* const cond, const char* const OrigFile, unsigned long OrigLine);

} // end of namespace CoCoA


// Handy macro for making assertions in test files (see src/tests/test-*.C).
#define CoCoA_ASSERT_ALWAYS(cond)   do { if (!(cond)) CoCoA::AssertionFailed(#cond, __FILE__, __LINE__); } while(false)

// The CoCoA_ASSERT macro checks conditions only if compilation is with the
// CoCoA_DEBUG flag set; if the flag is not set, the macro does nothing
// whatsoever (not even evaluating its argument).  After compiling with
// CoCoA_DEBUG set, the executable checks that the condition is true;
// if not, it calls CoCoA::AssertionFailed which prints out a useful message.
#ifndef CoCoA_DEBUG
#define CoCoA_ASSERT(cond) static_cast<void>(0)
#else
#define CoCoA_ASSERT(cond) CoCoA_ASSERT_ALWAYS(cond)
#endif

#endif
