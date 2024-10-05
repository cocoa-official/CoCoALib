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

// RCS header/log in the next few lines.
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/assert.H,v 1.9 2022/02/18 14:11:51 abbott Exp $
// $Log: assert.H,v $
// Revision 1.9  2022/02/18 14:11:51  abbott
// Summary: Updated copyright notice (now restrictive; see redmine 1555)
//
// Revision 1.8  2021/11/03 19:20:59  abbott
// Summary: Added comment
//
// Revision 1.7  2021/02/10 19:32:01  abbott
// Summary: Added noexcept (redmine 1572)
//
// Revision 1.6  2021/01/07 14:48:33  abbott
// Summary: Corrected copyright
//
// Revision 1.5  2019/03/18 17:56:49  abbott
// Summary: Changed 0 into nullptr where appropriate
//
// Revision 1.4  2017/03/13 12:17:09  abbott
// Summary: Include PREPROCESSOR_DEFNS.H (in header files which need the info)
//
// Revision 1.3  2016/09/16 16:27:09  abbott
// Summary: Created CoCoA_ASSERT_ALWAYS (was TEST_ASSERT in each test file)
//
// Revision 1.2  2007/10/30 17:14:11  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:11  abbott
// Imported files
//
// Revision 1.2  2006/12/06 17:17:44  cocoa
// -- removed #include "config.H"
//
// Revision 1.1  2006/10/06 14:04:57  cocoa
// The new assert header and implementation files.
// A new test.
//
//

#endif
