//   Copyright (c)  2008,2014  John Abbott,  Anna M. Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.

//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.

//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.

//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/BuildInfo.H"
#include "CoCoA/error.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <limits>
using std::numeric_limits;
#include <vector>
using std::vector;


#define NEVER_GET_HERE CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "Execution should never get here!")

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    try { symbol("_x"); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadSymbolHead); }

    const long MaxLong = numeric_limits<long>::max();
    const long MinLong = numeric_limits<long>::min();

    symbol a("a", -1, -2);
    try { subscript(a,2); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadSymbolSubscript); }
    try { subscript(a,MinLong); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadSymbolSubscript); }
    try { subscript(a,MaxLong); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadSymbolSubscript); }
    try { SymbolRange("a",MinLong,MaxLong); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ArgTooBig); }

    vector<long> inds(4);
    symbol sym1("x", inds); // x[0,0,0,0]
    {
      // Try largest single range -- should give an error.
      inds[0] = MinLong;
      symbol A("x", inds);
      inds[0] = MaxLong;
      symbol Z("x", inds);
      try { SymbolRange(A, Z); NEVER_GET_HERE; }
      catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ArgTooBig); }
      // Now create a range where upper subscript is MaxLong.
      inds[0] = MaxLong-1;
      A = symbol("x", inds);
      CoCoA_ASSERT_ALWAYS(SymbolRange(A, Z).size() == 2);
    }
    if (sizeof(long) == 4)
    {
      // Next 2 lines try to make a vector of 2^32+6 symbols -- will fail on 32-bit systems
      inds[0] = 40521; inds[1] = 105990;
      try { SymbolRange(sym1, symbol("x",inds)); NEVER_GET_HERE; }
      catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ArgTooBig); }
    }
    // The next lines try to make a vector of 2^64+4 symbols -- should fail on 32 & 64 bit systems.
    inds[0] = 55809; inds[1] = 43404;  inds[2] = 49476; inds[3] = 384772;
    symbol sym2("x", inds);
    try { SymbolRange(sym1, sym2); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::ArgTooBig); }
    try { SymbolRange(sym2, sym1); NEVER_GET_HERE; }
    catch (const ErrorInfo& err) { CoCoA_ASSERT_ALWAYS(err == ERR::BadSymbolRange); }
  }

} // end of namespace CoCoA


//----------------------------------------------------------------------
// Use main() to handle any uncaught exceptions and warn the user about them.
int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}

//----------------------------------------------------------------------
// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/tests/test-symbol1.C,v 1.13 2022/02/11 09:54:00 abbott Exp $
// $Log: test-symbol1.C,v $
// Revision 1.13  2022/02/11 09:54:00  abbott
// Summary: Updated copyright (redmine 855)
//
// Revision 1.12  2021/02/10 19:46:45  abbott
// Summary: Revised after updating symbol (redmine 925); removed tests involving large ulong
//
// Revision 1.11  2020/06/17 15:49:31  abbott
// Summary: Changed CoCoA_ERROR into CoCoA_THROW_ERROR
//
// Revision 1.10  2017/09/06 11:57:50  abbott
// Summary: Changed ERR::SERIOUS into ERR::ShouldNeverGetHere
//
// Revision 1.9  2016/09/16 16:36:42  abbott
// Summary: Changed TEST_ASSERT into CoCoA_ASSERT_ALWAYS; removed all include assert.H lines
//
// Revision 1.8  2015/05/11 13:00:53  abbott
// Summary: Put all code into namespace CoCoA (redmine #642)
// Author: JAA
//
// Revision 1.7  2014/10/27 15:37:55  abbott
// Summary: Experimental change: put program inside namespace CoCoa
// Author: JAA
//
// Revision 1.6  2014/06/14 19:29:44  abbott
// Summary: Removed test that "a1" was forbidden; reflowed source code to avoid long lines
// Author: JAA
//
// Revision 1.5  2012/05/24 14:49:22  bigatti
// -- changed symbol "index" into "subscripts"
//
// Revision 1.4  2011/08/23 06:40:30  bigatti
// -- fixed with new name "BigInt" for old "ZZ"
//
// Revision 1.3  2010/12/17 16:06:25  abbott
// Ensured that all i/o is on standard C++ streams (instead of GlobalInput, etc)
//
// Revision 1.2  2008/12/12 11:32:01  abbott
// Updated Makefiles to make the new test/example for symbol visible.
//
// Revision 1.1  2008/12/12 11:29:47  abbott
// Fixed a bug in SymbolRange.  Added example and test for symbols.
//
// Revision 1.2  2007/10/30 17:14:06  abbott
// Changed licence from GPL-2 only to GPL-3 or later.
// New version for such an important change.
//
// Revision 1.1.1.1  2007/03/09 15:16:12  abbott
// Imported files
//
// Revision 1.2  2007/03/07 11:34:44  bigatti
// -- minimized #include's
//
// Revision 1.1  2007/03/06 12:34:54  bigatti
// -- added test for GMPAllocator (problems before/after GlobalManager)
//
// Revision 1.4  2007/03/03 14:15:45  bigatti
// -- "foundations" renamed into "GlobalManager"
//
// Revision 1.3  2007/02/12 15:29:07  bigatti
// -- added strings ShortDescription and LongDescription for indexing
//
// Revision 1.2  2007/02/10 18:44:04  cocoa
// Added "const" twice to each test and example.
// Eliminated dependency on io.H in several files.
// Improved BuildInfo, and added an example about how to use it.
// Some other minor cleaning.
//
// Revision 1.1.1.1  2006/05/30 11:39:36  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/10/17 10:46:53  cocoa
// Imported files
//
// Revision 1.1.1.1  2005/05/03 15:47:30  cocoa
// Imported files
//
// Revision 1.1  2005/04/29 15:42:02  cocoa
// Improved documentation for GMPAllocator.
// Added example program for GMPAllocator.
// Added example program for simple ops on polynomials.
// Added two new ctors for (principal) ideals (from long, and from BigInt).
// Added (crude) printing for PPMonoids.
// Updated library.H (#included GMPAllocator.H).
//
// Revision 1.3  2005/04/27 16:14:56  cocoa
// Cleaned up example programs -- added "free use" permit.
// Changed a couple of ErrorInfo object names, and added
// ERR::NotTrueGCDDomain.
//
// Revision 1.2  2005/04/21 15:12:19  cocoa
// Revised NewPolyRing as Dag Arneson suggested (perhaps just an interim
// measure).
// Brought example programs up to date (new name for CoCoA error
// information objects).
//
// Revision 1.1.1.1  2005/01/27 15:12:13  cocoa
// Imported files
//
// Revision 1.3  2004/12/09 15:08:42  cocoa
// -- added log info
//
