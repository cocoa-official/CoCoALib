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
