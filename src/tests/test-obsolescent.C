//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidOv.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"

// *** EXTRA INCLUDE DIRECTIVE TO MAKE OBSOLESCENT FNS VISIBLE ***
#include "CoCoA/obsolescent.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test checks that the GlobalManager option "AllowObsolescentFns"
// actually does enable an obsolescent fn to be called.
// It makes 3 attempts: one should be successful, the other two not.
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    // First we make a "successful" call (because we told GlobalManager to allow it)
    {
      GlobalManager CoCoAFoundations(AllowObsolescentFns);
      PPMonoid PPM = NewPPMonoidOv(symbols("x,y"), lex);
      PPMonoidElem t = indet(PPM, 0);
      CoCoA_ASSERT_ALWAYS(IsRadical(t)); // the call to IsRadical will print a logging message.
      // other obsolescent functions
      ring P = NewPolyRing(RingQQ(), symbols("x,y"));
      CoCoA_ASSERT_ALWAYS(AreGensSquareFreeMonomial(ideal(indet(P,0))));
    }

    // Now we make an "unsuccessful" call (because by default GlobalManager forbids calling obsolescent fns)
    try
    {
      GlobalManager CoCoAFoundations;
      PPMonoid PPM = NewPPMonoidOv(symbols("x,y"), lex);
      PPMonoidElem t = indet(PPM, 0);
      CoCoA_ASSERT_ALWAYS(IsRadical(t)); // this should throw ERR::OBSOLESCENT
      CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "NEVER GET HERE");
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      CoCoA_ASSERT_ALWAYS(err == ERR::OBSOLESCENT);
    }

    // Now we make a second "unsuccessful" call (because by default GlobalManager forbids calling obsolescent fns)
    try
    {
      GlobalManager CoCoAFoundations(ForbidObsolescentFns);
      PPMonoid PPM = NewPPMonoidOv(symbols("x,y"), lex);
      PPMonoidElem t = indet(PPM, 0);
      CoCoA_ASSERT_ALWAYS(IsRadical(t)); // this should throw ERR::OBSOLESCENT
      CoCoA_THROW_ERROR(ERR::ShouldNeverGetHere, "NEVER GET HERE");
    }
    catch (const CoCoA::ErrorInfo& err)
    {
      CoCoA_ASSERT_ALWAYS(err == ERR::OBSOLESCENT);
    }
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
    cerr << "***ERROR***  UNCAUGHT CoCoA Error";
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
