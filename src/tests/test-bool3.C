//   Copyright (c)  2010  Anna Bigatti

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
#include "CoCoA/error.H"
#include "CoCoA/bool3.H"

#include <iostream>
using std::cerr;
using std::endl;

namespace CoCoA
{

  class SillyClass
  {
  public:
    SillyClass();

    bool3 myValue1;
    bool3 myValue2;
    bool3 myValue3;
    bool3 myValue4;
  };

  
  SillyClass::SillyClass():
      myValue1(true3),
      myValue2(uncertain3),
      myValue3(false3)
  { myValue4 = uncertain3; }


  void program()
  {
    GlobalManager CoCoAFoundations;

    SillyClass s;

    CoCoA_ASSERT_ALWAYS(IsTrue3(s.myValue1));
    CoCoA_ASSERT_ALWAYS(IsUncertain3(s.myValue2));
    CoCoA_ASSERT_ALWAYS(IsFalse3(s.myValue3));
    CoCoA_ASSERT_ALWAYS(IsUncertain3(s.myValue4));  
  }

} // end of namespace CoCoA


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
