//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/QBGenerator.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;


namespace CoCoA
{

  // A simple test of QBGenerator; result is checked by comparing printed output.

  void program()
  {
    GlobalManager CoCoAFoundations;

    ring QQxy = NewPolyRing(RingQQ(), SymbolRange("x", 1, 2));
    QBGenerator QBG(PPM(QQxy));
    const PPMonoidElem PP1 = QBG.myCorners().front();
    QBG.myCornerPPIntoQB(PP1);

    int counter = 0;
    while (!QBG.myCorners().empty()) 
    {
      ++counter;
      const PPMonoidElem t = QBG.myCorners().front();
      // Use a simple (pointless?) criterion for deciding whether to "avoid" or put in QB:
      if (counter%7 == 0)
        QBG.myCornerPPIntoAvoidSet(t);
      else
        QBG.myCornerPPIntoQB(t);
    }
    cout << "QBG = " << QBG << endl;
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
