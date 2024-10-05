//   Copyright (c)  2007  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/OpenMathXML.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/BigInt.H"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <utility>
using std::pair;
using std::make_pair;


namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations(UseNonNegResidues);

    OpenMathOutput OMOut(new OpenMathOutputXML(cout));

    const ring ZZ = RingZZ();

    RingElem n(ZZ);
    cout << "OpenMath for RingZZ is: " << endl;
    OMOut << ZZ;
    cout << endl << endl;

    cout << "OpenMath for RingElem(ZZ, 0) is " << endl;
    OMOut << n;
    cout << endl << endl;

    n = -12345;
    cout << "OpenMath for RingElem(ZZ, " << n << ") is " << endl;
    OMOut << n;
    cout << endl << endl;

    n = power(n, 100);
    cout << "OpenMath for RingElem(ZZ, " << n << ") is " << endl;
    OMOut << n;
    cout << endl << endl;

    const ring ZZ13 = NewZZmod(13);
    cout << "OpenMath for ZZ13 is: " << endl;
    OMOut << ZZ13;
    cout << endl << endl;

    RingElem r(ZZ13);
    cout << "OpenMath for RingElem(ZZ13, " << r << ") is " << endl;
    OMOut << r;
    cout << endl << endl;
    r = 12345;
    cout << "OpenMath for RingElem(ZZ13, " << r << ") is " << endl;
    OMOut << r;
    cout << endl << endl;

    ring QQ = RingQQ();
    PolyRing P = NewPolyRing(QQ, symbols("x,y,z"));
    cout << "OpenMath for " << P << " is" << endl;
    OMOut << P;
    RingElem x = indet(P, 0);
    RingElem f = (((2*x+4)*x+5)*x+4)*x+2;
    cout << "OpenMath for " << f << " is " << endl;
    OMOut << f;
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
