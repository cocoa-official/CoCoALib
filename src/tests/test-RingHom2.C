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
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingHom.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ring.H"


#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for RingHom
// Qxy->Qxy
//----------------------------------------------------------------------
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    const PolyRing Qxy = NewPolyRing(RingQQ(), symbols("x,y"));

    const vector<RingElem>& x = indets(Qxy);
    vector<RingElem> imx;
    imx.push_back(x[0]*x[0]);
    imx.push_back(x[1]*x[1]);
    RingHom Qxy2Qxy = PolyRingHom(Qxy, Qxy, CoeffEmbeddingHom(Qxy), imx);
    RingHom Qxy2QxyBis = PolyAlgebraHom(Qxy, Qxy, imx);

    cout << "Simple test on polynomial ring hom from " << Qxy << "  to  " << Qxy << endl
         << "sending " << indet(Qxy,0) << "  to  " << imx[0] << "  and" << endl
         << "sending " << indet(Qxy,1) << "  to  " << imx[1] << endl;

    RingElem f = x[0]+2*x[1];

    CoCoA_ASSERT_ALWAYS(Qxy2QxyBis(f) == Qxy2Qxy(f));

    cout << f << " in " << owner(f) << "  maps to " << Qxy2Qxy(f) << endl;
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
