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
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingTwinFloat.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ring.H"
#include "CoCoA/symbol.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


// This test computes a largish G-basis using twin-float coeffs; the computation
// will probably throw InsufficientPrecision if BuffBits is too small.
namespace CoCoA
{

  void solve(long BuffBits)
  {
    const ring RR = NewRingTwinFloat(20, BuffBits, 48);

    const PolyRing P = NewPolyRing(RR, symbols("x,y,z"), StdDegLex(3));
    const RingElem f(P, "x^41 - z^40*(x-z)");
    const RingElem g(P,"x^6-9*x^5*z+x*z^5-11*z^6-y^6");

    vector<RingElem> gens;
    gens.push_back(f);
    gens.push_back(g);
    const ideal I = ideal(gens);

    const vector<RingElem> GBasis = TidyGens(I); // might throw an exception
    CoCoA_ASSERT_ALWAYS(BuffBits >= 143);
    CoCoA_ASSERT_ALWAYS(GBasis.size() == 64);
    CoCoA_ASSERT_ALWAYS(GBasis[0] == g);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    for (long BuffBits = 32; BuffBits <= 160; BuffBits += 32)
    {
      try
      {
        solve(BuffBits);
      }
      catch (const RingTwinFloat::InsufficientPrecision&)
      {
        CoCoA_ASSERT_ALWAYS(BuffBits < 143);
      }
    }
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
