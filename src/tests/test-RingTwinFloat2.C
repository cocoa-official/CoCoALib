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


#include "CoCoA/BigRat.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/PPMonoidEv.H"
#include "CoCoA/PolyRing.H"
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

// This test checks that  (x-1/3)  divides  (3*x)^255-1  using 10-bit twin floats.
// Attempting the same computation using 8-bit precision fails (due to insufficient precision).
namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    const ring RR = NewRingTwinFloat(10);
    const PolyRing P = NewPolyRing(RR, symbols("x"));
    const RingElem x = indet(P, 0);

    // Create two polys f and g; f divides g in exact arithmetic
    const RingElem f = x - BigRat(1,3);
    const RingElem g = power(3*x, 255) - 1;

    vector<RingElem> gens;
    gens.push_back(f);
    gens.push_back(g);
    ideal I = ideal(gens);

    vector<RingElem> GBasis = TidyGens(I);
    CoCoA_ASSERT_ALWAYS(GBasis.size() == 1);
    CoCoA_ASSERT_ALWAYS(GBasis[0] == f);
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
