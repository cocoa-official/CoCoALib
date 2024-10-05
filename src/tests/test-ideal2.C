//   Copyright (c)  2007,2020  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingDistrMPolyClean.H"
#include "CoCoA/RingDistrMPolyInlFpPP.H"
#include "CoCoA/RingDistrMPolyInlPP.H"
#include "CoCoA/RingFp.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/ideal.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;

//----------------------------------------------------------------------
// Test for ideals in particular cases
// functions: IsElem, IsPrime
// environments: RingZZ,
//----------------------------------------------------------------------
namespace CoCoA
{

  void TestZZ()
  {
    ring ZZ = RingZZ();
  
    ideal Z = ideal(zero(ZZ));
    ideal J3 = ideal(RingElem(ZZ,3));
    ideal J35 = ideal(RingElem(ZZ,3), RingElem(ZZ,5));
    ideal J46 = ideal(RingElem(ZZ,4), RingElem(ZZ,6));
    ideal J6 = ideal(RingElem(ZZ,6));

    CoCoA_ASSERT_ALWAYS( !IsElem(RingElem(ZZ, 7), Z) );
    CoCoA_ASSERT_ALWAYS( IsPrime(Z) );
    CoCoA_ASSERT_ALWAYS( IsPrime(J3) );
    CoCoA_ASSERT_ALWAYS( IsPrime(J46) );
    CoCoA_ASSERT_ALWAYS( IsPrime(intersect(J35, J46)) );
    CoCoA_ASSERT_ALWAYS( !IsPrime(intersect(J3, J46)) );
    CoCoA_ASSERT_ALWAYS( IsPrime(intersect(J46, Z)) );
    CoCoA_ASSERT_ALWAYS( intersect(J3, J6) == J6 );
  }


  void TestIsPrimary_PID() // in QQ[x] which is PID
  {
    PolyRing P = NewPolyRing(RingQQ(), symbols("x")); // PID
    RingElem f1 = RingElem(P, "x");
    RingElem f2 = RingElem(P, "x+1");
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
      {
        if (i==0 && j==0) continue;
        ideal I(power(f1,i)*power(f2,j));
        CoCoA_ASSERT_ALWAYS(IsPrimary(I) == (i==0 || j==0));
      }
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    TestZZ();
    TestIsPrimary_PID();
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
