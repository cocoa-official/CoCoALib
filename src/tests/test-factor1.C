//   Copyright (c)  2006,2017  John Abbott,  Anna M. Bigatti

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



//----------------------------------------------------------------------
// WARNING: the user interface to the factorizer is not final!
// Of course, this is only a minimal test.
//----------------------------------------------------------------------

#include "CoCoA/BuildInfo.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/factor.H"


#include <vector>
using std::vector;
#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  void TestQQ()
  {
    const ring QQ = RingQQ();
    SparsePolyRing P = NewPolyRing(QQ, symbols("x,y"));

    const vector<RingElem>& x = indets(P);
    const RingElem f = 5*(power(2*x[0],96) - power(3*x[1],96));
    const factorization<RingElem> FacInfo = factor(f);
    const vector<RingElem>& facs = FacInfo.myFactors();   // handy alias
    const vector<long>& mult = FacInfo.myMultiplicities();// handy alias
    const int NumFacs = len(facs);
    CoCoA_ASSERT_ALWAYS(NumFacs == 12);
    RingElem prod = FacInfo.myRemainingFactor();
    for (int i = 0; i < NumFacs; ++i)
    {
      CoCoA_ASSERT_ALWAYS(mult[i] == 1);
      prod *= facs[i];
    }
    CoCoA_ASSERT_ALWAYS(prod == f);
  }

  void TestZZ()
  {
    const ring ZZ = RingZZ();
    SparsePolyRing P = NewPolyRing(ZZ, symbols("x,y"));

    const vector<RingElem>& x = indets(P);
    const RingElem f = 5*(power(2*x[0],96) - power(3*x[1],96));
    const factorization<RingElem> FacInfo = factor(f);
    const vector<RingElem>& facs = FacInfo.myFactors();   // handy alias
    const vector<long>& mult = FacInfo.myMultiplicities();// handy alias
    const int NumFacs = len(facs);
    CoCoA_ASSERT_ALWAYS(NumFacs == 12);
    RingElem prod = FacInfo.myRemainingFactor();
    for (int i = 0; i < NumFacs; ++i)
    {
      CoCoA_ASSERT_ALWAYS(mult[i] == 1);
      prod *= facs[i];
    }
    CoCoA_ASSERT_ALWAYS(prod == f);
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    TestQQ();
    TestZZ();
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
