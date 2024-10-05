//   Copyright (c)  2013  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/factor.H"
#include "CoCoA/symbol.H"
#include "CoCoA/error.H"

#include <vector>
using std::vector;
#include <iostream>
using std::cerr;
using std::endl;


namespace CoCoA
{

  RingElem product(const factorization<RingElem>& FacInfo)
  {
    RingElem ans = FacInfo.myRemainingFactor();
    const vector<RingElem>& facs = FacInfo.myFactors();   // handy alias
    const vector<long>& mult = FacInfo.myMultiplicities();// handy alias
    const long NumFacs = len(facs);
    for (long i=0; i < NumFacs; ++i)
    {
      ans *= power(facs[i], mult[i]);
    }
    return ans;
  }

  void TestCharP()
  {
    ring Fp = NewZZmod(5);
    PolyRing P = NewPolyRing(Fp, symbols("x,y,z"));

    // test (annoying) trivial case
    const RingElem f = one(P);
    CoCoA_ASSERT_ALWAYS(f == product(SqFreeFactor(f)));
    CoCoA_ASSERT_ALWAYS(2*f == product(SqFreeFactor(2*f)));

    const RingElem x = indet(P,0);
    const RingElem y = indet(P,1);
    const RingElem z = indet(P,2);

    const RingElem f1 = 3*x+2*y+1;
    const RingElem f2 = 4*x+3*z+3;
    const RingElem f3 = 2*y+4*z+3;
    const RingElem f4 = 3*x+2*y+4*z+2;

    vector<int> tbl; tbl.push_back(0); tbl.push_back(2); tbl.push_back(5); tbl.push_back(7); tbl.push_back(26);
    for (int e1=0; e1 < 4; ++e1) // skip last (o/w too slow)
      for (int e2=0; e2 < 5; ++e2)
        for (int e3=0; e3 < 4; ++e3) // skip last (o/w too slow)
          for (int e4=0; e4 < 3; ++e4) // skip last two (o/w too slow)
          {
            if (e1 == 3 && e3 == 3) continue; // skip some more slow cases
            const RingElem f = power(f1,tbl[e1])*power(f2,tbl[e2])*power(f3,tbl[e3])*power(f4,tbl[e4]);
            CoCoA_ASSERT_ALWAYS(f == product(SqFreeFactor(f)));
          }
  }

  void TestChar0()
  {
    ring QQ = RingQQ();
    PolyRing P = NewPolyRing(QQ, symbols("x,y,z"));

    // test (annoying) trivial case
    const RingElem f = one(P);
    CoCoA_ASSERT_ALWAYS(f ==  product(SqFreeFactor(f)));
    CoCoA_ASSERT_ALWAYS(2*f ==  product(SqFreeFactor(2*f)));

    const RingElem x = indet(P,0);
    const RingElem y = indet(P,1);
    const RingElem z = indet(P,2);

    const RingElem f1 = 3*x+2*y+1;
    const RingElem f2 = 4*x+3*z+3;
    const RingElem f3 = 2*y+4*z+3;
    const RingElem f4 = 3*x+2*y+4*z+2;

    for (int e1=0; e1 < 4; ++e1)
      for (int e2=0; e2 < 4; ++e2)
        for (int e3=0; e3 < 4; ++e3)
          for (int e4=0; e4 < 4; ++e4)
          {
            const RingElem f = power(f1,e1)*power(f2,e2)*power(f3,e3)*power(f4,e4);
            CoCoA_ASSERT_ALWAYS(f == product(SqFreeFactor(f)));
          }
  }


  void program()
  {
    GlobalManager CoCoAFoundations;

    TestCharP();
    TestChar0();
  }

}  // end of namespace CoCoA


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
