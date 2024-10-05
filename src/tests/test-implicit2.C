//   Copyright (c)  2014  Anna M Bigatti

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
#include "CoCoA/time.H"

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;


//----------------------------------------------------------------------
// Test for implicitization
// functions: ImplicitDirect, elim, SliceCore
// environments: ZZ/(32003)[x,y,z,w], ZZ/(32003)[s,t,u]
//----------------------------------------------------------------------

namespace CoCoA
{

  void test(const vector<RingElem>& L, long sl, std::string algo)
  {
    double T0 = CpuTime();
    RingElem f;
    std::cout << sl << "sl + "<< algo << std::flush;
    try { f = SliceCore(L, sl, algo); }
    catch (const CoCoA::ErrorInfo& err) {cerr<<"***ERROR***"; ANNOUNCE(cerr,err);}
    std::cout << "\t time: " << CpuTime()-T0;
    std::cout << "\t ("<< NumTerms(f) << ")" << std::endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    vector<symbol> syms = symbols("s,t,u");
  
    SparsePolyRing P = NewPolyRing(NewZZmod(32003), syms);
    vector<RingElem> T;
    T.push_back(indet(P, 0));
    T.push_back(indet(P, 1));
    T.push_back(indet(P, 2));

    std::cout << "*************************************************" << std::endl;
    vector<RingElem> L6;
    L6 = CoCoAVector(RingElem(P, "s+u^2"),
                     RingElem(P, "s^2-s*t-t*u"),
                     RingElem(P, "s^3-t^2+u"),
                     RingElem(P, "s*t^2-s*u")
                     );
    std::cout << " L6 = " << L6 << std::endl;
    {
      vector<RingElem> L = L6;

      test(L, 0, "IDWCLPP2");
      //    test(L, 1, "IDWCLPP");
      //    test(L, 1, "IDWCLPP2");
      //    test(L, 2, "IDWCLPP2");
      //    test(L, 2, "elim");
    }

    std::cout << "*************************************************" << std::endl;
    vector<RingElem> L7;
    L7 = CoCoAVector(RingElem(P, "s*t^2-s*u"),
                     RingElem(P, "s^3-t^2+u^2+s+u"),
                     RingElem(P, "s^3-t^2+u"),
                     RingElem(P, "s^5-t*u"));
    std::cout << " L7 = " << L7 << std::endl;
    {
      vector<RingElem> L = L7;

      test(L, 0, "IDWCLPP2");
      test(L, 1, "IDWCLPP2");
      test(L, 2, "IDWCLPP2");
      //    test(L, 2, "elim");
    }

  }

} // end of namespace CoCoA


//**********************************************************************
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
