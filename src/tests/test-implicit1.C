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
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/TmpImplicit.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/VectorOps.H"
#include "CoCoA/time.H"

#include <iostream>
using std::cerr;
using std::endl;

#include <vector>
using std::vector;


//----------------------------------------------------------------------
// Test for implicitization
// functions: ImplicitDirect, elim, SliceCore
// environments: ZZ/(32003)[x,y,z], ZZ/(32003)[s,t]
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
    
    //   T0 = CpuTime();
    //   RingElem g;
    //   std::cout << sl << "sl + "<< algo << std::flush;
    //   try { g = SliceCoreQQ(L, sl, algo); }
    //   catch (const CoCoA::ErrorInfo& err) {cerr<<"***ERROR***"; ANNOUNCE(cerr,err);}
    //   std::cout << "\t time: " << CpuTime()-T0;
    //   std::cout << "\t ("<< NumTerms(g) << ")" << std::endl;
  }

  void program()
  {
    GlobalManager CoCoAFoundations;

    const vector<symbol> syms = symbols("s,t");
  
    SparsePolyRing P = NewPolyRing(NewZZmod(32003), syms);
    //  SparsePolyRing P = NewPolyRing(RingQQ(), syms);
    vector<RingElem> T;
    T.push_back(indet(P, 0));
    T.push_back(indet(P, 1));

    std::cout << "*************************************************" << std::endl;
    vector<RingElem> L1;
    L1 = RingElems(P, "s*t^5-s*t^3-t,"
                      "s^3-s*t-t^2-1,"
                      "s^2*t^2-s");

    std::cout << " L1 = " << L1 << std::endl;
    {
      vector<RingElem> L = L1;

      test(L, 0, "IDWCLPP2");
      //test(L, 1, "elimth");
      //test(L, 1, "elim");
      //test(L, 1, "IDWCLPP2");
      //    test(L, 2, "IDWCLPP2");
      
    }

    std::cout << "*************************************************" << std::endl;
    vector<RingElem> L2;
    L2 = RingElems(P, "s^7-s*t^3-t,"
                       "s*t^3-s,"
                       "s^13-t^2");
    std::cout << " L2 = " << L2 << std::endl;
    {
      vector<RingElem> L = L2;

      test(L, 0, "IDWCLPP2");
    //    test(L, 1, "IDWCLPP2");
    //    test(L, 2, "IDWCLPP2");

    //    test(L, 1, "elim");
    //    test(L, 1, "elim1");
    //    test(L, 1, "elimt");
    //    test(L, 2, "elim");
    }

//   std::cout << "*************************************************" << std::endl;
//   vector<RingElem> L200;
//   L200 = RingElems(P, "s^7-t,"
//                       "s*t-s,"
//                       "s^10-t^2");
//   std::cout << " L200 = " << L200 << std::endl;
//   {
//     vector<RingElem> L = L200;

//     test(L, 0, "IDWCLPP2");
//     //    test(L, 1, "IDWCLPP2");
//     //    test(L, 2, "IDWCLPP2");

//     test(L, 1, "elim");
//     test(L, 1, "elimt");
//     //    test(L, 2, "elim");
//   }
  
//   std::cout << "*************************************************" << std::endl;
//   vector<RingElem> L4;
//   L4 = RingElem(P, "s^4-s*t^3,"
//                    "s^37-t^3-s,"
//                    "s*t^3-s");
//   std::cout << " L4 = " << L4 << std::endl;
//   {
//     vector<RingElem> L = L4;
//     RingElem f;
//     double T0=0;
    
//   }
  
//   std::cout << "*************************************************" << std::endl;
//   vector<RingElem> L5;
//   L5 = RingElems(P,"s^6+s^3+1,"
//                    "t^6+t^3+1,"
//                    "s^7+3*s^4*t^3+3*s^3*t^4+t^7");
//   std::cout << " L5 = " << L5 << std::endl;
//   {
//     vector<RingElem> L = L5;
//     RingElem f;
//     double T0=0;

//   }
  
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
