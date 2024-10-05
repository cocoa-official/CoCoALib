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

    //  vector<symbol> syms = symbols("s,t");
    //  vector<symbol> syms = symbols("s,t,u");
    vector<symbol> syms = symbols("s,t,u,w");
    SparsePolyRing P = NewPolyRing(NewZZmod(32003), syms);
    vector<RingElem> L;

    std::cout << "*************************************************" << std::endl;
    // vector<RingElem> LA;
    // LA = CoCoAVector(RingElem(P, "s^7-s*t^3-t"),
    //                  RingElem(P, "s*t^3-s"),
    //                  RingElem(P, "s^13-t^2"));
    // // LA = CoCoAVector(RingElem(P, "u^2+s"),
    // // 		   RingElem(P, "s^2-s*t-t*u"),
    // // 		   RingElem(P, "s^3-t^2+u"),
    // // 		   RingElem(P, "s*t^2-s*u")
    // // 		   );
  

    std::cout << "*************************************************" << std::endl;
    L.clear();
    //   std::cout << " ES_9 = ";
    //   L.push_back(RingElem(P, "s^2-t^2+w"));
    //   L.push_back(RingElem(P, "s^2-u-w"));
    //   L.push_back(RingElem(P, "s^2-t*u"));
    //   L.push_back(RingElem(P, "t^2-s*u"));
    //   L.push_back(RingElem(P, "s^3+t-u-w"));

    // std::cout << " ES_8 = ";
    // L.push_back(RingElem(P, "s^3-t^2+u"));
    // L.push_back(RingElem(P, "s^2-t-w"));
    // L.push_back(RingElem(P, "s^2-t*u"));
    // L.push_back(RingElem(P, "s*t^2-s*u"));
    // L.push_back(RingElem(P, "s+t-u^2-w"));

    std::cout << "Es homog = ";
    L.push_back(RingElem(P, "s^2-s*t-t*u"));
    L.push_back(RingElem(P, "s*t^4-s*u^4 "));
    L.push_back(RingElem(P, "s^3-t^3+u^3+s*t*u+s*t^2"));
    L.push_back(RingElem(P, "s^2+u^2-t*u"));
    L.push_back(RingElem(P, "s^2+w^2-t^2"));

    std::cout << L << std::endl;
    //  test(L, 1, "elim");
    //  test(L, 3, "IDWCLPP2");
    //   test(L, 1, "elim");
    test(L, 0, "elim");
    //   test(L, 1, "elim");
    //   test(L, 1, "elim1");
    //   test(L, 1, "IDWCLPP2");




    vector<symbol> syms2 = symbols("s,t,u,w");
    syms2.push_back(symbol("x"));
    syms2.push_back(symbol("y"));
    syms2.push_back(symbol("z"));
    syms2.push_back(symbol("a"));
    syms2.push_back(symbol("b"));
    SparsePolyRing R = NewPolyRing(NewRingFp(32003), syms2);
    vector<RingElem> LL;
    //   LL.push_back(RingElem(R, "-x+ s^2-t^2+w"));
    //   LL.push_back(RingElem(R, "-y+ s^2-u-w"));
    //   LL.push_back(RingElem(R, "-z+ s^2-t*u"));
    //   LL.push_back(RingElem(R, "-a+ t^2-s*u"));
    //   LL.push_back(RingElem(R, "-b+ s^3+t-u-w"));
    //   std::cout << " ES_9 = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x^2+ s^2-t-u+w  "));
    //   LL.push_back(RingElem(R, "-y^2+ t^2-u-w    "));
    //   LL.push_back(RingElem(R, "-z^2+ s-t*u      "));
    //   LL.push_back(RingElem(R, "-a^2+ u^2-s*w    "));
    //   LL.push_back(RingElem(R, "-b^2+ s^2+t-u-w^2"));
    //   std::cout << " ES_10 H = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x+ s^2-t-u+w  "));
    //   LL.push_back(RingElem(R, "-y+ t^2-u-w    "));
    //   LL.push_back(RingElem(R, "-z+ s-t*u      "));
    //   LL.push_back(RingElem(R, "-a+ u^2-s*w    "));
    //   LL.push_back(RingElem(R, "-b+ s^2+t-u-w^2"));
    //   std::cout << " ES_10 = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x^3+ s^3-t^2+u"));
    //   LL.push_back(RingElem(R, "-y^2+ s^2-t-w"));
    //   LL.push_back(RingElem(R, "-z^2+ s^2-t*u"));
    //   LL.push_back(RingElem(R, "-a^3+ s*t^2-s*u"));
    //   LL.push_back(RingElem(R, "-b^2+ s+t-u^2-w"));
    //   std::cout << " ES_8 = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x^15+ t^15-3*t^2-t+1"));
    //   LL.push_back(RingElem(R, "-y^23+ t^23+t^11+t^3-t-2"));
    //   std::cout << " ES_0 H = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x+ t^15-3*t^2-t+1"));
    //   LL.push_back(RingElem(R, "-y+ t^23+t^11+t^3-t-2"));
    //   std::cout << " ES_0 = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x+ s*t^5-s*t^3-t"));
    //   LL.push_back(RingElem(R, "-y+ s^3-s*t-t^2-1"));
    //   LL.push_back(RingElem(R, "-z+ s^2*t^2-s"));
    //   std::cout << " ES_1 = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x^6+ s*t^5-s*t^3-t"));
    //   LL.push_back(RingElem(R, "-y^3+ s^3-s*t-t^2-1"));
    //   LL.push_back(RingElem(R, "-z^2+ s^2*t^2-s"));
    //   std::cout << " ES_1 H = " << LL << std::endl;

    //   LL.push_back(RingElem(R, "-x^7+ s^7-s*t^3-t"));
    //   LL.push_back(RingElem(R, "-y^4+ s*t^3-s"));
    //   LL.push_back(RingElem(R, "-z^13+ s^13-t^2"));
    //   std::cout << " ES_2 = " << LL << std::endl;


    //   double T0 = CpuTime();
    //   RingElem ff = ComputeHypersurface(LL, LPP(RingElem(R, "s*t*u*w")));
    //   std::cout << "\t time: " << CpuTime()-T0;
    //   std::cout << "\t ("<< NumTerms(ff) << ")" << std::endl;
    //   std::cout << "\t ("<< ff << ")" << std::endl;

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
