//   Copyright (c)  2010 Anna Bigatti

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
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/FractionField.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/PolyRing.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/error.H"
#include "CoCoA/symbol.H"

#include <iostream>
using std::cerr;
using std::endl;
#include<sstream>
using std::ostringstream;


namespace CoCoA
{

  void program()
  {
    // This test is for testing printing functions
    GlobalManager CoCoAFoundations;


    PolyRing R = NewPolyRing(RingQQ(), symbols("a,b")); // QQ[a,b]
    ring K = NewFractionField(R); // QQ(a,b)
    PolyRing Kxy = NewPolyRing(K, symbols("x,y"));
    RingHom f = CanonicalHom(R,K);
    RingHom g = CanonicalHom(K,Kxy);
    RingHom h = g(f);
    RingElem a = RingElem(R, symbol("a"));
    RingElem b = RingElem(R, symbol("b"));
    RingElem Inv_a = 1/f(a);
    RingElem x = RingElem(Kxy, symbol("x"));
    RingElem Inv_a_x = g(Inv_a)*x;
  
    //  std::cout << 2*a/3-1 << std::endl;
    //  std::cout << 2*h(a)*x-1 << std::endl;
  
    { ostringstream s;  s << a/2;          CoCoA_ASSERT_ALWAYS(s.str() == "(1/2)*a"); }
    { ostringstream s;  s << -a/3;         CoCoA_ASSERT_ALWAYS(s.str() == "(-1/3)*a"); }
    { ostringstream s;  s << -x/3;         CoCoA_ASSERT_ALWAYS(s.str() == "(-1/3)*x"); }
    { ostringstream s;  s << h(-a)*x*x + h(1-a)*x + h(1-a);
      CoCoA_ASSERT_ALWAYS(s.str() == "-a*x^2 +(-a +1)*x -a +1"); }
    { ostringstream s;  s << x - h(a*b);   CoCoA_ASSERT_ALWAYS(s.str() == "x -a*b"); }
    { ostringstream s;  s << -x;           CoCoA_ASSERT_ALWAYS(s.str() == "-x"); }
    { ostringstream s;  s << -x+1;         CoCoA_ASSERT_ALWAYS(s.str() == "-x +1"); }
    { ostringstream s;  s << x-1;          CoCoA_ASSERT_ALWAYS(s.str() == "x -1"); }
    { ostringstream s;  s << x- h(a);      CoCoA_ASSERT_ALWAYS(s.str() == "x -a"); }
    { ostringstream s;  s << 1/f(a);       CoCoA_ASSERT_ALWAYS(s.str() == "1/a"); }
    { ostringstream s;  s << 1/h(a);       CoCoA_ASSERT_ALWAYS(s.str() == "1/a"); }
    { ostringstream s;  s << g(1/f(a))*x;  CoCoA_ASSERT_ALWAYS(s.str() == "(1/a)*x"); }
    { ostringstream s;  s << h(a)*x;       CoCoA_ASSERT_ALWAYS(s.str() == "a*x"); }
    { ostringstream s;  s << x/(-1);       CoCoA_ASSERT_ALWAYS(s.str() == "-x"); }
    { ostringstream s;  s << x-one(Kxy)/2; CoCoA_ASSERT_ALWAYS(s.str() == "x -1/2"); }
    { ostringstream s;  s << x-h(a/2);     CoCoA_ASSERT_ALWAYS(s.str() == "x -a/2"); }
    { ostringstream s;  s << x-h(a)/2;     CoCoA_ASSERT_ALWAYS(s.str() == "x -a/2"); }
    { ostringstream s;  s << x-h(a+1)/2;   CoCoA_ASSERT_ALWAYS(s.str() == "x +(-a -1)/2"); }
    { R->myOutputSelfLong(std::cout); }
    { std::cout << std::endl; }
    { NewPolyRing(R, SymbolRange("x",1,3))->myOutputSelfLong(std::cout); }
    { std::cout << std::endl; }
    { NewPolyRing(NewFractionField(R), SymbolRange("x",1,3))->myOutputSelfLong(std::cout); }
    { std::cout << std::endl; }
    { NewPolyRing(NewQuotientRing(R,"a^2-2"), SymbolRange("x",1,3))->myOutputSelfLong(std::cout); }
    { std::cout << std::endl; }
  
    //  std::cout << x-h(a)/2 << std::endl;

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
