//   Copyright (c)  2016  John Abbott,  Anna M. Bigatti

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
#include "CoCoA/QuotientRing.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/error.H"
#include "CoCoA/SparsePolyRing.H"


#include <iostream>
using std::cerr;
using std::endl;
#include <string>
using std::string;

//----------------------------------------------------------------------
// This test checks whether IsSqFree works correctly for polynomials.
//----------------------------------------------------------------------

namespace CoCoA
{
  // Put your code inside namespace CoCoA to avoid possibile
  // ambiguities with STL fns sharing names with CoCoALib fns.

  void SquareFree(ConstRefRingElem f)
  {
    CoCoA_ASSERT_ALWAYS(IsSqFree(f) == true);
  }

  void NotSquareFree(ConstRefRingElem f)
  {
    CoCoA_ASSERT_ALWAYS(IsSqFree(f) == false);
  }

  void program()
  {
    // you may use CoCoALib functions AFTER creating a GlobalManager:
    GlobalManager CoCoAFoundations;

    ring QQx = NewPolyRing(RingQQ(), symbols("x"));
    ring QQxyz = NewPolyRing(RingQQ(), symbols("x,y,z"));
    ring Fp = NewZZmod(3);
    ring Fpx = NewPolyRing(Fp, symbols("x"));
    ring Fpxyz = NewPolyRing(Fp, symbols("x,y,z"));

    // in QQ[x]
    SquareFree(RingElem(QQx, "1"));
    SquareFree(RingElem(QQx, "x"));
    SquareFree(RingElem(QQx, "x^2-1"));
    SquareFree(RingElem(QQx, "x^3-x"));

    NotSquareFree(RingElem(QQx, "x^3"));
    NotSquareFree(RingElem(QQx, "x^6-x^4"));
    NotSquareFree(RingElem(QQx, "(x^2+x+1)^2"));
    NotSquareFree(RingElem(QQx, "(x^5+x+1)^3*(x^2+x+1)^7"));

    // in QQ[x,y,z]
    SquareFree(RingElem(QQxyz, "1"));
    SquareFree(RingElem(QQxyz, "x"));
    SquareFree(RingElem(QQxyz, "y"));
    SquareFree(RingElem(QQxyz, "z"));
    SquareFree(RingElem(QQxyz, "x*y*z"));
    SquareFree(RingElem(QQxyz, "(x^2+x+1)*(y^3+y+1)*(z^4+z+1)"));
    SquareFree(RingElem(QQxyz, "(x+y)*(y+z)*(z+x)"));

    NotSquareFree(RingElem(QQxyz, "x^3"));
    NotSquareFree(RingElem(QQxyz, "y^4"));
    NotSquareFree(RingElem(QQxyz, "z^5"));
    NotSquareFree(RingElem(QQxyz, "x^3*y*z"));
    NotSquareFree(RingElem(QQxyz, "x*y^4*z"));
    NotSquareFree(RingElem(QQxyz, "x*y*z^5"));
    NotSquareFree(RingElem(QQxyz, "(x^2+x+1)^4*(y^3+y+1)*(z^4+z+1)"));
    NotSquareFree(RingElem(QQxyz, "(x^2+x+1)*(y^3+y+1)^3*(z^4+z+1)"));
    NotSquareFree(RingElem(QQxyz, "(x^2+x+1)*(y^3+y+1)*(z^4+z+1)^2"));
    NotSquareFree(RingElem(QQxyz, "(x+y)^5*(y+z)*(z+x)"));
    NotSquareFree(RingElem(QQxyz, "(x+y)*(y+z)^5*(z+x)"));
    NotSquareFree(RingElem(QQxyz, "(x+y)*(y+z)*(z+x)^5"));

    // in Fp[x]
    SquareFree(RingElem(Fpx, "1"));
    SquareFree(RingElem(Fpx, "x"));
    SquareFree(RingElem(Fpx, "x^2-1"));
    SquareFree(RingElem(Fpx, "x^3-x"));

    NotSquareFree(RingElem(Fpx, "x^3"));
    NotSquareFree(RingElem(Fpx, "x^3-1"));
    NotSquareFree(RingElem(Fpx, "x^27-1"));
//    NotSquareFree(RingElem(Fpx, "(x^2+1)*(x^27-x)"));

    // in Fp[x,y,z]
    SquareFree(RingElem(Fpxyz, "1"));
    SquareFree(RingElem(Fpxyz, "x"));
    SquareFree(RingElem(Fpxyz, "y"));
    SquareFree(RingElem(Fpxyz, "z"));
    SquareFree(RingElem(Fpxyz, "x*y*z"));
    SquareFree(RingElem(Fpxyz, "(x^3+x+1)*(y^3+y+1)*(z^3+z+1)"));
    SquareFree(RingElem(Fpxyz, "(x+y)*(y+z)*(z+x)"));

    NotSquareFree(RingElem(Fpxyz, "x^3"));
    NotSquareFree(RingElem(Fpxyz, "y^4"));
    NotSquareFree(RingElem(Fpxyz, "z^5"));
    NotSquareFree(RingElem(Fpxyz, "x^3*y*z"));
    NotSquareFree(RingElem(Fpxyz, "x*y^4*z"));
    NotSquareFree(RingElem(Fpxyz, "x*y*z^5"));
    NotSquareFree(RingElem(Fpxyz, "(x^2+x+1)^4*(y^3+y+1)*(z^4+z+1)"));
    NotSquareFree(RingElem(Fpxyz, "(x^2+x+1)*(y^3+y+1)^3*(z^4+z+1)"));
    NotSquareFree(RingElem(Fpxyz, "(x^2+x+1)*(y^3+y+1)*(z^4+z+1)^2"));
    NotSquareFree(RingElem(Fpxyz, "(x+y)^5*(y+z)*(z+x)"));
    NotSquareFree(RingElem(Fpxyz, "(x+y)*(y+z)^5*(z+x)"));
    NotSquareFree(RingElem(Fpxyz, "(x+y)*(y+z)*(z+x)^5"));
    NotSquareFree(RingElem(Fpxyz, "(x^9-x)*(y^9-y)*(z^9-1)"));
    NotSquareFree(RingElem(Fpxyz, "(x^9-x)*(y^9-1)*(z^9-z)"));
    NotSquareFree(RingElem(Fpxyz, "(x^9-1)*(y^9-y)*(z^9-z)"));
    NotSquareFree(RingElem(Fpxyz, "(x^9-x)*(y^9-1)*(z^9-1)"));
    NotSquareFree(RingElem(Fpxyz, "(x^9-1)*(y^9-y)*(z^9-1)"));
    NotSquareFree(RingElem(Fpxyz, "(x^9-1)*(y^9-1)*(z^9-z)"));
    NotSquareFree(RingElem(Fpxyz, "((x+y)^9-(x+y))*((y+z)^9-1)*((z+x)^9-1)"));
    NotSquareFree(RingElem(Fpxyz, "((x+y)^9-1)*((y+z)^9-(y+z))*((z+x)^9-1)"));
    NotSquareFree(RingElem(Fpxyz, "((x+y)^9-1)*((y+z)^9-1)*((z+x)^9-(z+x))"));
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
