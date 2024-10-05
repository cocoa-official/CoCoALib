//   Copyright (c)  2017  John Abbott,  Anna M. Bigatti
//   Author: Alice Moallemy, Anna M. Bigatti

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


//#include "CoCoA/verbose.H"
#include "CoCoA/BuildInfo.H"
#include "CoCoA/GlobalManager.H"
#include "CoCoA/QuotientRing.H" // for NewZZmod
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/SparsePolyOps-ideal.H"
#include "CoCoA/RingElemInput.H"
#include "CoCoA/RingQQ.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test comes Alice Moallemy translation of CoCoA-5 functions
// for computing radical/IsRadical/IsPrimary of a 0dim ideal
// (based on the computation of MinPoly)
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    //    SetVerbosityLevel(80);
    
    ring QQ = RingQQ();
    SparsePolyRing P2 = NewPolyRing(QQ, symbols("x, y"));
    SparsePolyRing P3 = NewPolyRing(QQ, symbols("x, y, z"));
    SparsePolyRing P4 = NewPolyRing(QQ, symbols("x, y, z, t"));

    //ideals in QQ[x,y]:
    //    ideal J1 = ideal(RingElem(P2, "y^2 +2"), RingElem(P2, "x^3 + x*y^2"));
    ideal J1 = ideal(RingElems(P2, "y^2 +2, x^3 + x*y^2"));
    ideal J2 = ideal(RingElems(P2, "x^3 +x^2 +1, y^3 + y"));
    ideal J3 = ideal(RingElems(P2, "x^3 +y^3 +1, y^3 + y^2"));
    ideal J4 = ideal(RingElems(P2, "x^3 +x +1, y^3 + y^2"));
    ideal J5 = ideal(RingElems(P2, "y^3 +x^2 +1, x*y^2 + y^2"));
    ideal J6 = ideal(RingElems(P2, "x^2 +y^2 +1, x*y^2 + y^2"));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(J1));
    CoCoA_ASSERT_ALWAYS(radical(J1) == ideal(RingElem(P2, "y^2 +2"), RingElem(P2, "x^3 +x*y^2")));
    CoCoA_ASSERT_ALWAYS(IsRadical_tmp(J1));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(J2));
    CoCoA_ASSERT_ALWAYS(radical(J2) == ideal(RingElem(P2, "x^3 +x^2 +1"), RingElem(P2, "y^3 +y")));
    CoCoA_ASSERT_ALWAYS(IsRadical_tmp(J2));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(J3));
    CoCoA_ASSERT_ALWAYS(radical(J3) == ideal(RingElem(P2, "y^3 +y^2 +x*y^2"), RingElem(P2, "x^3 -y^2 +1"), RingElem(P2, "y^2 +y")));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(J3));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(J4));
    CoCoA_ASSERT_ALWAYS(radical(J4) == ideal(RingElem(P2, "y^3 +y^2"), RingElem(P2, "x^3 +x +1"), RingElem(P2, "y^2 +y")));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(J4));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(J5));
    CoCoA_ASSERT_ALWAYS(radical(J5) == ideal(RingElem(P2, "y^3 +x^2 +1"), RingElem(P2, "x*y^2 +y^2"), RingElem(P2, "x^3 +x^2 +x +1"), RingElem(P2, "-x^2*y +y")));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(J5));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(J6));
    CoCoA_ASSERT_ALWAYS(radical(J6) == ideal(RingElem(P2, "x^2 +y^2 +1"), RingElem(P2, "x*y^2 +y^2"), RingElem(P2, "y^4 +2*y^2"), RingElem(P2, "y^3 +2*y")));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(J6));

    //ideals in QQ[x, y, z]:

    ideal A = ideal(RingElem(P3,"x"), RingElem(P3, "z^2"), RingElem(P3, "y"));
    ideal C = ideal(RingElem(P3,"y-1"), RingElem(P3, "z^2-z"), RingElem(P3, "x*z-x"), RingElem(P3, "x^2-x"));

    CoCoA_ASSERT_ALWAYS(IsPrimary(A));
    CoCoA_ASSERT_ALWAYS(radical(A) == ideal(RingElem(P3, "y"), RingElem(P3, "x"), RingElem(P3, "z^2"), RingElem(P3, "z")));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(A));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(C));
    CoCoA_ASSERT_ALWAYS(radical(C) == ideal(RingElems(P3, "y -1, z^2 -z, x*z -x, x^2 -x")));
    CoCoA_ASSERT_ALWAYS(IsRadical_tmp(C));

    //ideals in QQ[x, y, z, t]:
    
    ideal K2 = ideal(RingElem(P4,"x^2*t +z^2*t +1"), RingElem(P4,"x^3 +t^3"), RingElem(P4,"x*y^2+x"), RingElem(P4,"y*z^2 +y*z"));

    CoCoA_ASSERT_ALWAYS(!IsPrimary(K2));
    CoCoA_ASSERT_ALWAYS(radical(K2) == ideal(RingElem(P4, "x^2*t +z^2*t +1"), RingElem(P4, "x^3 +t^3"), RingElem(P4, "x*y^2 +x"), RingElem(P4, "y*z^2 +y*z")));
    CoCoA_ASSERT_ALWAYS(IsRadical_tmp(K2));

    //Examples from paper: 
    //positive characteristic:

    ring ZZMod101 = NewZZmod(101);
    ring ZZMod23 = NewZZmod(23);
    SparsePolyRing P8 = NewPolyRing(ZZMod101, symbols("x, y, z"));

    ideal L6 = ideal(RingElem(P8, "x"), RingElem(P8, "y"), RingElem(P8, "z"));

    CoCoA_ASSERT_ALWAYS(IsPrimary(L6));
    CoCoA_ASSERT_ALWAYS(radical(L6) == L6);
    CoCoA_ASSERT_ALWAYS(IsRadical_tmp(L6));
    
    //Characteristic is 0:

    ideal H2 = ideal(RingElem(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"), RingElem(P4, "y^3-x"), RingElem(P4, "z^3 +z -t"), RingElem(P4, "t^3-324*z^2 +94*y^2 +76*x"));

    std::vector<RingElem> vec;
    vec.push_back(RingElem(P4, "t^3 +94*y^2 -324*z^2 +76*x"));
    vec.push_back(RingElem(P4, "z^3 +z -t"));
    vec.push_back(RingElem(P4, "y^3 -x"));
    vec.push_back(RingElem(P4, "x^4 +83*x^3 +73*y^2 -85*z^2 -437*t"));
    vec.push_back(RingElem(P4, "139623900223990351168240621*x^3*y^2*z^2*t^2 +11588783718591199146963971543*x^2*y^2*z^2*t^2 +139623900223990351168240621*x^3*y^2*t^2 -45238143672572873778509961204*x^3*y^2*z +11588783718591199146963971543*x^2*y^2*t^2 -3754765924823548523616326779932*x^2*y^2*z +10192544716351295635281565333*y*z^2*t^2 +4637188974239167542999607504652*y^2*z^2 +901970395446977668546834411660*y^2*z +5735470573401075645288988229438*y*z^2 +10192544716351295635281565333*y*t^2 +4637188974239167542999607504652*y^2 -2186789525308136879996984606102*y*z +5735470573401075645288988229438*y"));
    
    CoCoA_ASSERT_ALWAYS(!IsPrimary(H2));
    CoCoA_ASSERT_ALWAYS(radical(H2) == ideal(vec));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(H2));

    //examples John

    ideal W1 = ideal(RingElem(P3,"x^2 +x*y +z^2"), RingElem(P3, "y^2 -y*z -x"), RingElem(P3, "y^2 +z^2 -x"));

    std::vector<RingElem> vec1;
    vec1.push_back(RingElem(P3, "y*z +z^2"));
    vec1.push_back(RingElem(P3, "x*z +(-1/2)*x +(-1/2)*y"));
    vec1.push_back(RingElem(P3, "y^2 +z^2 -x"));
    vec1.push_back(RingElem(P3, "x*y -2*z^2 +(3/2)*x +(1/2)*y"));
    vec1.push_back(RingElem(P3, "x^2 +3*z^2 +(-3/2)*x +(-1/2)*y"));
    vec1.push_back(RingElem(P3, "z^3 +(-1/4)*x +(-1/4)*y"));
    vec1.push_back(RingElem(P3, "(-1/2)*z^2 +(1/4)*x +(1/4)*y +(1/4)*z"));
    CoCoA_ASSERT_ALWAYS(!IsPrimary(W1));
    CoCoA_ASSERT_ALWAYS(radical(W1) == ideal(vec1));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(W1));

    ideal W2 = ideal(RingElem(P3,"y^2 +x*z +y*z"), RingElem(P3, "x^2 -y^2 +x*z"), RingElem(P3, "z^2 +x +y"));

    std::vector<RingElem> vec2;
    vec2.push_back(RingElem(P3, "z^2 +x +y"));
    vec2.push_back(RingElem(P3, "x*z -y*z -x"));
    vec2.push_back(RingElem(P3, "y^2 +2*y*z +x"));
    vec2.push_back(RingElem(P3, "x*y +3*y*z -y"));
    vec2.push_back(RingElem(P3, "x^2 +3*y*z +2*x"));
    vec2.push_back(RingElem(P3, "y*z -2*x -3*y -z"));
    CoCoA_ASSERT_ALWAYS(!IsPrimary(W2));
    CoCoA_ASSERT_ALWAYS(radical(W2) == ideal(vec2));
    CoCoA_ASSERT_ALWAYS(!IsRadical_tmp(W2));

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
