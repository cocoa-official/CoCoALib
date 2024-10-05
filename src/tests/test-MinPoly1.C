//   Copyright (c)  2017-2018  John Abbott and Anna M. Bigatti

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
#include "CoCoA/RingQQ.H"
#include "CoCoA/SparsePolyOps-MinPoly.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/VerificationLevel.H"

#include <iostream>
using std::cerr;
using std::endl;

//----------------------------------------------------------------------
// This test comes from redmine issue 1101: there was no check on the
// result of polynomial reconstruction in MinPolyQuot.
//----------------------------------------------------------------------

namespace CoCoA
{

  void program()
  {
    GlobalManager CoCoAFoundations;

    PolyRing Qt = RingQQt(1);
    
    PolyRing P = NewPolyRing(RingQQ(), symbols("x,y,z"));
    //-------------------------------------------------------
    ideal I1 = ideal(RingElem(P, "3*x^3*y +3*y*z^3 +1"),
                     RingElem(P, "2*x*y*z^2 +3*x"),
                     RingElem(P, "2*y^3*z +1"));
    RingElem x = indet(P,0);
    RingElem mp1 = MinPolyQuot(x, I1, x);
    CoCoA_ASSERT_ALWAYS(mp1 == RingElem(P, "x^16 +(15/2)*x^10 +(45/4)*x^4 +(1594195/5184)*x"));
    RingElem mp1t = MinPolyQuot(x, I1, indet(Qt,0));
    CoCoA_ASSERT_ALWAYS(mp1t == RingElem(Qt, "t^16 +(15/2)*t^10 +(45/4)*t^4 +(1594195/5184)*t"));
    RingElem mp1t3 = MinPolyQuot(x, I1, indet(Qt,0), VerificationLevel(3));
    CoCoA_ASSERT_ALWAYS(mp1t == mp1t3);

    //-------------------------------------------------------
    ideal I2 = ideal(RingElem(P, "3*x*y*z^3 +z^2 +1"),
                     RingElem(P, "3*y^3*z +z^2"),
                     RingElem(P, "2*x*y^2*z^2 +3*x*y*z^2"));
    RingElem mp2 = MinPolyQuot(x, I2, x);
    CoCoA_ASSERT_ALWAYS(mp2 == RingElem(P, "x^2 +(-106000/4782969)*x"));
    //-------------------------------------------------------
    ideal I3 = ideal(RingElem(P, "x^2*y^2 +3*y^3*z +1"),
                     RingElem(P, "x*y*z^2 +x*y^2"),
                     RingElem(P, "2*y^3 +3"));
    RingElem mp3 = MinPolyQuot(x, I3, x);
    CoCoA_ASSERT_ALWAYS(mp3 == RingElem(P, "x^13 +(-81/2)*x^9 +(8/9)*x^7 +(2187/4)*x^5 +36*x^3 +(-1594195/648)*x"));
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
